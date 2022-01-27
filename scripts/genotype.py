#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   genotype.py: genotypes from extracted reads.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#   This file is part of arcasHLA.
#
#   arcasHLA is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   arcasHLA is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with arcasHLA.  If not, see <https://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

import os
import sys
import re
import json
import pickle
import argparse
import logging as log

import numpy as np
import math
import pandas as pd

from datetime import date
from argparse import RawTextHelpFormatter
from textwrap import wrap
from collections import Counter, defaultdict
from itertools import combinations

from reference import check_ref
from arcas_utilities import *
from align import *

#-------------------------------------------------------------------------------

__version__     = '0.4.0'
__date__        = '2022-01-27'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir = os.path.dirname(os.path.realpath(__file__)) + '/../'

hla_json        = rootDir + 'dat/ref/hla.p.json'
hla_idx         = rootDir + 'dat/ref/hla.idx'
hla_freq        = rootDir + 'dat/info/hla_freq.tsv'
parameters_json = rootDir + 'dat/info/parameters.json'

#-----------------------------------------------------------------------------
# Genotype
#-----------------------------------------------------------------------------

def expectation_maximization(eqs, lengths, allele_idx, population, prior, 
                             tolerance, max_iterations, drop_iterations, 
                             drop_threshold):
    '''Quantifies allele transcript abundance. Based on the methods
       used in HISAT-genotype (http://dx.doi.org/10.1101/266197).
    '''

    # Divides raw counts between alleles for the first iteration 
    # of transcript quantification
    def initial_abundances(eqs, lengths, population):
    
        # Divides counts equally between alleles in the same compatibility class
        def divide_equally(alleles, count):
            n_alleles = len(alleles)
            for allele in alleles:
                counts[allele] += count / n_alleles
                
        # Divides counts proportionally to allele frequency
        def divide_prior(alleles, count, allele_prob):
            total_prob = sum(allele_prob.values())
            for allele in alleles:
                counts[allele] += count * (allele_prob[allele]/total_prob)
        
        counts = defaultdict(float)
        undivided_counts = defaultdict(float)
        for alleles, count in eqs:
            
            allele_prior = defaultdict(float)
            for idx in alleles:
                undivided_counts[idx] += count
                allele = process_allele(allele_idx[idx][0], 2)
                if population and allele in prior:
                    allele_prior[idx] = prior[allele][population]

            if population and allele_prior:
                divide_prior(alleles, count, allele_prior)
                continue
                    
            divide_equally(alleles, count)
                
        return counts_to_abundances(counts), undivided_counts
    
    # Normalizes counts by allele length and convert to abundances
    def counts_to_abundances(counts):
        abundances = defaultdict(float)
        
        for allele, count in counts.items():
            length = lengths[allele]
            abundances[allele] = count / length
        
        total_abundance = sum(abundances.values())
        
        for allele, abundance in abundances.items():
            abundances[allele] = abundance / total_abundance
            
        return abundances
    
    # Redistribute counts between alleles in the same compatibility 
    # class based on their overall abundance
    def update_abundances(eqs, abundances):
        counts = defaultdict(float)
        
        for alleles, count in eqs:
            alleles = [allele for allele in alleles if allele in abundances]
            total_abundance = sum([abundances[allele] for allele in alleles])
            
            if total_abundance == 0:
                continue
            
            for allele in alleles:
                counts[allele] += count * (abundances[allele]/total_abundance)
                
        return counts_to_abundances(counts)
    
    # Drop low support alleles after a specified number of iterations if their 
    # abundance is less than a specified proportion of the greatest abundance
    def drop_alleles(eqs, abundances, drop_iterations, drop_threshold, 
                     iterations, converged):
        if iterations == 1:
            abundances = {allele:abundance 
                          for allele, abundance in abundances.items()
                          if abundance > 0.0}
            
        elif iterations >= drop_iterations or converged:
            threshold = drop_threshold * max(abundances.values())
            abundances = {allele:abundance 
                          for allele, abundance in abundances.items() 
                          if abundance >= threshold}
        return abundances, eqs
    
    # Compute square root of sum of squares
    def SRSS(theta):
        square_sum = 0.0
        for i in theta:
            square_sum += i**2
        return math.sqrt(square_sum)

    # Check if sum difference between two iterations is below tolerance
    def check_convergence(theta0, theta_prime):
        diff = [theta_prime[allele] - theta0[allele] for allele in theta0]
        residual_error = SRSS(diff)
        return residual_error < tolerance

    converged = False
    iterations = 1

    theta0, undivided_counts = initial_abundances(eqs, lengths, population)
    
    log.info('[genotype] Top 10 alleles by undivided read count:')
    log.info('\t\t{: <20}    {: >10}\t'.format('allele', 'read count'))
    
    for idx, count in sorted(undivided_counts.items(), 
                             key=lambda x: x[1], 
                             reverse = True)[:10]:
                                     
        log.info('\t\t{: <20}    {: >10.0f}\t'
                 .format(process_allele(allele_idx[idx][0], 3), count))

    log.info(f'\n[genotype] Quantifying allele transcript abundance')
    
    # SQUAREM - accelerated EM
    # R. Varadhan & C. Roland (doi: 10.1 1 1 1/j. 1467-9469.2007.00585.X)
    # Used by HISAT-genotype, originaly used by Sailfish
    while iterations < max_iterations and not converged:
        # Get next two steps
        theta1 = update_abundances(eqs, theta0)
        theta2 = update_abundances(eqs, theta1)
        theta_prime = defaultdict(float)

        r = dict()
        v = dict()
        sum_r = 0.0
        sum_v = 0.0
        
        # Compute r and v
        for allele in theta1:
            r[allele] = theta1[allele] - theta0[allele]
            v[allele] = (theta2[allele] - theta1[allele]) - r[allele]
            
        srss_r = SRSS(r.values())
        srss_v = SRSS(v.values())

        if srss_v != 0:
            # Compute step length
            alpha = -(srss_r / srss_v)
            for allele in r:
                value =   theta0[allele] \
                        - 2*alpha*r[allele] \
                        + (alpha**2)*v[allele]
                        
                theta_prime[allele] = value
                
            step_min = min(theta_prime.values())
            step_max = max(theta_prime.values())

            # Adjust step rather than kicking out alleles with a negative result
            if step_min < 0:
                theta_prime = {allele:(value-step_min)/(step_max-step_min) 
                               for allele, value in theta_prime.items()}
                               
                total = sum(theta_prime.values())
                
                theta_prime = {allele:value/total 
                               for allele, value in theta_prime.items()}
     
            # Update abundances with given the new proportions
            theta_prime = update_abundances(eqs, theta_prime)
            
        else:
            theta_prime = theta1
            
        converged = check_convergence(theta0, theta_prime)

        theta0, eqs = drop_alleles(eqs, theta_prime, drop_iterations, 
                                   drop_threshold, iterations, converged)
        iterations += 1
         
    log.info(f'[genotype] EM converged after {iterations} iterations')
    
    return theta0

def predict_genotype(eqs, allele_idx, allele_eq, em_results, gene_count, 
                     population, prior, zygosity_threshold):
                     
    '''Predicts most likely genotype using scoring based on proportion of 
       explained reads, tie-breaking with allele priors.
    '''
    # Returns number of reads explained by an allele
    def get_count(a):
        observed_eqs = allele_eq[a]
        return sum([eqs[idx][1] for idx in observed_eqs])
    
    # Returns number of reads explained by a pair of alleles
    def get_pair_count(a1, a2):
        if type(a1) == tuple:
            a1_eqs = set.union(*[allele_eq[idx] for idx in a1])
        else:
            a1_eqs = allele_eq[a1]
        if type(a2) == tuple:
            a2_eqs =set.union(*[allele_eq[idx] for idx in a2])
        else:
            a2_eqs = allele_eq[a2]
        
        observed_eqs = a1_eqs | a2_eqs
        
        return sum([eqs[idx][1] for idx in observed_eqs])
        
    # Returns non-shared counts for a pair of alleles
    def get_nonshared_count(a1, a2):
        if type(a1) == tuple:
            a1_eqs = set.union(*[allele_eq[idx] for idx in a1])
        else:
            a1_eqs = allele_eq[a1]
        if type(a2) == tuple:
            a2_eqs =set.union(*[allele_eq[idx] for idx in a2])
        else:
            a2_eqs = allele_eq[a2]
        
        a1_nonshared_eqs = a1_eqs - a2_eqs
        a2_nonshared_eqs = a2_eqs - a1_eqs
        
        a1_count = sum([eqs[idx][1] for idx in a1_nonshared_eqs])
        a2_count = sum([eqs[idx][1] for idx in a2_nonshared_eqs])
        
        return a1_count, a2_count
        
    explained_reads = dict()
    if len(em_results) > 1:
        grouped_indices = defaultdict(set)
        for idx, alleles, abundances in em_results:
            allele = process_allele(alleles[0], 2)
            grouped_indices[allele].add(idx)

        grouped_indices = [tuple(v) for v in grouped_indices.values()]
        
        if len(grouped_indices) > 1:
            for a1, a2 in combinations(grouped_indices, 2):
                pair_count = get_pair_count(a1, a2)
                explained_reads[(a1, a2)] = pair_count/gene_count
        else:
            a1, a2 = sorted(list(grouped_indices)[0])[:2]
            pair_count = get_pair_count(a1, a2)             
            explained_reads[((a1,),(a2,))] = pair_count/gene_count
            
            
        # Print information
        log.info('\n[genotype] Pairs by % explained reads:')
        log.info('\t\t{: <28}    {: >7}\t'.format('allele pair', 'explained'))
        for (a1,a2), count in sorted(explained_reads.items(), 
                                 key=lambda x: x[1], 
                                 reverse = True):
            alleles = ', '.join([process_allele(allele_idx[a1[0]][0], 3),
                             process_allele(allele_idx[a2[0]][0], 3)])
            log.info('\t\t{: <28}    {: >9.2f}%\t'
                     .format(alleles, count*100))
            
            
        max_count = max(explained_reads.values())
        top_by_reads = {pair:count for pair,count in explained_reads.items() 
                       if count == max_count}
    
        # If more than one pair has the same number of explained reads
        # use allele frequency priors to break the tie
        pair_prior = dict()
        
        if len(top_by_reads) > 1 and population:
            for a1,a2 in top_by_reads.keys():
                allele1 = process_allele(allele_idx[a1[0]][0],2)
                allele2 = process_allele(allele_idx[a2[0]][0],2)
                if allele1 not in prior or allele2 not in prior:
                    continue

                pair_prior[(a1,a2)] =   prior[allele1][population] \
                                      * prior[allele2][population]
            
            
        if pair_prior:
            max_prior = max(pair_prior.values())
            pair_prior = {pair:prior for pair,prior in pair_prior.items() 
                          if prior >= max_prior}
                
            a1, a2 = sorted(pair_prior.keys(), key = lambda x: (x[0], x[1]))[0]
            
        else:
            a1, a2 = sorted(top_by_reads.items(), 
                              key = lambda x: x[1], 
                              reverse = True)[0][0]

        pair_count = get_pair_count(a1, a2)
        a1_count, a2_count = get_nonshared_count(a1, a2)
        
        a1 = process_allele(allele_idx[sorted(a1)[0]][0], 3)
        a2 = process_allele(allele_idx[sorted(a2)[0]][0], 3)
        # Zygosity check based on nonshared counts
        log.info(f'\n[genotype] Checking zygosity')
        if a1_count == a2_count == 0:
            log.info('[genotype] Unable to distinguish ' +
                     'between minor and major alleles')
            genotype = [a1, a2]
        elif a1_count == 0:
            log.info('[genotype] Likely heterozygous: minor allele has no '+
                     'nonshared reads')
            genotype = [a2, a2]
        elif a2_count == 0:
            log.info('[genotype] Likely heterozygous: minor allele has no '+
                     'nonshared reads')
            genotype = [a1, a1]
        elif min(a1_count/a2_count, a2_count/a1_count) < zygosity_threshold:
            log.info(f'[genotype] Likely homozygous: minor/major '+
                      'nonshared count {:.2f}'
                      .format(min(a1_count/a2_count, a2_count/a1_count)))
            if a1_count > a2_count:
                genotype = [a1, a1]
            else:
                genotype = [a2, a2]
        else:
            log.info(f'[genotype] Likely heterozygous: minor/major '+
                      'nonshared count {:.2f}'
                      .format(min(a1_count/a2_count, a2_count/a1_count)))
            genotype = [a1,a2]
        
        
        
    else:
        a1, alleles, _ = em_results[0]
        pair_count = get_count(a1)
        a1_count = pair_count
        a2_count = None
        genotype = [process_allele(alleles[0], 3),process_allele(alleles[0], 3)]
        
    return genotype, pair_count

def genotype_gene(gene, gene_count, eqs, lengths, allele_idx, population, 
                  prior, tolerance, max_iterations, drop_iterations, 
                  drop_threshold, zygosity_threshold):
    '''Calls transcript quantification and genotype prediction.'''
    
    if gene not in {'A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1'}:
        population = None

    allele_idx = json.loads(allele_idx)
    lengths = json.loads(lengths)
    lengths = dict([a, int(x)] for a, x in lengths.items())

    em_results = expectation_maximization(eqs, 
                                          lengths, 
                                          allele_idx, 
                                          population, 
                                          prior, 
                                          tolerance, 
                                          max_iterations, 
                                          drop_iterations, 
                                          drop_threshold)


    em_results = [[idx, allele_idx[idx], a] for idx, a in em_results.items()]
    
    log.info('\n[genotype] Top alleles by abundance:')
    log.info('\t\t{: <20}    {: >9}'.format('allele', 'abundance'))
    
    for _, alleles, abundance in sorted(em_results, 
                                          key=lambda x: x[2], 
                                          reverse = True):
                                     
        log.info('\t\t{: <20}    {: >8.2f}%'
                 .format(process_allele(alleles[0], 3), abundance*100))
    
    genotype, pair_count = predict_genotype(eqs, 
                                           allele_idx, 
                                           allele_eq, 
                                           em_results,
                                           gene_count,
                                           population, 
                                           prior, 
                                           zygosity_threshold)

    
    
    log.info('\n[genotype] Most likely genotype explaining {:.0f} reads:'
            .format(pair_count))
    
    for allele in genotype:
        log.info(f'\t\t{allele}')
    
    return em_results, genotype

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

def arg_check_files(parser, arg):
    accepted_formats = ('alignment.p','fq.gz','fastq.gz','fq','fastq')
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error('The file %s does not exist.' %file)
        elif not file.lower().endswith(accepted_formats):
            parser.error('The format of %s is invalid.' %file)
        return arg
        
def arg_check_genes(parser, arg):
    if arg.lower() == 'all':
        return sorted(genes)
    input_genes = {gene.upper() for gene in arg.split(',')} & genes
    if not input_genes:
        parser.error('The gene list %s is invalid.' %arg)
    return sorted(input_genes)

def arg_check_population(parser, arg):
    if arg.lower() == 'none':
        return None
    if arg not in populations:
        parser.error('The population %s is invalid.' %arg)
    return arg
    
def arg_check_tolerance(parser, arg):
    try:
        value = float(arg)
        if value > 1 or value < 0:
            parser.error('The tolerence must be between 0 and 1.')
        return value
    except:
        parser.error('The tolerence must be a floating point number.')
    
def arg_check_iterations(parser, arg):
    try:
        value = int(arg)
        if value < 0:
            parser.error('The number of iterations must be positive.')
        return value
    except:
        parser.error('The number of iterations must be an integer.')
    
def arg_check_threshold(parser, arg):
    try:
        value = float(arg)
        if value > 1 or value < 0:
            parser.error('The threshold must be between 0 and 1.')
        return value
    except:
        parser.error('The threshold is invalid.')

if __name__ == '__main__':
    
    with open(parameters_json, 'r') as file:
        genes, populations, _ = json.load(file)
        genes = set(genes)
        populations = set(populations)
    
    parser = argparse.ArgumentParser(prog='arcasHLA genotype',
                                 usage='%(prog)s [options] FASTQs ' + 
                                        'or alignment.p file',
                                 add_help=False,
                                 formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        help='list of fastq files (e.g. ' +
                              'sample.extracted.fq.gz) or alignment file ' + 
                              '(sample.alignment.p)', 
                        nargs='*',
                        type=lambda x: arg_check_files(parser, x))
                        
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
                        
    parser.add_argument('--log', 
                        type=str,
                        help='log file for run summary\n'+
                             'default: sample.genotype.log\n\n',
                        default=None, 
                        metavar='')
    
    parser.add_argument('-g',
                        '--genes',
                        help='comma separated list of HLA genes\n'+
                             'default: all\n' + '\n'.join(wrap('options: ' +
                             ', '.join(sorted(genes)), 60)) +'\n\n',
                        default='all', 
                        metavar='',
                        type=lambda x: arg_check_genes(parser, x))
    
    parser.add_argument('-p',
                        '--population', 
                        help= 'sample population\ndefault: prior\n' + 
                              '\n'.join(wrap('options: ' + 
                              ', '.join(sorted(populations)), 60)) +'\n\n',
                        default='prior', 
                        metavar='',
                        type=lambda x: arg_check_population(parser, x))
    '''
    parser.add_argument('-d',
                        '--database', 
                        type=str,
                        choices = databases,
                        help='frequency database\n  default: gold_smoothed10\n'+
                             '  options: ' + ', '.join(databases) + '\n\n',
                        default='freq', 
                        metavar='')
    '''
    
    parser.add_argument('--tolerance', 
                        type=lambda x: arg_check_tolerance(parser, x),
                        help='convergence tolerance\n  default: 10e-7\n\n',
                        default=10e-7, 
                        metavar='')
    
    parser.add_argument('--max_iterations', 
                        type=lambda x: arg_check_iterations(parser, x),
                        help='maximum # of iterations\n  default: 1000\n\n',
                        default=1000, 
                        metavar='')
    
    parser.add_argument('--drop_iterations', 
                        type=lambda x: arg_check_iterations(parser, x),
                        help='EM iteration to start dropping low-support ' + 
                             'alleles\n  default: 20\n  recommended paired:' + 
                             '20\n  recommended single: 4\n\n',
                        default=None, 
                        metavar='')
    
    parser.add_argument('--drop_threshold', 
                        type=lambda x: arg_check_threshold(parser, x),
                        help='proportion of max abundance allele needs to not '+
                             'be dropped\n  default: 0.1\n\n',
                        default=0.1, 
                        metavar='')
    
    parser.add_argument('--zygosity_threshold', 
                        type=lambda x: arg_check_threshold(parser, x),
                        help='proportion of major allele abundance needed to ' +
                             'be considered heterozygous\n  default: 0.1\n\n',
                        default=0.15, 
                        metavar='')
    
    parser.add_argument('--min_count', 
                        type=int,
                        help='minimum gene read count required for genotyping ' +
                             '\n  default: 75\n\n',
                        default=75, 
                        metavar='')
    
    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='out directory\n\n',
                        default='./', 
                        metavar='')
    
    parser.add_argument('--temp', 
                        type=str,
                        help='temp directory\n\n',
                        default='/tmp/', 
                        metavar='')
    
    parser.add_argument('--keep_files',
                        action = 'count',
                        help='keep intermediate files\n\n',
                        default=False)
                        
    parser.add_argument('-t',
                        '--threads', 
                        type = str,
                        default='1',
                        metavar='')

    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)

    parser.add_argument('-l',
                        '--avg', 
                        type=int,
                        help='Estimated average fragment length ' +
                                'for single-end reads\n  default: 200\n\n',
                        default=200)

    parser.add_argument('-s',
                        '--std', 
                        type=int,
                        help='Estimated standard deviation of fragment length ' +
                             'for single-end reads\n  default: 20\n\n',
                        default=20)

    parser.add_argument('--single',
                        action='store_true',
                        help='Include flag if single-end reads. Default is paired-end.\n\n',
                        default=False)


    args = parser.parse_args()
    
    if len(args.file) == 0:
        sys.exit('[genotype] Error: FASTQ or alignment.p file required.')
    
    # Set up temporary and output folders, log file
    sample = os.path.basename(args.file[0]).split('.')[0]
    outdir = check_path(args.outdir)
    temp = create_temp(args.temp)
    if args.log:
        log_file = args.log
    else:
        log_file = ''.join([outdir,sample,'.genotype.log'])  
    with open(log_file, 'w'):
        pass
    if args.verbose:
        handlers = [log.FileHandler(log_file), log.StreamHandler()]
        
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
    else:
        handlers = [log.FileHandler(log_file)]
            
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
   
    
    log.info('')
    hline()
    log.info(f'[log] Date: %s', str(date.today()))
    log.info(f'[log] Sample: %s', sample)
    log.info(f'[log] Input file(s): %s', f'\n\t\t     '.join(args.file))
        
    # Load HLA frequencies
    prior = pd.read_csv(hla_freq, delimiter='\t')
    prior = prior.set_index('allele').to_dict('index')
       
    # Checks if HLA reference exists
    check_path(rootDir + 'dat/ref')
    check_ref()
    
    # Loads reference information
    #with open(hla_p, 'rb') as file:
    #    reference_info = pickle.load(file)
    #    (commithash,(gene_set, allele_idx, 
    #     lengths, gene_length)) = reference_info
    with open(hla_json, 'r') as file:
        reference_info = json.load(file)
        (commithash,(gene_set, allele_idx, 
         lengths, gene_length)) = reference_info
        
    log.info(f'[log] Reference: %s', commithash)
    hline()
        
    if args.file[0].endswith('.alignment.p'):
        alignment_info = load_alignment(args.file[0], commithash)
    else:
        alignment_info = get_alignment(args.file, sample, hla_idx,
                                      reference_info, outdir, temp,
                                      args.threads, args.single,
                                      avg=args.avg, std=args.std)
        
    (commithash, eq_idx, allele_eq,
            paired, align_stats, gene_stats) = alignment_info

    # Set up EM parameters
    if not args.drop_iterations:
        if paired: args.drop_iterations = 20
        else: args.drop_iterations = 4
        
    em_results = dict()
    genotypes = dict()
    
    hline()
    log.info('[genotype] Genotyping parameters:')
    log.info(f'\t\tpopulation: %s', args.population)
    log.info(f'\t\tminimum count: %s', args.min_count)
    log.info(f'\t\tmax iterations: %s', args.max_iterations)
    log.info(f'\t\ttolerance: %s', args.tolerance)
    log.info(f'\t\tdrop iterations: %s', args.drop_iterations)
    log.info(f'\t\tdrop threshold: %s', args.drop_threshold)
    log.info(f'\t\tzygosity threshold: %s', args.zygosity_threshold)
    
    # For each HLA locus, perform EM then scoring
    for gene in args.genes:
        hline()
        log.info(f'[genotype] Genotyping HLA-{gene}')
        
        # Skips loci with not enough reads to genotype
        if gene not in gene_stats or gene_stats[gene][0] < args.min_count:
            log.info(f'[genotype] Not enough reads aligned to HLA-{gene} to genotype.') 
            continue
        gene_count, eq_count, abundance = gene_stats[gene]
        log.info(f'[genotype] {gene_count:.0f} reads aligned to HLA-{gene} '+
                 f'in {eq_count} classes')
            
            
        em, genotype = genotype_gene(gene,
                                     gene_count,
                                     eq_idx[gene], 
                                     lengths, 
                                     allele_idx,
                                     args.population, prior, 
                                     args.tolerance, args.max_iterations, 
                                     args.drop_iterations, args.drop_threshold,
                                     args.zygosity_threshold)
                                            
        em_results[gene] = em
        genotypes[gene] = genotype
            
    with open(''.join([outdir, sample, '.genotype.json']), 'w') as file:
            json.dump(genotypes, file)
            
    remove_files(temp, args.keep_files)
    
    hline()
    log.info('')

#-----------------------------------------------------------------------------
