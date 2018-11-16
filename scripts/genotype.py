#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   genotype.py: genotypes from extracted chromosome 6 reads.
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
from arcas_utilities import (process_allele, check_path, remove_files, 
                            run_command, get_gene)

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

hla_p    = 'database/hla.p'
hla_idx  = 'database/hla.idx'
hla_freq = 'database/hla_freq.tsv'

#-----------------------------------------------------------------------------
# Process and align FASTQ input
#-----------------------------------------------------------------------------

def analyze_reads(fqs, temp, reads_file, keep_files):
    '''Analyzes read length for single-end sampled, required by Kallisto.'''
    awk = "| awk '{if(NR%4==2) print length($1)}'"
    
    log.info('[align] analyzing read length')
    fq = fqs[0]

    command = ['zcat <', fq, awk, '>', reads_file]
    run_command(command)
        
    read_lengths = np.genfromtxt(reads_file)
    
    avg = round(np.mean(read_lengths), 4)
    std = round(np.std(read_lengths), 4)
    
    remove_files([reads_file], keep_files)
    
    return str(avg), str(std)

def pseudoalign(fqs, sample, paired, reference, outdir, temp, 
                         threads, keep_files, partial = False):
    '''Calls Kallisto to pseudoalign reads.'''
    file_list = []
    
    # Get read length stats for single-end sample
    if not paired:
        reads_file = ''.join([temp, sample, '.reads.txt'])
        avg, std = analyze_reads(fqs, temp, reads_file, keep_files)
        if std == 0: std = .00001

    temp2 = check_path(''.join([temp, sample]))
    command = ['kallisto pseudo -i', reference, '-t', threads, '-o', temp2]
        
    if paired:
        command.extend([fqs[0], fqs[1]])
    else:
        fq = fqs[0]
        command.extend(['--single -l', avg, '-s', std, fq])
        
    run_command(command, '[align] pseudoaligning with Kallisto: ')

    # Move and rename Kallisto output
    file_in = ''.join([temp2, 'pseudoalignments.tsv'])
    count_file = ''.join([temp, sample, '.counts.tsv'])
    file_list.append(file_in)
    run_command(['mv', file_in, count_file])
    
    file_in = ''.join([temp2, 'pseudoalignments.ec'])
    eq_file = ''.join([temp, sample, '.eq.tsv'])
    file_list.append(file_in)
    run_command(['mv', file_in, eq_file])

    run_command(['rm -rf', temp2])
    remove_files(file_list, keep_files)
           
    return count_file, eq_file

#-----------------------------------------------------------------------------
# Process transcript assembly output
#-----------------------------------------------------------------------------

def process_counts(count_file, eq_file, gene_list, allele_idx, 
                   allele_lengths, keep_files): 
    '''Processes pseudoalignment output, returning compatibility classes.'''
    log.info('[align] processing pseudalignment')
    # Process count information
    counts = dict()
    with open(count_file, 'r') as file:
        for line in file.read().splitlines():
            eq, count = line.split('\t')
            counts[eq] = float(count)

    
    # Process compatibility classes
    eqs = dict()
    with open(eq_file, 'r') as file:
        for line in file.read().splitlines():
            eq, indices = line.split('\t')
            eqs[eq] = indices.split(',')

    # Set up compatibility class index
    eq_idx = defaultdict(list)
    
    count_unique = 0
    count_multi = 0
    class_unique = 0
    class_multi = 0
    
    for eq, indices in eqs.items():
        if [idx for idx in indices if not allele_idx[idx]]:
            continue

        genes = list({get_gene(allele) for idx in indices 
                        for allele in allele_idx[idx]})
        count = counts[eq]
        
        if len(genes) == 1 and counts[eq] > 0:
            gene = genes[0]
            eq_idx[gene].append((indices, count))
            
            count_unique += count
            class_unique += 1
        else:
            count_multi += count
            class_multi += 1

    # Alleles mapping to their respective compatibility classes
    allele_eq = defaultdict(set)
    for eqs in eq_idx.values():
        for eq,(indices,_) in enumerate(eqs):
            for idx in indices:
                allele_eq[idx].add(eq)
            
    remove_files([count_file, eq_file], keep_files)
    
    align_stats = [count_unique, count_multi, class_unique, class_multi]
    
    return eq_idx, allele_eq, align_stats
           

def get_count_stats(eq_idx, gene_length):
    '''Returns counts and relative abundance of genes.'''
    stats = {gene:[0,0,0.] for gene in eq_idx}
    
    abundances = defaultdict(float)
    for gene, eqs in eq_idx.items():
        count = sum([count for eq,count in eqs])
        abundances[gene] = count / gene_length[gene]
        stats[gene][0] = count
        stats[gene][1] = len(eqs)
        
    total_abundance = sum(abundances.values())

    for gene, abundance in abundances.items():
        stats[gene][2] = abundance / total_abundance

    return stats
    
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
        for alleles, count in eqs:
            if population:
                allele_prior = defaultdict(float)
                for idx in alleles:
                    allele = process_allele(allele_idx[idx][0], 2)
                    if allele in prior:
                        allele_prior[idx] = prior[allele][population]

                if allele_prior:
                    divide_prior(alleles, count, allele_prior)
                    continue
                    
            divide_equally(alleles, count)
                
        return counts_to_abundances(counts)
    
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

    theta0 = initial_abundances(eqs, lengths, population)

    log.info(f'[genotype] quantifying allele transcript abundance')
    
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

def predict_genotype(eqs, allele_idx, allele_eq, em_results, population, 
                     prior, zygosity_threshold):
                     
    '''Predicts most likely genotype using scoring based on proportion of 
       explained reads, tie-breaking with allele priors.
    '''
    # Returns number of reads explained by an allele
    def get_count(a):
        observed_eqs = set.union(*[allele_eq[idx] for idx in a])
        return sum([eqs[idx][1] for idx in observed_eqs])
    
    # Returns number of reads explained by a pair of alleles
    def get_pair_count(a1, a2):
        
        a1_eqs = set.union(*[allele_eq[idx] for idx in a1])
        a2_eqs =set.union(*[allele_eq[idx] for idx in a2])
        
        observed_eqs = a1_eqs | a2_eqs
        
        return sum([eqs[idx][1] for idx in observed_eqs])
        
    # Returns non-shared counts for a pair of alleles
    def get_nonshared_count(a1, a2):
        a1_eqs = set.union(*[allele_eq[idx] for idx in a1])
        a2_eqs =set.union(*[allele_eq[idx] for idx in a2])
        
        a1_eqs = a1_eqs - a2_eqs
        a2_eqs = a2_eqs - a1_eqs
        
        a1_count = sum([eqs[idx][1] for idx in a1_eqs])
        a2_count = sum([eqs[idx][1] for idx in a2_eqs])
        
        return a1_count, a2_count
        
    if len(em_results) > 1:
        grouped_indices = defaultdict(set)
        for idx, alleles, abundances in em_results:
            allele = process_allele(alleles[0], 2)
            grouped_indices[allele].add(idx)

        grouped_indices = [tuple(v) for v in grouped_indices.values()]
        
        explained_reads = dict()
        if len(grouped_indices) > 1:
            for a1, a2 in combinations(grouped_indices, 2):
                pair_count = get_pair_count(a1, a2)
                explained_reads[(a1, a2)] = pair_count
        else:
            a1, a2 = sorted(list(grouped_indices)[0])[:2]
            pair_count = get_pair_count(a1, a2)             
            explained_reads[((a1,),(a2,))] = pair_count
            
        max_count = max(explained_reads.values())
        explained_reads = {pair:count for pair,count in explained_reads.items() 
                           if count == max_count}
    
        # If more than one pair has the same number of explained reads
        # use allele frequency priors to break the tie
        if len(explained_reads) > 1 and population:
            pair_prior = dict()
            for a1,a2 in explained_reads.keys():
                allele1 = process_allele(allele_idx[a1[0]][0],2)
                allele2 = process_allele(allele_idx[a2[0]][0],2)
                if allele1 not in prior or allele2 not in prior:
                    continue

                pair_prior[(a1,a2)] =   prior[allele1][population] \
                                      * prior[allele2][population]
                                          
            max_prior = max(pair_prior.values())
            pair_prior = {pair:prior for pair,prior in pair_prior.items() 
                          if prior >= max_prior}
                
            a1, a2 = sorted(pair_prior.keys(), key = lambda x: (x[0], x[1]))[0]
            pair_count = get_pair_count(a1, a2)
            
        else:
            (a1, a2), pair_count = sorted(explained_reads.items(), 
                                           key = lambda x: x[1], 
                                           reverse = True)[0]

        a1_count, a2_count = get_nonshared_count(a1, a2)
        
        a1 = [sorted(a1)[0]]
        a2 = [sorted(a2)[0]]
        
        # Zygosity check based on nonshared counts
        if a1_count == a2_count == 0:
            genotype = a1 + a2
        elif a1_count == 0:
            genotype = a2
        elif a2_count == 0:
            genotype = a1
        elif min(a1_count/a2_count, a2_count/a1_count) < zygosity_threshold:
            if a1_count > a2_count:
                genotype = a1
            else:
                genotype = a2
        else:
            genotype = a1 + a2
        
        genotype = [process_allele(allele_idx[idx][0], 3) for idx in genotype]
        
        
    else:
        a1, alleles, _ = em_results[0]
        pair_count = get_count(a1)
        genotype = [process_allele(alleles[0], 3),]
        
    return genotype, pair_count

def genotype_gene(gene, eqs, lengths, allele_idx, population, prior, tolerance, 
                  max_iterations, drop_iterations, drop_threshold, 
                  zygosity_threshold):
    '''Calls transcript quantification and genotype prediction.'''
    
    if gene not in {'A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1'}:
        population = None

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
    
    genotype, count = predict_genotype(eqs, 
                                       allele_idx, 
                                       allele_eq, 
                                       em_results, 
                                       population, 
                                       prior, 
                                       zygosity_threshold)

    log.info('[genotype] top alleles by abundance:')
    log.info('           {: <20}\t{: >9}\t'.format('allele', 'abundance'))
    
    for _, alleles, abundance in sorted(em_results, 
                                          key=lambda x: x[2], 
                                          reverse = True):
                                     
        log.info('           {: <20}\t{:>8.2f}%\t'
                 .format(process_allele(alleles[0], 3), abundance*100))
    
    log.info(f'[genotype] most likely genotype explaining {count:.0f} reads:')
    
    for allele in genotype:
        log.info(f'           {allele}')
    
    return em_results, genotype

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

def arg_check_files(parser, arg):
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error('The file %s does not exist.' %file)
        elif not (file.endswith('alignment.p') or file.endswith('.fq.gz')):
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
    
if __name__ == '__main__':
    
    with open('database/parameters.p', 'rb') as file:
        genes, populations, databases = pickle.load(file)
    
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
                        help='log file for run summary\n  '+
                             'default: sample.genotype.log\n\n',
                        default=None, 
                        metavar='')
    
    parser.add_argument('-g',
                        '--genes',
                        help='comma separated list of HLA genes\n  '+
                             'default: all\n' + '\n  '.join(wrap('  options: ' +
                             ', '.join(sorted(genes)), 60)) +'\n\n',
                        default='all', 
                        metavar='',
                        type=lambda x: arg_check_genes(parser, x))
    
    parser.add_argument('-p',
                        '--population', 
                        help= 'sample population\n  default: prior\n' + 
                              '\n  '.join(wrap('  options: ' + 
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
                        type=float,
                        help='convergence tolerance\n  default: 10e-7\n\n',
                        default=10e-7, 
                        metavar='')
    
    parser.add_argument('--max_iterations', 
                        type=int,
                        help='maximum # of iterations\n  default: 1000\n\n',
                        default=1000, 
                        metavar='')
    
    parser.add_argument('--drop_iterations', 
                        type=int,
                        help='EM iteration to start dropping low-support ' + 
                             'alleles\n  default: 20\n  recommended paired:' + 
                             '20\n  recommended single: 4\n\n',
                        default=None, 
                        metavar='')
    
    parser.add_argument('--drop_threshold', 
                        type=float,
                        help='proportion of max abundance allele needs to not '+
                             'be dropped\n  default: 0.1\n\n',
                        default=0.1, 
                        metavar='')
    
    parser.add_argument('--zygosity_threshold', 
                        type=float,
                        help='proportion of major allele abundance needed to ' +
                             'be considered heterozygous\n  default: 0.1\n\n',
                        default=0.15, 
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

    args = parser.parse_args()
    
    sample = os.path.basename(args.file[0]).split('.')[0]
    temp, outdir = [check_path(path) for path in [args.temp, args.outdir]]
    
    if args.log:
        log_file = args.log
    else:
        log_file = ''.join([outdir,sample,'.genotype.log'])
    
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
   
    
    if len(args.file) == 0:
        sys.exit('[genotype] Error: FASTQ or alignment.p file required')
    
    log.info('')
    log.info(f'[log] Date: %s', str(date.today()))
    log.info(f'[log] Sample: %s', sample)
    log.info(f'[log] Input file(s): %s', ', '.join(args.file))
        
    prior = pd.read_csv(hla_freq, delimiter='\t')
    prior = prior.set_index('allele').to_dict('index')
       
    # checks if HLA reference exists
    check_ref()
    
    # loads reference information
    with open(hla_p, 'rb') as file:
        (commithash,
        (gene_set, allele_idx, lengths, gene_length)) = pickle.load(file)
        
    log.info(f'[log] Reference: %s\n', commithash)
        
    # runs transcript assembly if intermediate json not provided
    reference = hla_idx
    if not args.file[0].endswith('.alignment.p'):
        
        paired = True if len(args.file) == 2 else False
        
        count_file, eq_file = pseudoalign(args.file,
                                          sample,
                                          paired,
                                          reference,
                                          outdir, 
                                          temp,
                                          args.threads,
                                          args.keep_files)
        
        eq_idx, allele_eq, align_stats = process_counts(count_file,
                                                         eq_file, 
                                                         gene_set, 
                                                         allele_idx, 
                                                         lengths,
                                                         args.keep_files)
        
        gene_stats = get_count_stats(eq_idx, gene_length)

        with open(''.join([outdir, sample, '.alignment.p']), 'wb') as file:
            output = [commithash, eq_idx, allele_eq, paired, align_stats, gene_stats]
            pickle.dump(output, file)
            
    else:
        log.info(f'\n[align] loading previous alignment %s', args.file[0])
        with open(args.file[0], 'r') as file:
            (commithash_alignment, eq_idx, allele_eq, paired, 
             align_stats, gene_stats) = pickle.load(file)
             
        if commithash != commithash_alignment:
            sys.exit('[genotype] Error: reference used for alignment ' +
                     'different than the one in the database')

    if not args.drop_iterations:
        if paired: args.drop_iterations = 20
        else: args.drop_iterations = 4
        
    em_results = dict()
    genotypes = dict()
    
    count_unique, count_multi, class_unique, class_multi = align_stats
    log.info('[align] {:.0f} reads aligned to HLA reference'
             .format(count_unique + count_multi))
              
    log.info('[align] {:.0f} reads map to a single HLA gene'
             .format(count_unique))
    
    # Print gene abundance information
    log.info('[align] observed HLA genes by overall abundance:')

    log.info('           {: <20}\t{: >9}\t{: >9}'
             .format('gene','abundance','count'))

    for g,(c,_,a) in sorted(gene_stats.items()):
        log.info('           HLA-{: <16}\t{:>8.2f}%\t{:>9.0f}'
                 .format(g, a*100, c))
        

    log.info(f'\n[genotype] population: %s', args.population)
    log.info(f'[genotype] max iterations: %s', args.max_iterations)
    log.info(f'[genotype] tolerance: %s', args.tolerance)
    log.info(f'[genotype] drop iterations: %s', args.drop_iterations)
    log.info(f'[genotype] drop threshold: %s', args.drop_threshold)
    log.info(f'[genotype] zygosity threshold: %s', args.zygosity_threshold)
    
    for gene in args.genes:
        log.info(f'\n[genotype] genotyping HLA-{gene}')
        
        if gene not in gene_stats or gene_stats[gene][0] == 0:
            log.info(f'[genotype] no reads aligned') 
            continue
        total_count, eq_count, abundance = gene_stats[gene]
        log.info(f'[genotype] {total_count:.0f} reads aligned')
            
        em, genotype = genotype_gene(gene, 
                                     eq_idx[gene], 
                                     lengths, 
                                     allele_idx,
                                     args.population, prior, 
                                     args.tolerance, args.max_iterations, 
                                     args.drop_iterations, args.drop_threshold,
                                     args.zygosity_threshold)
                                            
        em_results[gene] = em
        genotypes[gene] = genotype
        
    with open(''.join([outdir, sample, '.em.json']), 'w') as file:
            json.dump(em_results, file)
            
    with open(''.join([outdir, sample, '.genotype.json']), 'w') as file:
            json.dump(genotypes, file)
            
    log.info('')

#-----------------------------------------------------------------------------
