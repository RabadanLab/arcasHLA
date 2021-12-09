#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   partial.py: genotypes partial alleles from extracted reads.
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
import pandas as pd
import logging as log


from argparse import RawTextHelpFormatter
from itertools import combinations
from collections import Counter, defaultdict
from textwrap import wrap
from datetime import date

from reference import check_ref
from arcas_utilities import *
from align import *
from genotype import expectation_maximization

__version__     = '0.3.0'
__date__        = '2021-12-09'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir = os.path.dirname(os.path.realpath(__file__)) + '/../'

partial_json       = rootDir + 'dat/ref/hla_partial.p.json'
partial_idx     = rootDir + 'dat/ref/hla_partial.idx'
hla_freq        = rootDir + 'dat/info/hla_freq.tsv'

#-------------------------------------------------------------------------------
# Process transcript assembly output
#-------------------------------------------------------------------------------

def filter_eqs(complete_genotypes, allele_idx, eq_idx, partial_alleles):
    '''Filters compatibility classes if they contain partial alleles or
       at least one predicted complete allele.
    '''
    
    # All alleles returned from arcasHLA genotype
    all_predicted = {allele for alleles in complete_genotypes.values() 
                     for allele in alleles}
    
    # Find the indices of sequences derived from these alleles
    wanted_indices = {index for index, alleles in allele_idx.items() 
                      if alleles and 
                      (set(alleles) & (partial_alleles | set(all_predicted)))}
    
    # For each eq group, select only those that include the previous indices
    filtered_eqs = dict()
    for group, eq_list in eq_idx.items():
        filtered_eqs[group] = dict()
        for gene in complete_genotypes:
            if gene not in eq_list:
                continue

            filtered = []
            for indices, count in eq_list[gene]:
                indices = set(indices) & wanted_indices
                
                if not indices:
                    continue

                filtered.append([indices,count])
            filtered_eqs[group][gene] = filtered
            
    # Find indices given an allele
    allele_eq = {group:defaultdict(set) for group in filtered_eqs.keys()}
    for group, eq_list in filtered_eqs.items():
        for gene in eq_list:
            for i,(indices, count) in enumerate(eq_list[gene]):
                for idx in indices:
                    for allele in allele_idx[idx]:
                        allele = process_allele(allele,3)
                        allele_eq[group][allele].add(i)

    return filtered_eqs, allele_eq
    
#-------------------------------------------------------------------------------
# Partial genotyping
#-------------------------------------------------------------------------------

def type_partial(eqs, gene, partial_exons, complete_genotype, partial_alleles, 
                 population, prior, tolerance, max_iterations, 
                 drop_iterations, drop_threshold, zygosity_threshold):
    '''Types partial alleles.'''
    
    # Return count of a single allele
    def get_single_count(a):
        return sum([eqs[group][gene][idx][1] for idx in allele_eq[group][a]])

    # Return count of a pair of alleles
    def get_pair_count(a1, a2):
        indices = allele_eq[group][a1] | allele_eq[group][a2]
        return sum([eqs[group][gene][idx][1] for idx in indices])

    # Return nonshared count of a pair of alleles
    def get_nonshared_count(a1, a2):
        a1_indices = allele_eq[group][a1] - allele_eq[group][a2]
        a2_indices = allele_eq[group][a2] - allele_eq[group][a1]
        a1_count = sum([eqs[group][gene][idx][1] for idx in a1_indices])
        a2_count = sum([eqs[group][gene][idx][1] for idx in a2_indices])
        return a1_count, a2_count
    
    if gene not in {'A', 'B', 'C', 'DRB1', 'DQB1', 'DQA1'}:
        population = None

    # Set binding region by class
    if gene.startswith('D'): 
        binding_region = "['2']"
    else:
        binding_region = "['2', '3']"
        
    if gene not in eqs[binding_region]:
        log.info(f'[genotype] No reads aligned to HLA-{gene} binding region') 
        return complete_genotype
        
    # Get group of possible partial alleles by performing transcript 
    # quantification on the binding region exons
    results  = expectation_maximization(eqs[binding_region][gene], 
                                        lengths, 
                                        allele_idx, 
                                        population, 
                                        prior, 
                                        tolerance, 
                                        max_iterations, 
                                        drop_iterations, 
                                        drop_threshold)
    exon_groups = defaultdict(set)
    
    
    # Map partial alleles to their possible exon combinations
    for idx in results:
        alleles = {process_allele(allele,3) for allele in allele_idx[idx]}
        for allele in (alleles - set(complete_genotype)) & partial_alleles:
            for group in eqs.keys():
                if group[1:-1] in str(sorted(partial_exons[allele].keys())):
                    exon_groups[group].add(allele)
                    
    # Compare pairs of partial alleles and predicted alleles
    overall = []
    for group in sorted(exon_groups.keys(), 
                        key = lambda x: len(x),
                        reverse=False):
    
        # Skip just exon 2 for class I
        if not gene.startswith('D') and group == "['2']":
            continue
            
        # Only look at partial alleles that have a different sequence for 
        # this combination of exons than the complete alleles
        a1, a2 = complete_genotype
        possible_alleles = {allele for allele in exon_groups[group] 
                            if allele_eq[group][allele] != allele_eq[group][a1] 
                            and allele_eq[group][allele] != allele_eq[group][a2]}
        
        if not possible_alleles:
            continue

        explained_reads = dict()
        
        # Filter partial alleles that have only a few more reads 
        # than the complete alleles
        min_count = min(get_single_count(a1),get_single_count(a2))
        possible_alleles &= {allele for allele in exon_groups[group] 
                             if get_single_count(allele) > 10 + min_count}

        total_count = sum([count for _,count in eqs[group][gene]])

        # Get percent explained reads by complete genotype
        pair_count = get_pair_count(a1, a2)
        explained_percent = round(pair_count/total_count,8)
        explained_reads[(a1, a2)] = explained_percent

        # Only consider pairs with partial alleles if they explain 
        # a greater percentage of reads
        for a1, a2 in combinations(set(complete_genotype) | possible_alleles, 2):
            pair_count = get_pair_count(a1, a2)
            if pair_count/total_count > explained_percent:
                explained_reads[(a1,a2)] = round(pair_count/total_count, 8)
       
        if not explained_reads:
            continue 
        top_perc = max(explained_reads.values())
        explained_reads = {key:value for key,value in explained_reads.items() 
                           if value == top_perc}
        
        # If the top percentage of explained reads is shared by more than 
        # one pair, use priors to break the ties
        if population and len(explained_reads) > 1:
            pair_prior = dict()
            for a1,a2 in explained_reads.keys():
                if (process_allele(a1,2) not in prior 
                        or process_allele(a2,2) not in prior):
                    continue
                    
                pair_prior[(a1,a2)] =  prior[process_allele(a1,2)][population] \
                                     * prior[process_allele(a2,2)][population]
                                     
            if pair_prior:    
                a1,a2 = sorted(pair_prior.items(), 
                               key=lambda x:x[1], 
                               reverse=True)[0][0]
            else:
                a1,a2 = sorted(explained_reads.items(),
                                              key= lambda x:x[1],
                                              reverse=True)[0][0]
        else:
            a1,a2 = sorted(explained_reads.items(),
                                              key= lambda x:x[1],
                                              reverse=True)[0][0]
        group = re.sub('[\'\[\]]','',group)
        log.info('\t\texons {: <22}\t{: <28}\t{:.2f}%'
                 .format(group, ', '.join([a1,a2]), top_perc*100))
            
        overall.append([(a1,a2),top_perc])

    if overall:
        return sorted(overall,
                      key=lambda x: (x[1],x[0][0],x[0][1]), 
                      reverse = True)[0][0]
    
    return complete_genotype

#-------------------------------------------------------------------------------
# Runs partial genotyping
#-------------------------------------------------------------------------------

def arg_check_files(parser, arg):
    accepted_formats = ('alignment.p','fq.gz','fastq.gz','fq','fastq')
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error('The file %s does not exist.' %file)
        elif not file.lower().endswith(accepted_formats):
            parser.error('The format of %s is invalid.' %file)
        return arg
    
def arg_check_genotype(parser, arg):
    if not os.path.isfile(arg):
        parser.error('The genotype file %s does not exist.' %arg)
    elif not arg.endswith('.genotype.json'):
        parser.error('The genotype file %s is invalid' %arg)
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
    
    with open(rootDir + 'dat/info/parameters.json', 'r') as file:
        genes, populations, databases = json.load(file)
        genes = set(genes)
        populations = set(populations)
    
    parser = argparse.ArgumentParser(prog='arcasHLA partial',
                             usage='%(prog)s [options] -G genotype.json FASTQ',
                             add_help=False,
                             formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        type=lambda x: arg_check_files(parser, x), 
                        help='list of fastq files or ".partial.json" file', 
                        nargs='*')
                        
                        
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('-G',
                    '--genotype', 
                    type=lambda x: arg_check_genotype(parser, x),
                    help='"genotype.json" file from arcasHLA genotype'
                    , metavar='')
   
                        
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
                        help= 'sample population\n  default: prior\n' + 
                              '\n  '.join(wrap('  options: ' + 
                              ', '.join(sorted(populations)), 60)) +'\n\n',
                        default='prior', 
                        metavar='',
                        type=lambda x: arg_check_population(parser, x))
    
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
                             'alleles\n  default: 4\n\n',
                        default=4, 
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
    
    
    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='out directory\n\n',
                        default='./', metavar='')
    
    parser.add_argument('--temp', 
                        type=str,
                        help='temp directory\n\n',
                        default='/tmp/', metavar='')
    
    parser.add_argument('--keep_files',
                        action = 'count',
                        help='keep intermediate files\n\n',
                        default=False)
                        
    parser.add_argument('-t',
                        '--threads', 
                        type = str,
                        default='1')

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
        sys.exit('[genotype] Error: FASTQ or partial_alignment.p file required')
    
    # Set up directories and log file
    sample = os.path.basename(args.file[0]).split('.')[0]
    outdir = check_path(args.outdir)
    temp = create_temp(args.temp)
    if args.log:
        log_file = args.log
    else:
        log_file = ''.join([outdir,sample,'.partial_genotype.log'])    
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
    log.info(f'[log] Input file(s): %s', ', '.join(args.file))
        
    
    prior = pd.read_csv(hla_freq, delimiter='\t')
    prior = prior.set_index('allele').to_dict('index')
       
    # Checks if HLA reference exists
    check_ref()
    
    # Loads reference information
    #with open(partial_p, 'rb') as file:
    #    reference_info = pickle.load(file)
    #    (commithash, (gene_set, allele_idx, exon_idx, 
    #        lengths, partial_exons, partial_alleles)) = reference_info
    with open(partial_json, 'r') as file:
        reference_info = json.load(file)
        (commithash, (gene_set, allele_idx, exon_idx, 
            lengths, partial_exons, partial_alleles)) = reference_info
        gene_set = set(gene_set)
        allele_idx = json.loads(allele_idx)
        exon_idx = json.loads(exon_idx)
        lengths = json.loads(lengths)
        lengths = dict([a, int(x)] for a, x in lengths.items())
        partial_exons = json.loads(partial_exons)
        partial_alleles = set(partial_alleles)
        
    log.info(f'[log] Reference: %s', commithash)
    hline()  
          
    # Runs transcript assembly if intermediate json not provided
    if args.file[0].endswith('.partial_alignment.p'):
        alignment_info = load_alignment(args.file[0], commithash, True)
    else:
        alignment_info = get_alignment(args.file, sample, partial_idx,
                                       reference_info, outdir, temp, 
                                       args.threads, args.single, True, args.avg, args.std)
    commithash, eq_idx, _, paired, align_stats, _ = alignment_info
     
    # Load alleles from arcasHLA genotype
    with open(args.genotype,'r') as file:
        complete_genotypes = json.load(file)
    
    genes = set(args.genes) & set(complete_genotypes.keys())
    
    # Filter out alleles not returned by arcasHLA genotype
    eq_idx, allele_eq = filter_eqs(complete_genotypes, allele_idx, 
                                   eq_idx, partial_alleles)
    
    partial_results = dict()

    # For each HLA locus, check for possible partial alleles
    for gene in sorted(genes):
        hline()
        log.info(f'[genotype] Partial genotyping HLA-{gene}')
        complete_genotype = complete_genotypes[gene]
        if len(complete_genotype) == 1:
            complete_genotype.append(complete_genotype[0])
                                     
                                        
        genotype = type_partial(eq_idx, 
                                gene,
                                partial_exons,
                                complete_genotype, 
                                partial_alleles, 
                                args.population, 
                                prior, 
                                args.tolerance, 
                                args.max_iterations, 
                                args.drop_iterations, 
                                args.drop_threshold,
                                args.zygosity_threshold)
                                
        partial_results[gene] = genotype
        log.info('\n[genotype] Most likely genotype:')
    
        for allele in genotype:
            log.info(f'\t\t{allele}')
    
            
    with open(''.join([outdir, sample,'.partial_genotype.json']),'w') as file:
            json.dump(partial_results, file)
            
    remove_files(temp, args.keep_files)
            
    hline()
    log.info('')

#-----------------------------------------------------------------------------
