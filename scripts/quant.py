#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   quant.py: genotypes from extracted chromosome 6 reads.
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
import time

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
                            get_gene, hline)

from subprocess import PIPE, run

__version__     = '1.0'
__date__        = 'November 2018'

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

def analyze_reads(fqs, paired, reads_file, keep_files):
    '''Analyzes read length for single-end sampled, required by Kallisto.'''
    awk = "| awk '{if(NR%4==2) print length($1)}'"
    
    #log.info('[alignment] Analyzing read length')
    if paired:
        fq1, fq2 = fqs

        command = ['zcat <', fq1, awk, '>' , reads_file]
        run_command(command)
        
        command = ['zcat <',fq2,awk, '>>', reads_file]
        run_command(command)
        
    else:
        fq = fqs[0]
        command = ['zcat <',fq,awk,'>',reads_file]
        run_command(command)
        
    read_lengths = np.genfromtxt(reads_file)
    
    num = len(read_lengths)
    avg = round(np.mean(read_lengths), 4)
    std = round(np.std(read_lengths), 4)
    
    remove_files([reads_file], keep_files)
    
    return num, avg, std

def run_command(command, message = ''):
    '''Outputs message and command to log, runs command and returns output.'''
    if type(command) == list:
        command = ' '.join(command)

    print(command)
    output = run(command, shell=True, stdout=PIPE, stderr=PIPE)
        
    return output

def arg_check_files(parser, arg):
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error('The file %s does not exist.' %file)
        elif not (file.endswith('alignment.p') or file.endswith('.fq.gz') or file.endswith('.fastq.gz')):
            parser.error('The format of %s is invalid.' %file)
        return arg
    
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
        
        counts = defaultdict(float)
        for alleles, count in eqs:
                    
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
                #print('dropping eq')
                continue

            for allele in alleles:
                counts[allele] += count * (abundances[allele]/total_abundance)
                
        return counts_to_abundances(counts)
    
    def abundance_to_counts(eqs, abundances):
        
        counts = defaultdict(float)
        
        for alleles, count in eqs:
            alleles = [allele for allele in alleles if allele in abundances]
            total_abundance = sum([abundances[allele] for allele in alleles])
            
            if total_abundance == 0:
                continue
            
            for allele in alleles:
                counts[allele] += count * (abundances[allele]/total_abundance)
                
        return counts
    
    # Drop low support alleles after a specified number of iterations if their 
    # abundance is less than a specified proportion of the greatest abundance
    def drop_alleles(eqs, abundances, drop_iterations, drop_threshold, 
                     iterations, converged):
        
        drop_threshold = 1e-6
        eqs_ = []
        #if iterations >= drop_iterations or converged:
        if min(abundances.values()) < drop_threshold:
            #print('dropping alleles')
            #threshold = drop_threshold * max(abundances.values())
            abundances = {allele:abundance 
                          for allele, abundance in abundances.items() 
                          if abundance >= drop_threshold}
            
            kept_indices = set(abundances.keys())
            for indices,count in eqs:
                indices = kept_indices & set(indices)
                if indices:
                    eqs_.append([tuple(indices),count])
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
    
    #log.info('[genotype] Top 10 alleles by undivided read count:')
    #log.info('\t\t{: <20}    {: >10}\t'.format('allele', 'read count'))
    
    #for idx, count in sorted(undivided_counts.items(), 
    #                         key=lambda x: x[1], 
    #                         reverse = True)[:10]:
                                     
        #log.info('\t\t{: <20}    {: >10.0f}\t'
        #         .format(process_allele(allele_idx[idx][0], 3), count))

    #log.info(f'\n[genotype] Quantifying allele transcript abundance')
    
    # SQUAREM - accelerated EM
    # R. Varadhan & C. Roland (doi: 10.1 1 1 1/j. 1467-9469.2007.00585.X)
    # Used by HISAT-genotype, originaly used by Sailfish
    while iterations < max_iterations and not converged:
        #print(iterations)
        print(iterations, end='\r')
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
                        
                theta_prime[allele] = max(value, 0.0)
                
            #step_min = min(theta_prime.values())
            #step_max = max(theta_prime.values())

            # Adjust step rather than kicking out alleles with a negative result
            #if step_min < 0:
            #    theta_prime = {allele:(value-step_min)/(step_max-step_min) 
            #                   for allele, value in theta_prime.items()}
            #                   
            #    total = sum(theta_prime.values())
            #    
            #    theta_prime = {allele:value/total 
            #                   for allele, value in theta_prime.items()}
     
            # Update abundances with given the new proportions
            theta_prime = update_abundances(eqs, theta_prime)
            
        else:
            theta_prime = theta1
            
        converged = check_convergence(theta0, theta_prime)

        theta0, eqs = drop_alleles(eqs, theta_prime, drop_iterations, 
                                   drop_threshold, iterations, converged)

        #theta0 = theta_prime
        iterations += 1
         
    #log.info(f'[genotype] EM converged after {iterations} iterations')
    
    return abundance_to_counts(eqs, theta0)
    
if __name__ == '__main__':
    
    with open('dat/info/parameters.p', 'rb') as file:
        genes, populations, databases = pickle.load(file)
    
    parser = argparse.ArgumentParser(prog='arcasHLA quant',
                                 usage='%(prog)s [options] FASTQs',
                                 add_help=False,
                                 formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        help='list of fastq files', 
                        nargs='*',
                        type=lambda x: arg_check_files(parser, x))
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('--sample',
                    help = 'sample name',
                    type = str,
                    default = None)

    parser.add_argument('--ref', 
                        type=str,
                        help='arcasHLA quant_ref path (e.g. "/path/to/ref/sample")\n  ',
                        default=None, 
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

    parser.add_argument('--pseudo', action = 'count', default=False)

    args = parser.parse_args()
    
    paired = False
    if len(args.file) == 0:
        sys.exit('[genotype] Error: FASTQ required')
    elif len(args.file) == 2:
        paired = True
    
    if args.sample == None:
        sample = os.path.basename(args.file[0]).split('.')[0]
    else:
        sample = args.sample
        
    temp = args.temp + '/' + sample
    temp, outdir = [check_path(path) for path in [temp, args.outdir]]
    
    reads_file = ''.join([temp, sample, '.reads.txt'])
    if not paired:
        num, avg, std = analyze_reads(args.file, paired, reads_file, False)
        if std == 0.0: std = .00000001

    indv_idx = args.ref + '.idx'
    indv_p = args.ref + '.p'
    indv_abundance =  outdir + sample + '.quant.tsv'
    
    
    with open(indv_p, 'rb') as file:
        genes,hla_idx,allele_idx,lengths = pickle.load(file)
    
    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        idx_allele[gene].add(idx)

    if args.pseudo:
        indv_results = outdir + sample + '.arcas.quant.json'
        results = dict()
        count_file = outdir + sample + '.counts.tsv'
        eq_file = outdir + sample + '.eq.tsv'
        
        
        command = ['kallisto pseudo', '-i', indv_idx, '-o', temp, '-t', args.threads]
        if len(args.file) == 1:
            command.extend('--single -l', str(avg), '-s', str(std))
        command.extend(args.file)
        output = run_command(command)
        
        
        n_reads = int(''.join(str(output.stderr).split('\\n')[-3].split()[2].split(',')))
        aligned_reads = int(''.join(str(output.stderr).split('\\n')[-3].split()[4].split(',')))
        
        results['total_count'] = n_reads
        results['aligned_reads'] = aligned_reads
        
        
        run_command(['mv',temp + 'pseudoalignments.ec', eq_file])
        run_command(['mv', temp + 'pseudoalignments.tsv', count_file])
        run_command(['rm -rf', temp])
        
        
        idx_allele = dict()
        for idx, alleles in allele_idx.items():
            if alleles:
                idx_allele[process_allele(sorted(alleles)[0], 2)] = idx

        
        #log.info('[alignment] processing pseudoalignment')
        # Process count information
        counts = dict()
        with open(count_file, 'r') as file:
            for line in file.read().splitlines():
                eq, count = line.split('\t')
                counts[eq] = float(count)
                
        hla_indices = set.union(*[set(indices) for indices in hla_idx.values()])
        
        idx_to_allele_idx = dict()
        for allele_id,indices in hla_idx.items():
            for idx in indices:
                idx_to_allele_idx[idx] = indices[0]
                
        # Process compatibility classes
        allele_to_eq = defaultdict(set)
        eqs = dict()
        i = 0
        with open(eq_file, 'r') as file:
            for line in file.read().splitlines():
                eq, indices = line.split('\t')
                indices_ = []
                for idx in indices.split(','):
                    if idx in hla_indices:
                        indices_.append(idx_to_allele_idx[idx])
                        allele_to_eq[allele_idx[idx]].add(i)
                    else:
                        indices_.append(idx)
                eqs[eq] = list(set(indices_))
                i += 1
                
        eq_idx = []

        for eq, indices in eqs.items():
            count = counts[eq]
            if count == 0:
                continue
                
            eq_idx.append((indices, count))
            
        counts = expectation_maximization(eq_idx, lengths, allele_idx, None, None, 
                             1e-6, 1000, 1000, 
                             0)
        
        
        scale = n_reads/1e6
        rpm = {idx:count/scale for idx, count in counts.items()}
        rpkm = {idx: idx_rpm/(lengths[idx]/1000) for idx, idx_rpm in rpm.items()}

        rpk = {idx: count/(lengths[idx]/1000) for idx, count in counts.items()}
        scale = sum(rpk.values())/1e6
        tpm = {idx: idx_rpk / scale for idx, idx_rpk in rpk.items()}

        
        for gene, allele_ids in genes.items():
            for allele_id in set(allele_ids):
                idx = hla_idx[allele_id][0]
                results[allele_id + '_count'] = counts[idx]
                results[allele_id + '_rpkm'] = rpkm[idx]
                results[allele_id + '_tpm'] = tpm[idx]
        for gene, allele_ids in genes.items():
            for allele_id in set(allele_ids):
                a1_a2 = results[gene + '1_count'] / results[gene + '2_count'] 
                a2_a1 = results[gene + '2_count'] / results[gene + '1_count'] 
                results[gene + '_ratio'] = min(a1_a2, a2_a1)
                
        with open(indv_results, 'w') as file:
            json.dump(results, file)
        
        sys.exit()
    
    indv_results = outdir + sample + '.kallisto.quant.json'
    command = ['kallisto quant', '-i', indv_idx, '-o', temp, '-t', args.threads]
    
    if len(args.file) == 1:
        command.extend('--single -l', str(avg), '-s', str(std))
        
    command.extend(args.file)
        
    run_command(command)
    run_command(['mv',temp + '/abundance.tsv',indv_abundance])
    

    kallisto_results = pd.read_csv(indv_abundance, sep = '\t')

    genes = set(genes.keys())
    
    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        if gene[:-1] in genes:
            idx_allele[gene].add(idx)

    results = defaultdict(float)
    for gene, indices in idx_allele.items():
        for idx in indices:
            results[gene + '_count'] += kallisto_results.loc[int(idx)]['est_counts']
            results[gene + '_tpm'] += kallisto_results.loc[int(idx)]['tpm']
    for gene in genes:
        c1 = results[gene + '1_count']
        c2 = results[gene + '2_count']

        if c1 == 0 or c2 == 0:
            results[gene + '_ratio'] = 0.0
        else:
            results[gene + '_ratio'] = min(c1/c2, c2/c1)
            
    with open(indv_results, 'w') as file:
        json.dump(results, file)
        
    run_command(['rm -rf', temp])
#-----------------------------------------------------------------------------
