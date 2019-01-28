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
    indv_results = outdir + sample + '.quant.json'

    if args.pseudo:
        command = ['kallisto pseudo', '-i', indv_idx, '-o', temp, '-t', args.threads]
        if len(args.file) == 1:
            command.extend('--single -l', str(avg), '-s', str(std))
        command.extend(args.file)
        run_command(command)
        run_command(['mv',temp + 'pseudoalignments.ec', outdir + sample + '.classes.tsv'])
        run_command(['mv', temp + 'pseudoalignments.tsv', outdir + sample + '.counts.tsv'])
        run_command(['rm -rf', temp])
        sys.exit()
    
    
    command = ['kallisto quant', '-i', indv_idx, '-o', temp, '-t', args.threads]
    
    if len(args.file) == 1:
        command.extend('--single -l', str(avg), '-s', str(std))
        
    command.extend(args.file)
        
    run_command(command)
    run_command(['mv',temp + '/abundance.tsv',indv_abundance])
    
    
    with open(indv_p, 'rb') as file:
        genes,hla_idx,allele_idx,lengths = pickle.load(file)
    
    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        idx_allele[gene].add(idx)

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
