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
from arcas_utilities import *

#-------------------------------------------------------------------------------

__version__     = '0.2'
__date__        = '2019-04-02'

#-------------------------------------------------------------------------------

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
        
    outdir = check_path(args.outdir)
    temp = create_temp(args.temp)
    
    reads_file = ''.join([temp, sample, '.reads.txt'])
    if not paired:
        num, avg, std = analyze_reads(args.file, paired, reads_file, False)
        if std == 0.0: std = .00000001

    indv_idx = args.ref + '.idx'
    indv_p = args.ref + '.p'
    indv_abundance =  outdir + sample + '.quant.tsv'
    
    
    with open(indv_p, 'rb') as file:
        genes,genotype,hla_idx,allele_idx,lengths = pickle.load(file)
    
    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        idx_allele[gene].add(idx)
    
    indv_results = outdir + sample + '.quant.json'
    command = ['kallisto quant', '-i', indv_idx, '-o', temp, '-t', args.threads]
    
    if len(args.file) == 1:
        command.extend('--single -l', str(avg), '-s', str(std))
        
    command.extend(args.file)
        
    output = run_command(command).stderr.decode()

    total_reads = re.findall('(?<=processed ).+(?= reads,)',output)[0]
    total_reads = int(re.sub(',','',total_reads))
    aligned_reads = re.findall('(?<=reads, ).+(?= reads pseudoaligned)',output)[0]
    aligned_reads = int(re.sub(',','',aligned_reads))
    
    
    run_command(['mv',temp + '/abundance.tsv',indv_abundance])
    

    kallisto_results = pd.read_csv(indv_abundance, sep = '\t')

    idx_allele = defaultdict(set)
    hla_indices = set()
    for idx, gene in allele_idx.items():
        if gene[:-1] in genes:
            idx_allele[gene].add(idx)
            hla_indices.add(int(idx))

    lengths = defaultdict(float)
    counts = defaultdict(float)
    tpm = defaultdict(float)
    for gene, indices in idx_allele.items():
        for idx in indices:
            counts[gene] += kallisto_results.loc[int(idx)]['est_counts']
            lengths[gene] += kallisto_results.loc[int(idx)]['length']
            tpm[gene] += kallisto_results.loc[int(idx)]['tpm']

    scale = aligned_reads/1e6
    rpm = {idx:count/scale for idx, count in counts.items()}
    rpkm = {idx: idx_rpm/(lengths[idx]/1000) for idx, idx_rpm in rpm.items()}

    results = defaultdict(float)
    results['total_count'] = total_reads
    results['aligned_reads'] = aligned_reads

    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            results[allele_id + '_count'] = counts[allele_id]
            results[allele_id + '_rpkm'] = rpkm[allele_id]
            results[allele_id + '_tpm'] = tpm[allele_id]
    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            baf = results[gene + '1_count'] / (results[gene + '1_count'] + results[gene + '2_count'])
            results[gene + '_baf'] = min(baf, 1-baf)
            
    with open(indv_results, 'w') as file:
        json.dump(results, file)
        
    if not args.keep_files: run_command(['rm -rf', temp])
#-----------------------------------------------------------------------------
