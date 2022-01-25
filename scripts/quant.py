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
rootDir = os.path.dirname(os.path.realpath(__file__)) + '/../'
parameters = rootDir + 'dat/info/parameters.p'
#-------------------------------------------------------------------------------
def arg_check_files(parser, arg):
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error('The file %s does not exist.' %file)
        elif not (file.endswith('alignment.p') or file.endswith('.fq.gz') or file.endswith('.fastq.gz') or file.endswith('.tsv') or file.endswith('.json')):
            parser.error('The format of %s is invalid.' %file)
        return arg

def analyze_reads(fqs, paired, reads_file):
    '''Analyzes read length for single-end sampled, required by Kallisto.'''
    
    awk = "| awk '{if(NR%4==2) print length($1)}'"
    
    if fqs[0].endswith('.gz'):
        cat = 'zcat'
    else:
        cat = 'cat'
    
    log.info('[alignment] Analyzing read length')
    if paired:
        fq1, fq2 = fqs
        
        command = [cat, '<', fq1, awk, '>' , reads_file]
        run_command(command)
        
        command = [cat, '<', fq2, awk, '>>', reads_file]
        run_command(command)
        
    else:
        fq = fqs[0]
        command = [cat, '<', fq, awk, '>', reads_file]
        run_command(command)
        
    read_lengths = np.genfromtxt(reads_file)
    
    if len(read_lengths) == 0:
        sys.exit('[genotype] Error: FASTQ files are empty; check arcasHLA extract for issues.')
    
    num = len(read_lengths)
    avg = round(np.mean(read_lengths), 6)
    std = round(np.std(read_lengths), 6)
    
    
    return num, avg, std
    
if __name__ == '__main__':
    
    with open(parameters, 'rb') as file:
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
        
    indv_idx = args.ref + '.idx'
    indv_p = args.ref + '.p'
    indv_abundance =  outdir + sample + '.quant.tsv'
    allele_results_json = outdir + sample + '.quant.alleles.json'
    gene_results_json = outdir + sample + '.quant.genes.json'
    allele_results_tsv = outdir + sample + '.quant.alleles.tsv'
    gene_results_tsv = outdir + sample + '.quant.genes.tsv'


    with open(indv_p, 'rb') as file:
        genes,genotype,hla_idx,allele_idx,lengths = pickle.load(file)

    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        idx_allele[gene].add(idx)
        
    if args.file[0].endswith('.fq.gz'):
    

        reads_file = ''.join([temp, sample, '.reads.txt'])
        if not paired:
            num, avg, std = analyze_reads(args.file, paired, reads_file)
            if std == 0.0: std = .00000001


        command = ['kallisto quant', '-i', indv_idx, '-o', temp, '-t', args.threads]

        if len(args.file) == 1:
            command.extend(['--single -l', str(avg), '-s', str(std)])

        command.extend(args.file)

        output = run_command(command).stderr.decode()

        total_reads = re.findall('(?<=processed ).+(?= reads,)',output)[0]
        total_reads = int(re.sub(',','',total_reads))
        aligned_reads = re.findall('(?<=reads, ).+(?= reads pseudoaligned)',output)[0]
        aligned_reads = int(re.sub(',','',aligned_reads))


        run_command(['mv',temp + '/abundance.tsv',indv_abundance])
        kallisto_results = pd.read_csv(indv_abundance, sep = '\t')
        
    else:
        with open(args.file[1], 'r') as file:
            previous_results = json.load(file)
            
        #total_reads = previous_results['total_count']
        #aligned_reads = previous_results['aligned_reads']
            
        kallisto_results = pd.read_csv(args.file[0], sep = '\t')

    

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
    '''
    scale = aligned_reads/1e6
    rpm = {idx:count/scale for idx, count in counts.items()}
    rpkm = {idx: idx_rpm/(lengths[idx]/1000) for idx, idx_rpm in rpm.items()}

    
    gene_results = defaultdict(float)
    gene_results['total_count'] = total_reads
    gene_results['aligned_reads'] = aligned_reads
    
    allele_results = defaultdict(float)
    allele_results['total_count'] = total_reads
    allele_results['aligned_reads'] = aligned_reads
    
    for allele_id, allele in genotype.items():
        allele_results[allele_id] = allele

    total_hla_count = 0.
    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            gene_results[gene + '_count'] += counts[allele_id]
            total_hla_count += counts[allele_id]
            #gene_results[gene + '_rpkm'] += rpkm[allele_id]
            gene_results[gene + '_tpm'] += tpm[allele_id]
            
            allele_results[allele_id + '_count'] = counts[allele_id]
            #allele_results[allele_id + '_rpkm'] = rpkm[allele_id]
            allele_results[allele_id + '_tpm'] = tpm[allele_id]
    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            baf = allele_results[gene + '1_count'] / (allele_results[gene + '1_count'] + allele_results[gene + '2_count'])
            allele_results[gene + '_baf'] = min(baf, 1-baf)
    '''
    gene_results = {gene:defaultdict(int) for gene in genes}
    
    allele_results = {gene:defaultdict(float) for gene in genes}
    
    total_hla_count = 0
    for allele_id, allele in genotype.items():
        allele_results[allele_id[:-1]]['allele' + allele_id[-1]] = allele
        total_hla_count += counts[allele_id]

    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            gene_results[gene]['count'] += round(counts[allele_id])
            gene_results[gene]['tpm'] += round(tpm[allele_id])
            if counts[allele_id]:
                gene_results[gene]['abundance'] += counts[allele_id]/total_hla_count
            
            allele_results[gene]['allele' + allele_id[-1] + '_count'] = round(counts[allele_id])
            allele_results[gene]['allele' + allele_id[-1] + '_tpm'] = round(tpm[allele_id])
    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            baf = allele_results[gene]['allele1_count'] / (allele_results[gene]['allele1_count'] + allele_results[gene]['allele2_count'])
            allele_results[gene]['baf'] = round(min(baf, 1-baf),2)
    for gene in genes:
        gene_results[gene]['abundance'] = str(round(gene_results[gene]['abundance']*100,2)) + '%'
    '''
    print('-'*80)
    print('[quant] Observed HLA genes:')

    print('        {: <10}    {}    {}    {}'
             .format('gene','count','   TPM',' abundance'))
    for gene in sorted(genes):
        print('        HLA-{: <6}    {: >5.0f}    {: >6.0f}    {: >9.2f}%'
                 .format(gene, gene_results[gene + '_count'], gene_results[gene + '_tpm'], (gene_results[gene + '_count']/total_hla_count)*100))
            
    print('-'*80)
    
    
    print('-'*80)
    print('[quant] Observed HLA genes:')

    print('        {: <10}    {: <12}    {}    {: <12}    {}'
          .format('gene','allele 1','count','allele 2','count'))
    
    for gene in sorted(genes):
        a1 = allele_results[gene + '1']
        a2 = allele_results[gene + '2']
        c1 = allele_results[gene+'1_count']
        c2 = allele_results[gene+'2_count']
        if c2 == 0: baf = 0
        else: baf = min(c1/(c1+c2),c2/(c1+c2))
        if not a2: a2 = ''
        print('        HLA-{: <6}    {: <12}    {: >5.0f}    {: <12}    {: >5.0f}     {: > 5.2f}'
                 .format(gene, a1, c1, a2, c2, baf))
    
    '''
    
    df = pd.DataFrame(allele_results).T
    df.index.names = ['gene']
    df = df[['allele1','allele2', 'allele1_count',  'allele2_count', 'allele2_tpm','allele1_tpm', 'baf']]
    df.to_csv(allele_results_tsv,sep='\t')
    
    df = pd.DataFrame(gene_results).T
    df.index.names = ['gene']
    df = df[['count','tpm','abundance']]
    df.to_csv(gene_results_tsv,sep='\t')
    
    with open(allele_results_json, 'w') as file:
        json.dump(allele_results, file)
        
    with open(gene_results_json, 'w') as file:
        json.dump(gene_results, file)
        
    if not args.keep_files: run_command(['rm -rf', temp])
#-----------------------------------------------------------------------------
