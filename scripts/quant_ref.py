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

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

cDNA            = 'dat/ref/cDNA.p'
GRCh38_chr6     = 'dat/ref/GRCh38.chr6.noHLA.fasta'
GRCh38          = 'dat/ref/GRCh38.all.noHLA.fasta'
HLA_transcripts = 'dat/ref/hla_transcripts.json'
dummy_HLA       = 'dat/ref/GRCh38.chr6.HLA.fasta'

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

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from reference import check_ref
from arcas_utilities import (process_allele, check_path, remove_files, 
                            get_gene, hline)

from subprocess import PIPE, run

__version__     = '1.0'
__date__        = 'November 2018'

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

def run_command(command, message = ''):
    '''Outputs message and command to log, runs command and returns output.'''
    if type(command) == list:
        command = ' '.join(command)

    ##if message: log.info(''.join([message,'\n\n\t', command,'\n']))
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
        
def arg_check_genotype(parser, arg):
    input_genotype = defaultdict(list)
    for allele in arg.split(','):
        gene = get_gene(allele)
        if gene not in genes:
            parser.error('The gene {} is invalid.', gene)
        input_genotype[gene].append(allele)
        
    genotype = dict()
    for gene, alleles in input_genotype.items():
        if len(alleles) == 1:
            genotype[gene + '1'] = genotype[gene + '2'] = alleles[0]
        elif len(alleles) == 2:
            genotype[gene + '1'], genotype[gene + '2'] = alleles
        else:
            parser.error('A gene can have a max of 2 alleles')
    return genotype

    
if __name__ == '__main__':
    
    with open('dat/info/parameters.p', 'rb') as file:
        genes, populations, databases = pickle.load(file)
    
    parser = argparse.ArgumentParser(prog='arcasHLA quant',
                                 usage='%(prog)s [options] FASTQs',
                                 add_help=False,
                                 formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('individual', 
                        help='individual ID\n',
                        type=str)
    
    parser.add_argument('genotype', 
                        help='comma separated list of HLA alleles\n',
                        type=lambda x: arg_check_genotype(parser, x))
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('--chr6', 
                        action='count',
                        help='restrict transcriptome to chr 6\n  default: False\n\n',
                        default=False)

    
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
        
    temp, outdir = [check_path(path) for path in [args.temp, args.outdir]]

    dummy_HLA_dict = SeqIO.to_dict(SeqIO.parse(dummy_HLA, 'fasta'))
    
    if args.chr6:
        transcriptome = list(SeqIO.parse(GRCh38_chr6, 'fasta'))
    else:
        transcriptome = list(SeqIO.parse(GRCh38, 'fasta'))
        
    with open(HLA_transcripts,'r') as file:
        HLA_transcripts = json.load(file)
    
    for gene in set(HLA_transcripts) - set(args.genotype.keys()):
        for transcript in HLA_transcripts[gene]:
            transcriptome.append(dummy_HLA_dict[transcript])
            

    with open(cDNA,'rb') as file:
        cDNA = pickle.load(file)
        
    ID = args.individual

    genotype = args.genotype
    
    indv_fasta = ''.join([outdir,ID,'.fasta'])
    indv_idx  = ''.join([outdir,ID,'.idx'])
    indv_p = ''.join([outdir,ID,'.p'])

    indv_records = []

    allele_idx = dict()
    lengths = defaultdict(list)
    id_to_gene = dict()

    idx = 0
    for allele_id, allele in genotype.items():
        gene = get_gene(allele)
        for seq in cDNA[allele]:
            id_to_gene[allele_id] = gene
            allele_idx[str(idx)] = allele_id
            lengths[allele_id].append(len(seq))

            record = SeqRecord(Seq(seq),
                               id=str(idx),
                               description='')

            indv_records.append(record)
            idx += 1
            
    for transcript in transcriptome:
        allele_idx[str(idx)] = transcript.id
        lengths[transcript.id].append(len(seq))

        record = SeqRecord(transcript.seq,
                           id=str(idx),
                           description='')

        indv_records.append(record)
        idx += 1

    SeqIO.write(indv_records, indv_fasta, 'fasta')

    with open(indv_p, 'wb') as file:
        pickle.dump([id_to_gene,allele_idx,lengths], file)

    run_command(['kallisto', 'index','-i', indv_idx, indv_fasta])

#-----------------------------------------------------------------------------
