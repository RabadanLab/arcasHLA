#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   quant_ref.py: genotypes from extracted chromosome 6 reads.
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

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from reference import check_ref
from arcas_utilities import *

#-------------------------------------------------------------------------------

__version__     = '0.2'
__date__        = '2019-04-02'

#-------------------------------------------------------------------------------
#   Paths and filenames
#-------------------------------------------------------------------------------

rootDir            = os.path.dirname(os.path.realpath(__file__)) + '/../'
cDNA_p             = rootDir + 'dat/ref/cDNA.p'
cDNA_single_p      = rootDir + 'dat/ref/cDNA.single.p'
GRCh38_chr6        = rootDir + 'dat/ref/GRCh38.chr6.noHLA.fasta'
GRCh38             = rootDir + 'dat/ref/GRCh38.all.noHLA.fasta'
HLA_json           = rootDir + 'dat/ref/hla_transcripts.json'
dummy_HLA_fa       = rootDir + 'dat/ref/GRCh38.chr6.HLA.fasta'

#-------------------------------------------------------------------------------

def build_custom_reference(subject, genotype, mode, chr6, temp):
    
    dummy_HLA_dict = SeqIO.to_dict(SeqIO.parse(dummy_HLA_fa, 'fasta')) 
    
    if chr6:
        transcriptome = list(SeqIO.parse(GRCh38_chr6, 'fasta'))
    else:
        transcriptome = list(SeqIO.parse(GRCh38, 'fasta'))
        
    with open(HLA_json,'r') as file:
        HLA_transcripts = json.load(file)
    
    genes = {allele_id[:-1] for allele_id in genotype.keys()}
    for gene in set(HLA_transcripts) - genes:
        for transcript in HLA_transcripts[gene]:
            transcriptome.append(dummy_HLA_dict[transcript])
     
    with open('dat/ref/allele_groups.p','rb') as file:
        groups = pickle.load(file)
    with open(cDNA_p,'rb') as file:
        cDNA = pickle.load(file)    
    with open(cDNA_single_p,'rb') as file:
        cDNA_single = pickle.load(file)
    
    indv_fasta = ''.join([temp,subject,'.fasta'])
    indv_idx  = ''.join([outdir,subject,'.idx'])
    indv_p = ''.join([outdir,subject,'.p'])

    indv_records = []

    allele_idx = dict()
    lengths = dict()
    genes = defaultdict(list)
    hla_idx = defaultdict(list)

    idx = 0
    for allele_id, allele in genotype.items():
        gene = get_gene(allele)
        
        if mode == 'single':
            sequences = [cDNA_single[allele]]
        elif mode == 'G-group':
            sequences = [seq for a in groups[allele] for seq in cDNA[a]]
        else:
            sequences = [seq for seq in cDNA[allele]]
            
        for seq in sequences:
            hla_idx[allele_id].append(str(idx))
            genes[gene].append(allele_id)
            allele_idx[str(idx)] = allele_id
            lengths[str(idx)] = len(seq)

            record = SeqRecord(Seq(seq),
                               id=str(idx),
                               description='')

            indv_records.append(record)
            idx += 1
     
    for transcript in transcriptome:
        allele_idx[str(idx)] = transcript.id
        lengths[str(idx)] = len(seq)

        record = SeqRecord(transcript.seq,
                           id=str(idx),
                           description='')

        indv_records.append(record)
        idx += 1

    SeqIO.write(indv_records, indv_fasta, 'fasta')

    with open(indv_p, 'wb') as file:
        pickle.dump([genes,genotype,hla_idx,allele_idx,lengths], file)
        

    run_command(['kallisto', 'index','-i', indv_idx, indv_fasta])
    
def process_json_genotype(input_genotype):
    genotype = dict()
    for gene, alleles in input_genotype.items():
        alleles = [process_allele(allele,2) for allele in alleles]
        if len(alleles) == 2:
            genotype[gene + '1'], genotype[gene + '2'] = alleles
        else:
            genotype[gene + '1'], genotype[gene + '2'] = alleles[0]
    return genotype

def process_str_genotype(input_genotype):
    genotype = dict()
    for allele in input_genotype.split(','):
        gene = get_gene(allele)
        if gene + '1' not in genotype:
            genotype[gene + '1'] = process_allele(allele,2)
        elif gene + '2' not in genotype:
            genotype[gene + '1'] = process_allele(allele,2)
        else:
            sys.exit('[quant] Error: more than 2 alleles provided for a gene.')
            run_command(['rm -rf', temp])
            
    return genotype
    
    
if __name__ == '__main__':
    
    with open('dat/info/parameters.p', 'rb') as file:
        genes, populations, databases = pickle.load(file)
    
    parser = argparse.ArgumentParser(prog='arcasHLA quant',
                                 usage='%(prog)s [options] FASTQs',
                                 add_help=False,
                                 formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('input',
                        type=str,
                        nargs='*',
                        help='arcasHLA output genotype.json or genotypes.json \nor tsv with format specified in README.md\n\n')
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('--chr6', 
                        action='count',
                        help='restrict transcriptome to chr 6\n  default: False\n\n',
                        default=False)

    parser.add_argument('--resolution', 
                        type = int,
                        help='genotype resolution, only use >2 when typing performed with assay or Sanger sequencing\n  default: 2\n\n',
                        default='2')
    
    parser.add_argument('--mode', 
                        type=str,
                        help='number of transcripts to include per allele\n single: one 3-field resolution transcript per allele (e.g. A*01:01:01)\nG-group (all transcripts with identical binding regions)\n  default: protein-group\n\n',
                        default='protein-group')
    
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
    
    if args.resolution != 2:
        sys.exit('[quant] only 2-field resolution supported at this time.')
        
    outdir = check_path(args.outdir)
    temp = args.temp
        
    if args.input[0].endswith('genotype.json'):
        subject = os.path.basename(args.input[0]).split('.')[0]
        if args.verbose: print('[quant] Building reference for', subject)
        
        with open(args.input[0], 'r') as file:
            input_genotype = json.load(file)
            
        genotype = process_json_genotype(input_genotype)
        
        build_custom_reference(subject, genotype, args.mode, args.chr6, temp)
        
    elif len(args.input) == 2:
        subject, filepath = args.input
        if args.verbose: print('[quant] Building reference for', subject)
        if filepath.endswith('.json'):
            with open(filepath, 'r') as file:
                input_genotype = json.load(file)[subject]
            genotype = process_json_genotype(input_genotype)
        else:
            input_genotypes = pd.read_csv(args.input[1], sep='\t').set_index('subject').to_dict('index')
            genotype = input_genotypes[subject]
            
        build_custom_reference(subject, genotype, args.mode, args.chr6, temp)
        
                
    elif args.input[0].endswith('genotypes.json') or args.input[0].endswith('.tsv'):
        temp = create_temp(temp)
        
        if args.verbose: print('[quant] Building references from',os.path.basename(args.input[0]))
        
        if args.input[0].endswith('genotypes.json'):
            with open(args.input[0], 'r') as file:
                subjects = json.load(file).keys()   
        else:
            subjects = pd.read_csv(args.input[0], sep='\t')['subject']
            
            
        subject_file = temp + 'subjects.txt'
        with open(subject_file,'w') as file:
            file.write('\n'.join(subjects))
            
        command = ['cat', subject_file, '|','parallel', '-j', args.threads,
                   rootDir + '/arcasHLA', 'quant_ref', 
                   '{}', args.input[0],
                   '--resolution', args.resolution,
                   '--outdir', outdir,
                   '--temp', temp]
        
        if args.chr6: command.append('--chr6')
        if args.mode: command.extend(['--mode', args.mode])
        if args.verbose: command.append('--verbose')
        
        run_command(command)
            
    
        if not args.keep_files: run_command(['rm -rf', temp])

#-----------------------------------------------------------------------------
