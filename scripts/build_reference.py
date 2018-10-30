#This file is part of arcasHLA.
#
#    arcasHLA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    arcasHLA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with arcasHLA.  If not, see <https://www.gnu.org/licenses/>.


import os
import sys
import re
import json
import argparse
from subprocess import PIPE, run

from collections import defaultdict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from arcas_common import process_allele, run_command


def check_ref(verbose):
    if not os.path.isfile('database/IMGTHLA/hla.dat'):
        fetch_hla_dat(verbose)
        
    if not os.path.isfile('database/reference/hla.fasta') or not os.path.isfile('database/reference/hla.json'): 
        build_fasta(filter_reference, verbose)
        
def fetch_hla_dat(verbose):
    if os.path.isdir('database/IMGTHLA/'):
        
        run('rm -rf database/IMGTHLA', shell=True, stdout=PIPE, stderr=PIPE)
        
    command = 'git clone https://github.com/ANHIG/IMGTHLA.git database/IMGTHLA'

    run_command(command,'[reference] cloning IMGT/HLA database:',verbose)
    
def checkout_version(verbose, commithash):
    command = 'git -C database/IMGTHLA checkout ' + commithash
    
    run_command(command,'[reference] checking out IMGT/HLA:',verbose)
        
def hla_dat_version():
    results = run('git -C database/IMGTHLA show',
                  shell=True, stdout=PIPE, stderr=PIPE)
    
    lines = results.stdout.decode().split('\n')
    commithash = lines[0].split()[1]
    
    return commithash

def process_hla_dat(verbose):  
    sequences = dict()
    utrs = defaultdict(dict)
    exons = defaultdict(dict)

    sequence = partial = utr = exon = False

    complete_alleles = set()
    partial_alleles = set()
    gene_set = set()

    with open('database/IMGTHLA/hla.dat', 'r') as file:
        lines = file.read().splitlines()
        
    for line in lines:
        if line.startswith('//'):
            if sequence and allele in exons:
                gene_set.add(gene)  
                sequences[allele] = seq
                if not partial:
                    complete_alleles.add(allele)
                else:
                    partial_alleles.add(allele)
            partial = False

        elif line.startswith('FT') and 'partial' in line:
            partial = True

        elif line.startswith('FT') and re.search('allele\="HLA-', line):   
            allele = re.split('HLA-', re.sub('["\n]','',line))[1]
            gene = allele.split('*')[0]
            
            exon = sequence = False
            seq = ''

        elif line.startswith('FT') and re.search('exon',line):
            info = re.split('\s+', line)
            exon_coord = [int(info[2].split('..')[0]) - 1 , int(info[2].split('..')[1])]
            exon = True

        elif exon:
            number = re.split('"', line)[1]
            exons[allele][number] = exon_coord
            exon = False

        elif line.startswith('FT') and (re.search('\sUTR\s',line)):
            info = re.split('\s+', line)
            utr_coord = [int(info[2].split('..')[0]) - 1 , int(info[2].split('..')[1])]

            if allele not in exons:
                utrs[allele]['utr5'] = utr_coord
            else:
                utrs[allele]['utr3'] = utr_coord

        elif line.startswith('SQ'):
            sequence = True

        elif sequence and line.startswith(' '):
            seq += ''.join(line.split()[:-1]).upper()
        
    return complete_alleles, partial_alleles, gene_set, sequences, utrs, exons

def index_fasta(verbose):
    command = 'kallisto index -i database/reference/hla.idx database/reference/hla.fasta'
    run_command(command, '[reference] indexing reference with Kallisto:',verbose)
    
def write_complete(sequences_out,sequence_info,verbose):
    with open('database/reference/hla.fasta','w') as file:
        SeqIO.write(sequences_out, file, 'fasta')
        
    commithash = hla_dat_version()
    with open('database/reference/hla.json','w') as file:
        json.dump([commithash,sequence_info],file)
        
    command = 'kallisto index -i database/reference/hla.idx database/reference/hla.fasta'
    run_command(command, '[reference] indexing reference with Kallisto:',verbose)

def write_partial(sequences_out,sequence_info,verbose):
    with open('database/reference/hla.partial.fasta','w') as file:
        SeqIO.write(sequences_out, file, 'fasta')

    commithash = hla_dat_version()
    with open('database/reference/hla.partial.json','w') as file:
        json.dump([commithash,sequence_info],file)
        
    command = 'kallisto index -i database/reference/hla.partial.idx database/reference/hla.partial.fasta'
    run_command(command, '[reference] indexing partial reference with Kallisto:',verbose)

def build_fasta(verbose):
    
    if verbose: print('[reference] processing IMGT/HLA database')
    
    def build_complete_sequence(allele):
        gene = allele.split('*')[0]
        coords = [[start,stop] for n,(start,stop) in sorted(exons[allele].items())]
        sequence = ''.join([sequences[allele][start:stop] for start,stop in coords])
        
        if allele in stop_loss_alleles:
            sequence = sequence[:stop_loss_alleles[allele]]

        
        mRNA_sequences[sequence].add(allele)
        
        if allele in utrs:
            for (start,stop) in utrs[allele].values():
                sequence = sequences[allele][start:stop]
                other_sequences.add(sequence)
                
    def build_combination_sequences(allele):
        for exon_group in exon_combinations:
            if set(exon_group) & set(exons[allele]) != set(exon_group):
                continue
 
            coords = []
            for n in exon_group: coords.append(exons[allele][n])
            sequence = ''.join([sequences[allele][start:stop] for start,stop in coords])
            combo_sequences[str(exon_group)][sequence].add(allele)

    def add_complete_records(mRNA_sequences, other_sequences):
        sequences_out = []
        allele_index = dict()
        allele_lengths = dict()
    
        for index,(sequence,alleles) in enumerate(sorted(mRNA_sequences.items(),key=lambda x:x[1])):
            length = len(sequence)
            
            sequences_out.append(SeqRecord(Seq(sequence),id=str(index),description=''))
            
            allele_index[index] = sorted(alleles)
            allele_lengths[index] = length
        
        offset = index + 1
        for index, sequence in enumerate(other_sequences):
            sequences_out.append(SeqRecord(Seq(sequence),id=str(index + offset),description=''))
            allele_index[index + offset] = None
            
        return sequences_out, allele_index, allele_lengths
    
    def add_partial_records(combo_sequences, other_sequences):
        sequences_out = []
        exon_index = dict()
        allele_index = dict()
        allele_lengths = dict()
    
        offset = 0
        for exon in sorted(combo_sequences):
            for index,(sequence,alleles) in enumerate(sorted(combo_sequences[exon].items(),key=lambda x:x[1])):
                length = len(sequence)

                sequences_out.append(SeqRecord(Seq(sequence),id=str(index + offset),description=''))

                allele_index[index + offset] = sorted(alleles)
                allele_lengths[index + offset] = length
                exon_index[index + offset] = exon
            offset += index + 1
            
        for index, sequence in enumerate(other_sequences):
            sequences_out.append(SeqRecord(Seq(sequence),id=str(index + offset),description=''))
            allele_index[index + offset] = None
        
        return sequences_out, allele_index, allele_lengths, exon_index

    
    complete_alleles, partial_alleles, gene_set, sequences, utrs, exons = process_hla_dat(verbose)
    
    mRNA_sequences = defaultdict(set)
    combo_sequences = {str(i):defaultdict(set) for i in exon_combinations}
    other_sequences = set()

    for allele in complete_alleles:
        build_complete_sequence(allele)
        build_combination_sequences(allele)
        
    for allele in partial_alleles:
        build_complete_sequence(allele)
        build_combination_sequences(allele)
        
    mRNA_sequences = {seq:sorted(alleles) for seq, alleles in mRNA_sequences.items()}
    other_sequences = list(other_sequences)

    if verbose: print('[reference] building HLA database')
    genes = sorted(gene_set)
    sequences_out, allele_index, allele_lengths = add_complete_records(mRNA_sequences, other_sequences)
    write_complete(sequences_out,[genes, allele_index, allele_lengths],verbose)
    
    if verbose: print('[reference] building partial HLA database')
    sequences_out, allele_index, allele_lengths, exon_index = add_partial_records(combo_sequences, other_sequences)
    partial_alleles = {process_allele(allele,3) for allele in partial_alleles}
    partial_exons = {process_allele(allele,3):exons[allele] for allele in partial_alleles}
    partial_index = {index:[process_allele(allele,3) for allele in alleles] if alleles else None for index,alleles in allele_index.items()}
    write_partial(sequences_out,[genes, allele_index, allele_lengths, partial_index, exon_index, partial_exons, list(partial_alleles)],verbose)

stop_loss_alleles = {'C*04:09N':1100}
exon_combinations = [['2'],['2','3'],['1','2'],['1','2','3'],['2','3','4'],['1','2','3','4'],['2','3','4','5'],['1','2','3','4','5'],['2','3','4','5','6'],['1','2','3','4','5','6'],['2','3','4','5','6','7'],['1','2','3','4','5','6','7']]

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--update', action = 'count', help='force update', default=False)
    parser.add_argument('--version', type = str, help='checkout IMGT/HLA version using commithash', default=False)
    parser.add_argument('-v','--verbose', action = 'count',default=False)
    
    args = parser.parse_args()

    if args.update:
        if args.verbose: print('[reference] updating HLA reference')
        #fetch_hla_dat(args.verbose)
        build_fasta(args.verbose)
        
    if args.version:
        checkout_version(args.verbose, args.version)
        build_fasta(args.verbose)

    check_ref(args.verbose)
        
    commithash = hla_dat_version()
    if args.verbose: print('[reference] HLA database version:',commithash)
        
