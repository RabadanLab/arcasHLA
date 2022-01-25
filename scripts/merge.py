#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   merge.py: merges genotype files into a single table.
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
import json
import argparse
import pandas as pd

from collections import defaultdict
from argparse import RawTextHelpFormatter

from arcas_utilities import check_path


#-------------------------------------------------------------------------------

__version__     = '0.3.0'
__date__        = '2021-12-09'

#-------------------------------------------------------------------------------

def get_paths(indir):
    '''Get all file paths that match arcasHLA output.'''
    partial_files = []
    genotype_files = []
    gene_count_files = []
    quant_files = []
    
    for file in os.listdir(indir):
        if file.endswith('.partial_genotype.json'):
            partial_files.append(file)
        elif file.endswith('.genotype.json'):
            genotype_files.append(file)
        elif file.endswith('.genes.json'):
            gene_count_files.append(file)
        elif file.endswith('.alleles.json'):
            quant_files.append(file)
            
    return genotype_files, partial_files, gene_count_files, quant_files
    
def process_json(json_files, indir, outdir, run, suffix):
    '''Merge genotype.json files into single tsv.'''
    file_out = ''.join([outdir, run, suffix])
    
    genotypes = dict()
    for file in json_files:
        sample = file.split('.')[0]
        file_path = indir + file

        with open(file_path,'r') as file:
            genotypes[sample] = json.load(file)

    #with open(file_out + '.json', 'w') as file:
    #    json.dump(genotypes, file)
        
    genotypes_ = defaultdict(dict)
    for sample, genotype in genotypes.items():
        for gene, alleles in genotype.items():
            if len(alleles) == 2:
                genotypes_[sample][gene + '1'] = alleles[0]
                genotypes_[sample][gene + '2'] = alleles[1]
            else:
                genotypes_[sample][gene + '1'] = alleles[0]
                genotypes_[sample][gene + '2'] = alleles[0]

    pd.DataFrame(genotypes_).T.rename_axis('subject').to_csv(file_out + '.tsv', sep = '\t')
    

def process_quant(json_files, indir, outdir, run, suffix):
    file_out = [outdir]
    if run: file_out.append(run)
    file_out = ''.join(file_out)

    all_results = []
    for file in json_files:
        sample = file.split('.')[0]
        file_path = indir + file

        results = {'sample':sample}
        with open(file_path,'r') as file:
            results_ = json.load(file)
            for gene, gene_results in results_.items():
                for key,value in gene_results.items():
                    results[gene + '_' + key] = value

                c1 = gene_results['allele1_count']
                c2 = gene_results['allele2_count']
                if c1 == 0 or c2 == 0:
                    results[gene + '_baf'] = 0
                else:
                    results[gene + '_baf'] = min(c1/(c1+c2),c2/(c1+c2))

            all_results.append(results)

    pd.DataFrame(all_results).set_index('sample').to_csv(file_out + '.tsv', sep = '\t')
    with open(file_out + '.json', 'w') as file:
        json.dump(all_results, file)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA merge',
                                     usage='%(prog)s [options]',
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)

    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)

    parser.add_argument('-i',
                        '--indir',
                        type=str,
                        help='directory containing arcasHLA files\n\n',
                        default='.',
                        metavar='')

    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='out directory\n\n',
                        default='.',
                        metavar='')

    parser.add_argument('--run',
                        type=str,
                        help='run name\n\n',
                        default='',
                        metavar='')


    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    if args.run: args.run += '.'
    
    indir, outdir = [check_path(path) for path in [args.indir, args.outdir]]
    
    (genotype_files, partial_files,
            gene_count_files, quant_files) = get_paths(indir)

    if genotype_files:
        process_json(genotype_files, indir, outdir,
                args.run, 'genotypes')
        
    if partial_files:
        process_json(partial_files, indir, outdir,
                args.run, 'partial_genotypes')
        
    if gene_count_files:
        process_json(gene_count_files, indir, outdir,
                args.run, 'genes')

    if quant_files:
        process_quant(quant_files, indir, outdir,
                args.run, 'quant')
        
#-------------------------------------------------------------------------------
