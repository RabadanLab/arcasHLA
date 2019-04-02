#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   merge.py: merges genotype files into a single json or table.
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

from argparse import RawTextHelpFormatter

from arcas_utilities import check_path

#-------------------------------------------------------------------------------

__version__     = '0.2'
__date__        = '2019-04-02'

#-------------------------------------------------------------------------------

def get_paths(indir):
    '''Get all file paths that match arcasHLA output.'''
    partial_files = []
    genotype_files = []
    
    log_files = []
    quant_files = []
    
    for file in os.listdir(indir):
        if file.endswith('.partial_genotype.json'):
            partial_files.append(file)
        elif file.endswith('.genotype.json'):
            genotype_files.append(file)
        elif file.endswith('.quant.json'):
            quant_files.append(file)
        elif file.endswith('.genotype.log'):
            log_files.append(file)
            
    return genotype_files, partial_files, log_files, quant_files
    
def process_json(json_files, indir, outdir, run, suffix):
    '''Merge genotype.json files.'''
    file_out = [outdir]
    if run: file_out.extend([run,'.'])
    file_out.append(suffix)
    
    genotypes = dict()
    for file in json_files:
        sample = file.split('.')[0]
        file_path = indir + file

        with open(file_path,'r') as file:
            genotypes[sample] = json.load(file)
            
    with open(''.join(file_out), 'w') as file:
        json.dump(genotypes, file)
        
def process_count(count_files, indir, outdir, run, suffix):
    file_out = [outdir]
    if run: file_out.extend([run,'.'])
    file_out.append(suffix)
    file_out = ''.join(file_out)
    
    all_counts = []
    for file in count_files:
        sample = file.split('.')[0]
        file_path = indir + file

        with open(file_path, 'r') as file:
            lines = file.read()

        lines = lines.split('-'*80)[2].split('\n')
        counts = {'Sample':sample}
        for line in lines:
            if line.startswith('[alignment] Processed '):
                _,_,counts['total_count'],_,counts['aligned_count'],_,_,_,_ = line.split()
            elif line.endswith(' reads mapped to a single HLA gene'):
                counts['single_count'] = line.split()[1]
            elif line.startswith('\t\tHLA-'):
                line = line.split()
                gene = line[0].split('-')[1]
                counts[gene + '_read_count'] = line[2]
        all_counts.append(counts)
            
    pd.DataFrame(all_counts).set_index('Sample').to_csv(file_out, sep = '\t')
        
def process_quant(json_files, indir, outdir, run, suffix):
    file_out = [outdir]
    if run: file_out.extend([run,'.'])
    file_out.append(suffix)
    file_out = ''.join(file_out)
    
    all_results = []
    for file in json_files:
        sample = file.split('.')[0]
        file_path = indir + file

        with open(file_path,'r') as file:
            results = json.load(file)
            results['Sample'] = sample
            
            all_results.append(results)
            
    pd.DataFrame(all_results).set_index('Sample').to_csv(file_out, sep = '\t')
    

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
                        '--input',
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
    
    outdir = check_path(args.outdir)
    
    if os.path.isdir(args.input):
        
        # add convert polysolver output
        
        indir = check_path(args.input)
    
        genotype_files, partial_files, count_files, quant_files = get_paths(indir)
    
        if genotype_files:
            process_json(genotype_files, 
                         indir, 
                         outdir, 
                         args.run, 
                         'genotypes.json')

        if partial_files:
            process_json(partial_files, 
                         indir, 
                         outdir,  
                         args.run, 
                         'partial_genotypes.json')
        if count_files:
            process_count(count_files,
                         indir,
                         outdir,
                         args.run,
                         'counts.tsv')

        if quant_files:
            process_quant(quant_files, 
                         indir, 
                         outdir,  
                         args.run, 
                         'quant.tsv')
    else:
        # check input type: current, pre-2009 (only if Cw), polysolver
        # check input resolution: number of fields
        # warn if number of fields is inconsistent
        
        # output type is default: output current nomenclature, same number of fields as input
        # change resolution: warn if output resolution > input resolution (take the first allele)
        # change output type:
        #    polysolver - increase resolutions with warning, warn if allele is missing, 
        #                 take nearest neighbor by ambiguity then by numerical order
        #                 option to output to individual files (input for LOHHLA)
        #    ambiguity - output with "g" ending
        #    super-groups - output based on supergrouping
        pass
        
#-------------------------------------------------------------------------------
