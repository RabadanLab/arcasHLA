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

from argparse import RawTextHelpFormatter

from arcas_utilities import check_path

__version__     = '0.1.1'
__date__        = '2019-05-07'

#-------------------------------------------------------------------------------

def get_paths(indir):
    '''Get all file paths that match arcasHLA output.'''
    partial_files = []
    genotype_files = []
    
    for file in os.listdir(indir):
        if file.endswith('.partial_genotype.json'):
            partial_files.append(file)
        elif file.endswith('.genotype.json'):
            genotype_files.append(file)
            
    return genotype_files, partial_files 
    
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
    
    indir, outdir = [check_path(path) for path in [args.indir, args.outdir]]
    
    genotype_files, partial_files = get_paths(indir)
    
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
        
#-------------------------------------------------------------------------------
