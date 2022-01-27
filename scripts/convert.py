#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   convert.py: changes HLA resolution and converts nomenclature to p/g-groups.
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

import sys
import os
import json
import argparse
import pandas as pd
import pickle

from argparse import RawTextHelpFormatter
from os.path import isfile, isdir, dirname, realpath

from arcas_utilities import check_path, process_allele

#-------------------------------------------------------------------------------

__version__     = '0.4.0'
__date__        = '2022-01-27'

#-------------------------------------------------------------------------------
#   Paths and fileames
#-------------------------------------------------------------------------------

rootDir = dirname(realpath(__file__)) + '/../'

hla_convert_json = rootDir + 'dat/ref/hla.convert.json'

#-------------------------------------------------------------------------------
        
def convert_allele(allele, resolution):
    '''Checks nomenclature of input allele and returns converted allele.'''
    i = len(allele.split(':'))
        
    # Input: P-group allele
    if allele[-1] == 'P':
        if resolution == 'g-group': 
            sys.exit('[convert] Error: p-group cannot be converted ' +
                     'to g-group.')

        # Output: 1-field allele unless forced
        elif type(resolution) == int:
            if resolution > 1 and not args.force:
                sys.exit('[convert] Error: p-group cannot be ' +
                         'converted to %.0f fields.' %resolution)
            allele = process_allele(allele[:-1], resolution)

    # Input: G-group allele
    elif allele[-1] == 'G':

        # Output: 1-field allele unless forced
        if type(resolution) == int:
            if resolution > 1 and not args.force:
                sys.exit('[convert] Error: g-group cannot be converted' +
                         'to %.0f fields.' %resolution)
            allele = process_allele(allele[:-1], resolution)
            
        # Output: P-group allele
        elif resolution == 'p-group': 
            if allele[:-1] in p_group[i]:
                allele = p_group[i][allele[:-1]]

            elif process_allele(allele[:-1], i - 1) in p_group[i]:
                allele = p_group[i][process_allele(allele[:-1], i -1)]

    # Input: ungrouped allele
    # Output: G-group allele
    elif resolution == 'g-group':
        if allele in g_group[i]:
            allele = g_group[i][allele]
        elif allele[-1] != 'N':
            allele = process_allele(allele,3)

    # Input: ungrouped allele
    # Output: P-group allele
    elif resolution == 'p-group':
        if allele in p_group[i]:
            allele = p_group[i][allele]
            
    # Input: ungrouped allele
    # Output: reduced resolution, ungrouped allele
    elif type(resolution) == int:
        allele = process_allele(allele, resolution)
        
    return allele

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA convert',
                                     usage='%(prog)s [options]',
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        help='tsv containing HLA genotypes, see github for ' +
                              'example file structure.\n\n',
                        type=str)
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('-r',
                        '--resolution',
                        help = 'output resolution (1,2,3) or grouping ' + 
                               '(g-group, p-group)\n\n',
                        metavar='')

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        help='output file\n  default: ' +
                             './file_basename.resolution.tsv\n\n',
                        default='',
                        metavar='')
    
    parser.add_argument('-f',
                        '--force',
                        help = 'force conversion for grouped alleles even if ' +
                                'it results in loss of resolution',
                        action = 'count',
                        default=False)
    
    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    #p_group, g_group = pickle.load(open(hla_convert,'rb'))
    #to do, test this
    with open(hla_convert_json, 'r') as file:
        p_group,g_group = json.load(file)
    
    # Check input resolution
    accepted_fields = {'1','2','3','4'}
    accepted_groupings = {'g-group','p-group'}
        
    resolution = None

    if args.resolution in accepted_fields:
        resolution = int(args.resolution)
    elif args.resolution.lower() in accepted_groupings:
        resolution = args.resolution.lower()
        
    if not resolution:
        sys.exit('[convert] Error: output resolution is needed ' +
                 '(1, 2, 3, g-group, p-group).')
    
    
    # Create outfile name
    if not args.outfile:
        outfile = [os.path.splitext(os.path.basename(args.file))[0],
                   args.resolution.lower(),
                   'tsv']
        outfile = '.'.join(outfile)
    else:
        outfile = args.outfile
        
    # Load input genotypes
    df_genotypes = pd.read_csv(args.file, sep = '\t').set_index('subject')
    genotypes = df_genotypes.to_dict('index')
    

    for subject, genotype in genotypes.items():
        for gene, allele in genotype.items():
            if type(allele) != str:
                continue

            genotypes[subject][gene] = convert_allele(allele, resolution)

    pd.DataFrame(genotypes).T.rename_axis('subject').to_csv(outfile, sep = '\t')
        
#-------------------------------------------------------------------------------
