#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   convert.py: converts HLA nomenclature and resolution.
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
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA convert',
                                     usage='%(prog)s [options]',
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        help='tsv or json containing HLA genotypes\n\n', 
                        nargs='*',
                        type=str)
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        help='out directory\n\n',
                        default='',
                        metavar='')

    parser.add_argument('--resolution',
                    type = int,
                    help = 'output resolution (1,2,3,4), must be equal or lower to input resolution',
                    default = 3,
                    metavar = '')
    
    parser.add_argument('--grouping',
                        type = str,
                        help = 'allele grouping (g-group, p-group, supertype) [default = False]',
                        default = False,
                        metavar = '')
    
    parser.add_argument('--input_nomenclature',
                        type = str,
                        help = 'input nomenclature (current, pre-2009, polysolver), automatically detects current vs polysolver',
                        default = '',
                        metavar = '')
    
    parser.add_argument('--output_nomenclature',
                        type = str,
                        help = 'input nomenclature (current, polysolver)',
                        default = 'current',
                        metavar = '')
    
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
