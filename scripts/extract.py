#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   extract.py: extracts chromosome 6 reads from a BAM file for HLA genotyping.
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
import pickle
import json
import argparse
import logging as log

from datetime import date
from os.path import isfile
from argparse import RawTextHelpFormatter
from arcas_utilities import *

__version__     = '0.2.5'
__date__        = '2021-04-12'

#-------------------------------------------------------------------------------
#   Extract Reads
#-------------------------------------------------------------------------------

def index_bam(bam):
    '''Attempts to index BAM if .bai file is not found.'''
    if not isfile(''.join([bam, '.bai'])):
        run_command(['samtools', 'index', bam], 
                    '[extract] indexing bam: ')
                    
    if not isfile(''.join([bam, '.bai'])):
        sys.exit('[extract] Error: unable to index bam file.')
        

def extract_reads(bam, outdir, paired, unmapped, alts, temp, threads, sample):
    '''Extracts reads from chromosome 6 and alts/decoys if applicable.'''
    
    log.info(f'[extract] Extracting reads from {bam}')
    
    file_list = []
    #sample = os.path.splitext(os.path.basename(bam))[0]
    
    # Index bam
    index_bam(bam)
    
    hla_filtered = ''.join([temp, sample, '.hla.sam'])
    file_list.append(hla_filtered)
    hla_filtered_bam = ''.join([temp, sample, '.hla.bam'])
    file_list.append(hla_filtered_bam)
        
    # Get bam header to check for chromosome nomenclature
    output = run_command(['samtools', 'view', '-@'+threads, '-H', bam])
    header = output.stdout.decode('utf-8')
    
    if 'SN:chr' in header: 
        chrom = 'chr6'
    else: 
        chrom = '6'

    # Extract BAM header
    message = '[extract] Extracting chromosome 6: '
    command = ['samtools', 'view', '-H', '-@'+threads]
    command.extend([bam, '-o', hla_filtered])
    run_command(command, message)
    
    # Extracted reads mapped to chromosome 6
    message = '[extract] Extracting chromosome 6: '
    command = ['samtools', 'view', '-@'+threads]
    if paired: command.append('-f 2')
    else: command.append('-F 4')
    command.extend([bam, chrom, '>>', hla_filtered])
    run_command(command, message)
    
    # Extract unmapped reads
    if unmapped:
        message = '[extract] Extracting unmapped reads: '
        command = ['samtools', 'view', '-@'+threads]
        
        if paired: command.append('-f 12')
        else: command.append('-f 4')
        
        command.extend([bam, '>>', hla_filtered])
        run_command(command, message)
    
    # Check for alts in header and extract reads if present
    for alt in alts:
        if alt in header:
            command = ['samtools', 'view', '-@'+threads]
            
            if paired: command.append('-f 2')
            else: command.append('-F 4')
            
            command.extend([bam, alt+':', '>>', hla_filtered])
            run_command(command)


    # Convert SAM to BAM
    message = '[extract] Converting SAM to BAM: '
    command = ['samtools', 'view', '-Sb', '-@'+threads,
                hla_filtered, '>', hla_filtered_bam]    
    run_command(command, message)
            

    # Sort BAM
    hla_sorted = ''.join([temp, sample, '.hla.sorted.bam'])
    file_list.append(hla_sorted)
    message = '[extract] Sorting bam: '
    command = ['samtools', 'sort', '-n', '-@'+threads, 
                hla_filtered_bam, '-o', hla_sorted]
    run_command(command, message)

    # Convert BAM to FASTQ and compress
    message = '[extract] Converting bam to fastq: '
    command = ['bedtools', 'bamtofastq', '-i', hla_sorted]
    if paired:
        fq1 = ''.join([outdir, sample, '.extracted.1.fq'])
        fq2 = ''.join([outdir, sample, '.extracted.2.fq'])
        command.extend(['-fq', fq1, '-fq2', fq2])
        run_command(command, message)
        
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq1])
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq2])
        
    else:
        fq = ''.join([outdir, sample, '.extracted.fq'])
        command.extend(['-fq', fq])
        run_command(command, message)
        run_command(['pigz', '-f', '-p', threads, '-S', '.gz', fq])

#-------------------------------------------------------------------------------
#   Main
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA extract',
                                     usage='%(prog)s [options] BAM file',
                                     add_help=False,
                                     formatter_class=RawTextHelpFormatter)
    
    
    parser.add_argument('bam', 
                        type=str, 
                        help='/path/to/sample.bam')
    
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
                        
    parser.add_argument('--log', 
                    type=str,
                    help='log file for run summary\n  '+
                         'default: sample.extract.log\n\n',
                    default=None, 
                    metavar='')
                        
    parser.add_argument('--single', 
                        action = 'count',
                        help='single-end reads\n  default: False\n\n',
                        default=False)    
    
    parser.add_argument('--unmapped', 
                        action = 'count',
                        help='include unmapped reads\n  default: False\n\n',
                        default=False)    
                        
    parser.add_argument('-o', '--outdir', 
                        type=str,
                        help='out directory\n\n',
                        default='./', metavar='')
                        
    parser.add_argument('--temp', 
                        type=str,
                        help='temp directory\n\n',
                        default='/tmp/', metavar='')
                        
    parser.add_argument('--keep_files',
                        action = 'count',
                        help='keep intermediate files\n\n',
                        default=False)
   
    parser.add_argument('--sample',
                        type = str,
                        help = 'User defined sample id\n',
                        default='')

    parser.add_argument('-t',
                        '--threads', 
                        type = str,
                        default='1')
    
    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    outdir = check_path(args.outdir)
    temp = create_temp(args.temp)
    
    sample = os.path.basename(args.bam).split('.')[0]
    if (args.sample != ""):
      sample = args.sample

    datDir = os.path.dirname(os.path.realpath(__file__)) + '/../dat/'
    
    # Set up log file
    if args.log:
        log_file = args.log
    else:
        log_file = ''.join([outdir,sample,'.extract.log'])
    with open(log_file, 'w'):
        pass
    if args.verbose:
        handlers = [log.FileHandler(log_file), log.StreamHandler()]
        
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
    else:
        handlers = [log.FileHandler(log_file)]
            
        log.basicConfig(level=log.DEBUG, 
                        format='%(message)s', 
                        handlers=handlers)
        
    log.info('')
    hline()
    log.info(f'[log] Date: %s', str(date.today()))
    log.info(f'[log] Sample: %s', sample)
    log.info(f'[log] Input file: %s', args.bam)
    log.info('[log] Read type: {}-end'
             .format( 'paired' if not args.single else 'single'))
    hline()
    
    # Load names of regions outside chr6 with HLA loci
    #with open(datDir + '/info/decoys_alts.p', 'rb') as file:
    #    alts = pickle.load(file)
    with open(datDir + 'info/decoys_alts.json', 'r') as file:
        alts = json.load(file)

    extract_reads(args.bam,
                  outdir, 
                  not args.single,
                  args.unmapped,
                  alts,
                  temp,
                  args.threads,
									sample)
    
    remove_files(temp, args.keep_files)
    
    hline()
    log.info('')
#-------------------------------------------------------------------------------
