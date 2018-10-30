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
import re
import json
import argparse
from argparse import RawTextHelpFormatter
from arcas_common import process_allele, check_path, remove_files, run_command


def extract_reads(bam, outdir, paired, temp, threads, verbose, keep_files):
    if verbose: print('[extract] extracting reads from', bam)
    
    file_list = []
    sample = os.path.splitext(os.path.basename(bam))[0]

    if not os.path.isfile(''.join([bam,'.bai'])):
        command = ' '.join(['samtools','index',bam])
        run_command(command, '[extract] indexing bam: ', verbose)
    if not os.path.isfile(''.join([bam,'.bai'])):
        print('extract: error - unable to index file')
        return

        
    output = run_command(['samtools','view','-@'+threads,'-H',bam],'',False).stdout.decode('utf-8')
    
    if 'SN:chr' in output: chrom = 'chr6'
    else: chrom = '6'

    hla_filtered = ''.join([temp,sample,'.hla.sam'])
    file_list.append(hla_filtered)
    
    hla_filtered_bam = ''.join([temp,sample,'.hla.bam'])
    file_list.append(hla_filtered_bam)
    
    command = ['samtools','view','-H','-@'+threads]
    if paired: command.append('-f 2')
    command.extend([bam,'-o',hla_filtered])
    run_command(command, '[extract] extracting chromosome 6: ', verbose)
    
    
    command = ['samtools','view','-@'+threads]
    if paired: command.append('-f 2')
    command.extend([bam,chrom,'>>',hla_filtered])
    run_command(command, '[extract] extracting chromosome 6: ', verbose)
    
    for chrom in ['HSCHR6_MHC_APD','HSCHR6_MHC_COX','HSCHR6_MHC_DBB','HSCHR6_MHC_MANN','HSCHR6_MHC_MCF','HSCHR6_MHC_QBL','HSCHR6_MHC_SSTO']:
        if chrom in output:
            command = ['samtools','view','-@'+threads]
            if paired: command.append('-f 2')
            command.extend([bam,chrom,'>>',hla_filtered])
            run_command(command, ' '.join(['[extract] extracting decoy ',chrom,': ']), verbose)


    command = ['samtools','view','-Sb','-@'+threads,hla_filtered,'>',hla_filtered_bam]    
    run_command(command, '[extract] converting SAM to BAM: ', verbose)
            
    hla_sorted = ''.join([temp,sample,'.hla.sorted.bam'])
    file_list.append(hla_sorted)
    command = ['samtools','sort','-n','-@'+threads,hla_filtered_bam,'-o',hla_sorted]
    #command = ['bamsort','SO=queryname','threads='+threads,'<',hla_filtered_bam,'>',hla_sorted]
    run_command(command, '[extract] sorting bam: ', verbose)

    if paired:
        fq1 = ''.join([outdir,sample,'.1.fq'])
        fq2 = ''.join([outdir,sample,'.2.fq'])
        command = ['bedtools','bamtofastq','-i',hla_sorted,'-fq',fq1,'-fq2',fq2]
        #command = ['bamtofastq','threads='+threads,'filename='+hla_sorted,'F='+fq1,'F2='+fq2]
        run_command(command, '[extract] converting bam to fastq: ', verbose)
        
        run_command(' '.join(['pigz','-f','-p',threads,'-S','.gz',fq1]), '', False)
        run_command(' '.join(['pigz','-f','-p',threads,'-S','.gz',fq2]), '', False)
        
    else:
        fq = ''.join([outdir,sample,'.fq'])
        command = ['bedtools','bamtofastq','-i',hla_sorted,'-fq',fq]
        run_command(command, '[extract] converting bam to fastq: ', verbose)
        run_command(' '.join(['pigz','-f','-p',threads,'-S','.gz',fq]),'', False)

    remove_files(file_list,keep_files)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog='arcasHLA extract',
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
                        
    parser.add_argument('--paired', 
                        action = 'count',
                        help='paired-end reads\n  default: False\n\n',
                        default=False)               
                        
    parser.add_argument('-o','--outdir', 
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
    
    parser.add_argument('-t',
                        '--threads', 
                        type = str,
                        default='1')
    
    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    if args.temp[-1] != '/': args.temp += '/'
    if args.outdir[-1] != '/': args.outdir += '/'
    
    extract_reads(args.bam, 
                  args.outdir, 
                  args.paired, 
                  args.temp,
                  args.threads,
                  args.verbose, 
                  args.keep_files)
