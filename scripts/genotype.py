#-----------------------------------------------------------------------------
import os
import sys
import re
import json
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from subprocess import PIPE, run
from textwrap import wrap
import math
import pandas as pd

from collections import Counter, defaultdict
from itertools import combinations

from build_reference import check_ref
from arcas_common import process_allele, check_path, remove_files, run_command


#-----------------------------------------------------------------------------
# not mine
class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)
#
#-----------------------------------------------------------------------------
# Process BAM input
#-----------------------------------------------------------------------------

# Returns number of reads, average read length, read length standard deviation
def analyze_reads(fqs,temp,reads_file,keep_files,verbose):
    
    awk = "| awk '{if(NR%4==2) print length($1)}'"
    
    # this breaks if files aren't gzipped
    if verbose: print('[align] analyzing read length')
    if len(fqs) == 2:
        fq1, fq2 = fqs

        command = ['zcat <', fq1, awk, '>' , reads_file]
        run_command(command, '', False)
        
        command = ['zcat <',fq2,awk, '>>', reads_file]
        run_command(command, '', False)
        
    else:
        fq = fqs[0]
        
        command = ['zcat <',fq,awk,'>',reads_file]
        run_command(command, '', False)
        
    read_lengths = np.genfromtxt(reads_file)
    
    n = len(read_lengths)
    l = np.mean(read_lengths)
    s = np.std(read_lengths)
    
    if verbose: print(f'[align] {n} reads, average length of {l:.2f} w/ std {s:.2f}')
    
    remove_files([reads_file],keep_files)
    
    return n,l,s

# Runs transcript pseudo-alignment/transcript assembly
def assemble_transcripts(fqs, reference, outdir,temp,threads,verbose,keep_files, partial = False):
    file_list = []
    sample = os.path.splitext(os.path.basename(fqs[0]))[0].split('.')[0]
    
    reads_file = ''.join([temp,sample,'.reads.txt'])
    n,l,s = analyze_reads(fqs,temp,reads_file,keep_files,verbose)
    
    # kallisto breaks if single-end read length standard deviation is 0
    if s == 0: s = .00001

    temp2 = ''.join([temp,sample,'/'])
    command = ['kallisto pseudo -i',reference,'-t',threads,'-o',temp2]
    
    if not os.path.isdir(temp2): 
        run(' '.join(['mkdir',temp2]), shell=True, stdout=PIPE, stderr=PIPE)
        
    if len(fqs) == 2:
        command.extend([fqs[0],fqs[1]])
    else:
        fq = fqs[0]
        command.extend(['--single -l',str(l),'-s',str(s),fq])
        
    run_command(command, '[align] aligning reads with Kallisto:', verbose)

    
    file_in = ''.join([temp2,'pseudoalignments.tsv'])
    count_file = ''.join([temp,sample,'.counts.tsv'])
    file_list.append(file_in)
    
    command = ['mv',file_in,count_file]
    run_command(command, '', False)

    
    file_in = ''.join([temp2,'pseudoalignments.ec'])
    eq_file = ''.join([temp,sample,'.eq.tsv'])
    file_list.append(file_in)
    
    command = ['mv',file_in,eq_file]
    run_command(command, '', False)

    
    run_command(['rm -rf',temp2], '', False)
    
    remove_files(file_list,keep_files)
           
    return [count_file,eq_file],[n,l,s]

#-----------------------------------------------------------------------------
# Process transcript assembly output
#-----------------------------------------------------------------------------

def process_counts(count_file, eq_file, gene_list, allele_index, allele_lengths, verbose, keep_files): 
    if verbose: print('[align] processing counts')
    counts_index = dict()
    with open(count_file,'r') as file:
        for line in file.read().splitlines():
            eq, count = line.split('\t')
            counts_index[eq] = float(count)

    eq_index = dict()
    with open(eq_file,'r') as file:
        for line in file.read().splitlines():
            eq, indices = line.split('\t')
            eq_index[eq] = indices.split(',')

    eqs = defaultdict(list)

    for eq,indices in eq_index.items():
        if [index for index in indices if not allele_index[index]]:
            continue

        genes = list({allele.split('*')[0] for index in indices for allele in allele_index[index]})
        if len(genes) == 1 and counts_index[eq] > 0:
            gene = genes[0]
            eqs[gene].append((indices,counts_index[eq]))

    remove_files([count_file, eq_file], keep_files)
    return eqs

#-----------------------------------------------------------------------------
# Genotype
#-----------------------------------------------------------------------------

def expectation_maximization(eqs, lengths, allele_index, 
                             population, prior,
                             tolerance, max_iterations, 
                             drop_iterations, drop_threshold, verbose):
    

    def initial_abundances(eqs,lengths,population):
        def divide_equally(alleles, count):
            n_alleles = len(alleles)
            for allele in alleles:
                counts[allele] += count / n_alleles
                
        def divide_prior(alleles, count, allele_prob):
            total_prob = sum(allele_prob.values())
            for allele in alleles:
                counts[allele] += count * (allele_prob[allele]/total_prob)
        
        counts = defaultdict(float)
        for alleles,count in eqs:
            if population:
                allele_prior = defaultdict(float)
                for index in alleles:
                    allele = process_allele(allele_index[index][0],2)
                    if allele in prior:
                        allele_prior[index] = prior[allele][population]

                if allele_prior:
                    divide_prior(alleles, count, allele_prior)
                    continue
                    
            divide_equally(alleles, count)
                
        return counts_to_abundances(counts)
    
    def counts_to_abundances(counts):
        abundances = defaultdict(float)
        for allele,count in counts.items():
            length = lengths[allele]
            abundances[allele] = count / length
        
        total_abundance = sum(abundances.values())
        
        for allele,abundance in abundances.items():
            abundances[allele] = abundance / total_abundance
        
        return abundances
    
    def update_abundances(eqs, abundances):
        counts = defaultdict(float)
        for alleles, count in eqs:
            alleles = [allele for allele in alleles if allele in abundances]
            total_abundance = sum([abundances[allele] for allele in alleles])
            
            if total_abundance == 0:
                continue
            
            for allele in alleles:
                counts[allele] += count * (abundances[allele]/total_abundance)
        return counts_to_abundances(counts)
    
    def drop_alleles(eqs, abundances, drop_iterations, drop_threshold, iterations, converged):
        #if iterations <= 1:
        #    return abundances, eqs
        if iterations < drop_iterations and not converged:
            abundances = {allele:abundance for allele,abundance in abundances.items() if abundance > 0.0}
            return abundances, eqs


        threshold = drop_threshold * max(abundances.values())
        abundances = {allele:abundance for allele,abundance in abundances.items() if abundance >= threshold}
        #if iterations != drop_iterations:
        return abundances, eqs
        
        eqs_ = []
        for alleles, count in eqs:
            alleles = [allele for allele in alleles if allele in abundances]
            if not alleles:
                continue
            threshold = drop_threshold * max([abundances[allele] for allele in alleles])
            alleles = [allele for allele in alleles if abundances[allele] >= threshold]
            eqs_.append((alleles,count))
        return abundances, eqs_
    
    def SRSS(theta):
        square_sum = 0.0
        for i in theta:
            square_sum += i**2
        return math.sqrt(square_sum)

    def check_convergence(theta0, theta_prime):
        diff = [theta_prime[allele] - theta0[allele] for allele in theta0.keys()]
        residual_error = SRSS(diff)
        return residual_error < tolerance

    converged = False
    iterations = 0

    theta0 = initial_abundances(eqs, lengths, population)

    
    while iterations < max_iterations and not converged:
        #theta_prime = update_abundances(eqs, theta0)
        #converged = check_convergence(theta0, theta_prime)
        #theta0, eqs = drop_alleles(eqs, theta_prime, drop_iterations, drop_threshold, iterations, converged)
        #iterations += 1
        #continue
        if verbose: print(f'[genotype] EM iterations: {iterations}', end="\r")

        theta1 = update_abundances(eqs, theta0)
        theta2 = update_abundances(eqs, theta1)
        theta_prime = defaultdict(float)

        r = dict()
        v = dict()
        sum_r = 0.0
        sum_v = 0.0
        for allele in theta1:
            r[allele] = theta1[allele] - theta0[allele]
            v[allele] = (theta2[allele] - theta1[allele]) - r[allele]
            
        srss_r = SRSS(r.values())
        srss_v = SRSS(v.values())

        if srss_v != 0:
            x = 0
            alpha = -(srss_r / srss_v)
            for allele in r:
                value = theta0[allele] - 2 * alpha * r[allele] + (alpha**2) * v[allele]
                theta_prime[allele] = value
                
            step_min = min(theta_prime.values())
            step_max = max(theta_prime.values())

            if step_min < 0:
                theta_prime = {allele:(value-step_min)/(step_max-step_min) for allele,value in theta_prime.items()}
                total = sum(theta_prime.values())
                theta_prime = {allele:value/total for allele,value in theta_prime.items()}
     
            theta_prime = update_abundances(eqs, theta_prime)
        else:
            theta_prime = theta1
            
        converged = check_convergence(theta0, theta_prime)

        theta0, eqs = drop_alleles(eqs, theta_prime, drop_iterations, drop_threshold, iterations, converged)
        #print([[allele_index[index],count] for index,count in sorted(theta0.items(),key=lambda x:x[1],reverse=True)])
        iterations += 1
                
    return theta0


def genotype_gene(gene, eqs, lengths, allele_index,
                  population, prior, 
                  tolerance, max_iterations, 
                  drop_iterations, drop_threshold, verbose):

    if gene not in {'A','B','C','DRB1','DQB1','DQA1'}: population = None
        
    if verbose: print('\n[genotype] genotyping HLA-' + gene)
    em_results = expectation_maximization(eqs[gene], lengths, allele_index, population, prior, tolerance, 
                                          max_iterations, drop_iterations, drop_threshold, verbose)


    em_results = [[index,allele_index[index],abundance] for index,abundance in em_results.items()]
    
    if verbose:
        print('\n[genotype] top alleles by abundance:')
        print('           {: <20}\t{: >9}\t'.format('allele','abundance'))
        for _, alleles, abundance in sorted(em_results,key=lambda x: x[2], reverse = True):
            print('           {: <20}\t{:>8.2f}%\t'.format(process_allele(alleles[0],3),abundance*100))
    
    if len(em_results) > 1:
        grouped_indices = defaultdict(set)
        for index, alleles, abundances in em_results:
            allele = process_allele(alleles[0],2)
            grouped_indices[allele].add(index)

        grouped_indices = [tuple(indices) for indices in grouped_indices.values()]

        if len(grouped_indices) > 1:

            grouped = dict()

            for a1,a2 in combinations(grouped_indices,2):

                total = sum([count for indices,count in eqs[gene]])
                a1_or_a2 = sum([count for indices,count in eqs[gene] if set(a1) & set(indices) or set(a2) & set(indices)])

                grouped[(a1,a2)] = (a1_or_a2/total)
        else:
            grouped = dict()

            a1,a2 = sorted(list(grouped_indices)[0])[:2]

            total = sum([count for indices,count in eqs[gene]])
            a1_or_a2 = sum([count for indices,count in eqs[gene] if set(a1) & set(indices) or set(a2) & set(indices)])
            
            grouped[((a1,),(a2,))] = (a1_or_a2/total)
            
        a1,a2 = sorted(grouped.items(),key = lambda x: x[1],reverse = True)[0][0]
        
        total = sum([count for indices,count in eqs[gene]])

        a1_unique = sum([count for indices,count in eqs[gene] if set(a1) & set(indices) and not set(a2) & set(indices)])
        a2_unique = sum([count for indices,count in eqs[gene] if set(a2) & set(indices) and not set(a1) & set(indices)])
        
        a1 = [sorted(a1)[0]]
        a2 = [sorted(a2)[0]]
        
        #check for 0
        if a1_unique == a2_unique == 0:
            genotype = a1 + a2
        elif a1_unique == 0:
            genotype = a2
        elif a2_unique == 0:
            genotype = a1
        elif min(a1_unique/a2_unique,a2_unique/a1_unique) < 1.5*drop_threshold:
            if a1_unique > a2_unique:
                genotype = a1
            else:
                genotype = a2
        else:
            genotype = a1 + a2
        
        genotype = [process_allele(allele_index[index][0],3) for index in genotype]
        
        
        
    else:
        index, alleles, _ = em_results[0]
        genotype = [process_allele(alleles[0],3),]
        
        
    if verbose:
        print('[genotype] most likely genotype:')
        for allele in genotype:
            print(f'           {allele}')
    
    return em_results, genotype

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

if __name__ == '__main__':
    
    genes = ['all','A', 'B', 'C', 'DMA', 'DMB', 'DOA', 'DOB', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRA', 'DRB1', 'DRB3', 'DRB5', 'E', 'F', 'G', 'H', 'J', 'K', 'L']
    populations = ['asian_pacific_islander','black','caucasian','hispanic','native_american','prior']
    databases = ['unsmoothed','gold_smoothed10']
    
    parser = argparse.ArgumentParser(prog='arcasHLA genotype',add_help=False,formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        type=str, 
                        help='list of fastq files or json file', 
                        nargs='*')
                        
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('-g',
                        '--genes', 
                        type=str,
                        help='comma separated list of HLA genes\n  default: all\n' + '\n  '.join(wrap('  options: ' + ', '.join(genes),60)) +'\n\n',
                        default='all', metavar='')
    
    parser.add_argument('-p',
                        '--population', 
                        type=str,
                        choices = populations,
                        help= 'sample population\n  default: None\n' + '\n  '.join(wrap('  options: ' + ', '.join(populations),60)) +'\n\n',
                        default='prior', metavar='')
    
    parser.add_argument('-d',
                        '--database', 
                        type=str,
                        choices = databases,
                        help='frequency database\n  default: gold_smoothed10\n  options: ' + ', '.join(databases) + '\n\n',
                        default='gold_smoothed10', metavar='')
   
    
    parser.add_argument('--tolerance', 
                        type=float,
                        help='convergence tolerance\n  default: 10e-7\n\n',
                        default=10e-7, metavar='')
    
    parser.add_argument('--max_iterations', 
                        type=int,
                        help='maximum # of iterations\n  default: 1000\n\n',
                        default=1000, metavar='')
    
    parser.add_argument('--drop_iterations', 
                        type=int,
                        help='# of iterations before dropping low-support alleles\n  default: 20\n  recommended paired: 20\n  recommended single: 2\n\n',
                        default=20, metavar='')
    
    parser.add_argument('--drop_threshold', 
                        type=float,
                        help='%% of max abundance allele needs to not be dropped\n  default: 0.1\n\n',
                        default=0.1, metavar='')
    
    parser.add_argument('--zygosity_threshold', 
                        type=float,
                        help='%% of unique reads needed to call genotype heterozygous\n  default: 0.15\n\n',
                        default=0.15, metavar='')
    
    
    '''
    parser.add_argument('--bootstrap', 
                        type = int,
                        help = 'number of bootstrap iterations to perform\n  default: 0\n\n',
                        default=False, metavar='')
    '''
    
    parser.add_argument('-o',
                        '--outdir',
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
                        default='1', metavar='')

    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    
    args = parser.parse_args()
    
    if len(args.file) == 0:
        print('genotype.py: error: the following arguments are required: file')
        sys.exit()
    
    sample = os.path.basename(args.file[0]).split('.')[0]
    
    
    temp, outdir = [check_path(path) for path in [args.temp,args.outdir]]
        
    prior = pd.read_csv('database/frequency/' + args.database + '.tsv', delimiter='\t')
    prior = prior.set_index('allele').to_dict('index')
        
    # checks if HLA reference exists
    check_ref(args.verbose)
    
    # loads reference information
    with open('database/reference/hla.json','r') as file:
        commithash,(gene_list, allele_index, allele_lengths) = json.load(file)
          
    # runs transcript assembly if intermediate json not provided
    reference = 'database/reference/hla.idx'
    if not args.file[0].endswith('.json'):
        [count_file,eq_file],[n,l,s] = assemble_transcripts(args.file,
                                                            reference,
                                                            outdir, 
                                                            temp,
                                                            args.threads,
                                                            args.verbose, 
                                                            args.keep_files)
        
        eq_index = process_counts(count_file,
                                                         eq_file,gene_list, 
                                                         allele_index, 
                                                         allele_lengths,
                                                         args.verbose,
                                                         args.keep_files)

        with open(''.join([outdir,sample,'.json']),'w') as file:
            json.dump([eq_index,n,l,s],file,cls=SetEncoder)
            
    else:
        if args.verbose: print('[genotype] loading intermediate json')
        with open(args.file[0],'r') as file:
            eq_index,n,l,s = json.load(file)
    
    if args.genes == 'none':
        sys.exit(0)
    
    # sets up gene list
    if args.genes == 'all': 
        genes = sorted(gene_list)
    else:                                                                              
        genes = sorted(set(args.genes.split(',')) & set(gene_list))
        
    em_results = dict()
    genotypes = dict()
    
    for gene in genes:
        total_count = sum([count for _,count in eq_index[gene]])
        if not total_count:
            continue
            
        em_results[gene],genotypes[gene] = genotype_gene(gene, eq_index, allele_lengths, allele_index,
                  args.population, prior, 
                  args.tolerance, args.max_iterations, 
                  args.drop_iterations, args.drop_threshold,
                  args.verbose)
        
    with open(''.join([outdir,sample,'.em.json']),'w') as file:
            json.dump(em_results,file,cls=SetEncoder)
            
    with open(''.join([outdir,sample,'.genotype.json']),'w') as file:
            json.dump(genotypes,file,cls=SetEncoder)

#-----------------------------------------------------------------------------
