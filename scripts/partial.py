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
from itertools import combinations

from collections import Counter, defaultdict
from itertools import combinations

from build_reference import check_ref
from arcas_common import process_allele, check_path, remove_files, run_command
from genotype import assemble_transcripts, genotype_gene


#-----------------------------------------------------------------------------
# not mine
class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)

#-----------------------------------------------------------------------------
# Process transcript assembly output
#-----------------------------------------------------------------------------

def process_counts(count_file, eq_file, gene_list, allele_index, allele_lengths, exon_index, keep_files): 
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

    eqs = {str(i):defaultdict(list) for i in [['2'],['2','3'],['1','2'],['1','2','3'],['2','3','4'],['1','2','3','4'],['2','3','4','5'],['1','2','3','4','5'],['2','3','4','5','6'],['1','2','3','4','5','6'],['2','3','4','5','6','7'],['1','2','3','4','5','6','7']]}

    for eq,indices in eq_index.items():
        if [index for index in indices if not allele_index[index]]:
            continue

        genes = list({allele.split('*')[0] for index in indices for allele in allele_index[index]})
        exons = list({exon_index[index] for index in indices})
        if len(genes) == 1 and counts_index[eq] > 0:
            gene = genes[0]
            for exon in exons:
                exon_indices = list({index for index in indices if exon_index[index] == exon})
                eqs[exon][gene].append((exon_indices,counts_index[eq]))

    remove_files([count_file, eq_file], keep_files)
    return eqs

#-----------------------------------------------------------------------------
# Partial genotyping
#-----------------------------------------------------------------------------

def get_single_count(a,exon_group):
    #print(a)
    return sum([filtered_eqs[exon_group][gene][index][1] for index in allele_to_eq[exon_group][a]])

def get_pair_count(a1,a2,exon_group):
    indices = allele_to_eq[exon_group][a1] | allele_to_eq[exon_group][a2]
    return sum([filtered_eqs[exon_group][gene][index][1] for index in indices])

def get_nonshared_count(a1,a2,exon_group):
    a1_indices = allele_to_eq[exon_group][a1] - allele_to_eq[exon_group][a2]
    a2_indices = allele_to_eq[exon_group][a2] - allele_to_eq[exon_group][a1]
    return sum([filtered_eqs[exon_group][gene][index][1] for index in a1_indices]), sum([filtered_eqs[exon_group][gene][index][1] for index in a2_indices])

def type_partial(filtered_eqs,partial_exons,gene,prediction,population):
    #print('gene',gene)
    gene_exons = {'A':"['2', '3']",'B':"['2', '3']",'C':"['2', '3']",'DRB1':"['2']",'DQB1':"['2']"}

    results = genotype_gene(gene, filtered_eqs[gene_exons[gene]], allele_lengths, allele_index, 
                            population, prior,10e-7, 1000, 5, .1, None)
    #print(results)
    exon_groups = defaultdict(set)

    for _,alleles,_ in results[0]:
        #print((set(alleles) - set(prediction)) & partial_alleles)
        for allele in (set(alleles) - set(prediction)) & partial_alleles:
            for exon_group in filtered_eqs.keys():
                if exon_group[1:-1] in str(sorted(partial_exons[allele].keys())):
                    exon_groups[exon_group].add(allele)
                    
    overall_explained_counts = dict()
    #print(exon_groups)
    for exon_group in sorted(exon_groups.keys(),key = lambda x: len(x), reverse=False):
        if gene in ['A','B','C'] and exon_group == "['2']":
            continue

        possible_alleles = {allele for allele in exon_groups[exon_group] if allele_to_eq[exon_group][allele] != allele_to_eq[exon_group][prediction[0]] and allele_to_eq[exon_group][allele] != allele_to_eq[exon_group][prediction[1]]}
        if not possible_alleles:
            continue

        explained_counts = dict()

        a1,a2 = prediction
        min_count = min(get_single_count(prediction[0],exon_group),get_single_count(prediction[1],exon_group))

        possible_alleles &= {allele for allele in exon_groups[exon_group] if get_single_count(allele,exon_group) > 10 + min_count}

        total_counts = sum([count for _,count in filtered_eqs[exon_group][gene]])

        expained_count = get_pair_count(a1,a2,exon_group)
        explained_perc = round(expained_count/total_counts,8)
        explained_counts[(a1,a2)] = explained_perc


        for a1,a2 in combinations(set(prediction) | possible_alleles,2):
            expained_count = get_pair_count(a1,a2,exon_group)
            if expained_count/total_counts > explained_perc:
                explained_counts[(a1,a2)] = round(expained_count/total_counts,8)
        if not explained_counts:
            continue 
        top_perc = max(explained_counts.values())
        explained_counts = {key:value for key,value in explained_counts.items() if value == top_perc}
        
        
        if len(explained_counts) > 1:
            pair_freq = dict()
            for a1,a2 in explained_counts.keys():
                f1 = f2 = 0.
                if process_allele(a1,2) in prior: f1 = prior[process_allele(a1,2)][population]
                if process_allele(a2,2) in prior: f2 = prior[process_allele(a2,2)][population]
                pair_freq[(a1,a2)] = f1 + f2
                
            a1,a2 = sorted(pair_freq.items(), key=lambda x:x[1], reverse=True)[0][0]
            
            overall_explained_counts[explained_counts[(a1,a2)]] = (a1,a2)
            
            
        else:
            overall_explained_counts[sorted(explained_counts.items(), key= lambda x:x[1],reverse=True)[0][1]] = sorted(explained_counts.items(), key= lambda x:x[1],reverse=True)[0][0]

    if overall_explained_counts:

        returned = sorted(overall_explained_counts.items(),key=lambda x: x[0], reverse = True)[0][1]
        if returned == prediction: prediction
        return returned
    return prediction

#-----------------------------------------------------------------------------
# Runs genotyping
#-----------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='arcasHLA partial',add_help=False,formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('file', 
                        type=str, 
                        help='list of fastq files or ".partial.json" file', 
                        nargs='*')
                        
    parser.add_argument('-h',
                        '--help', 
                        action = 'help',
                        help='show this help message and exit\n\n',
                        default=argparse.SUPPRESS)
    
    parser.add_argument('-g',
                        '--genes', 
                        type=str,
                        help='comma separated list of HLA genes\n  default: A,B,C,DRB1,DQB1\n',
                        default='A,B,C,DRB1,DQB1', metavar='')
    
    parser.add_argument('-G',
                        '--genotype', 
                        type=str,
                        help='"genotype.json" file from arcasHLA genotype'
                        , metavar='')
    
    parser.add_argument('-p',
                        '--population', 
                        type=str,
                        help= 'sample population\n  default: Prior\n',
                        default='Prior', metavar='')
    
    parser.add_argument('-d',
                        '--database', 
                        type=str,
                        help='frequency database\n  default: gold_smoothed10\n',
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
                        default=5, metavar='')
    
    parser.add_argument('--drop_threshold', 
                        type=float,
                        help='%% of max abundance allele needs to not be dropped\n  default: 0.1\n\n',
                        default=0.1, metavar='')
    
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
                        default='1')

    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    args = parser.parse_args()
    
        
    sample = os.path.basename(args.file[0]).split('.')[0]
    
    temp, outdir = [check_path(path) for path in [args.temp,args.outdir]]
        
    prior = pd.read_csv('database/frequency/' + args.database + '.tsv', delimiter='\t')
    prior = prior.set_index('allele').to_dict('index')

    
    # loads reference information
    with open('database/reference/hla.partial.json','r') as file:
        commithash,(gene_list, allele_index, allele_lengths, partial_index, exon_index, partial_exons, partial_alleles) = json.load(file)
        partial_alleles = set(partial_alleles)
          
    # runs transcript assembly if intermediate json not provided
    reference = 'database/reference/hla.partial.idx'
    if not args.file[0].endswith('.json'):
        [count_file,eq_file],[n,l,s] = assemble_transcripts(args.file,
                                                            reference,
                                                            outdir, 
                                                            temp,
                                                            args.threads,
                                                            args.verbose, 
                                                            args.keep_files)
        
        eqs = process_counts(count_file,
                                                         eq_file,gene_list, 
                                                         allele_index, 
                                                         allele_lengths,
                                                         exon_index,
                                                         args.keep_files)

        with open(''.join([outdir,sample,'.partial.json']),'w') as file:
            json.dump([eqs,n,l,s],file,cls=SetEncoder)
            
    else:
        if args.verbose: print('[genotype] loading intermediate json')
        with open(args.file[0],'r') as file:
            eqs,n,l,s = json.load(file)
            
    with open(args.genotype,'r') as file:
        predictions = json.load(file)
    
    if args.genes == 'none':
        sys.exit(0)
    
    # sets up gene list
    if args.genes == 'all': 
        genes = sorted(gene_list)
    else:                                                                              
        genes = sorted(set(args.genes.split(',')) & set(gene_list))
        
    all_predicted = {allele[0] if type(allele) == list else allele for alleles in predictions.values() for allele in alleles}
    wanted_indices = {index for index, alleles in allele_index.items() if alleles and (set(alleles) & (partial_alleles | set(all_predicted)))}
    filtered_eqs = dict()
    for exon_group,eq_list in eqs.items():
        filtered_eqs[exon_group] = dict()
        for gene in predictions:
            if gene not in eq_list:
                continue

            filtered = []
            for indices, count in eq_list[gene]:
                indices = set(indices) & wanted_indices
                if not indices:
                    continue

                filtered.append([indices,count])
            filtered_eqs[exon_group][gene] = filtered


    allele_to_eq = {exon_group:defaultdict(set) for exon_group in filtered_eqs.keys()}
    for exon_group,eq_list in filtered_eqs.items():
        for gene in predictions:
            if gene not in eq_list:
                continue

            for i,(indices, count) in enumerate(eq_list[gene]):
                for index in indices:
                    for allele in allele_index[index]:
                        allele = process_allele(allele,3)
                        allele_to_eq[exon_group][allele].add(i)



    individual_results = dict()

    for gene in genes:
        print(gene)
        if type(predictions[gene][0][0]) == str:
            prediction = [predictions[gene][0][0],predictions[gene][0][0]]
        else:
            prediction = [alleles[0] for alleles in predictions[gene][0]]
            
        if len(prediction) == 1:
            prediction.append(prediction[0])
            
        individual_results[gene] = type_partial(filtered_eqs,partial_exons,gene,prediction,args.population)
    
            
    with open(''.join([outdir,sample,'.genotype.partial.json']),'w') as file:
            json.dump(individual_results,file,cls=SetEncoder)

#-----------------------------------------------------------------------------
