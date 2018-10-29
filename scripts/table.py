#-----------------------------------------------------------------------------
import os
import json
import pandas as pd
import argparse

from arcas_common import check_path

def get_paths(indir):
    partial_json_files = []
    json_files = []
    
    for file in os.listdir(indir):
        if file.endswith('genotype.partial.json'):
            partial_json_files.append(file)
        elif file.endswith('genotype.json'):
            json_files.append(file)
            
    return json_files, partial_json_files
    
def process_json(json_files, indir, outdir, run, partial = False):
    file_out = [outdir]
    if run:
        file_out.extend([run,'.'])
        
    if partial: file_out.append('genotypes.partial.json')
    else: file_out.append('genotypes.json')
    
    genotypes = dict()
    for file in json_files:
        sample = file.split('.')[0]
        file_path = indir + file

        with open(file_path,'r') as file:
            genotypes[sample] = json.load(file)
            
    with open(''.join(file_out), 'w') as file:
        json.dump(genotypes,file)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i',
                        '--indir',
                        type=str,
                        help='directory containing arcasHLA files')

    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='out directory',
                        default='./')
    
    parser.add_argument('--run',
                        type=str,
                        help='run name',
                        default='')
    
    parser.add_argument('-v',
                        '--verbose', 
                        action = 'count',
                        default=False)
    
    args = parser.parse_args()
    
    indir, outdir = [check_path(path) for path in [args.indir, args.outdir]]
    
    
    
    json_files, partial_json_files = get_paths(indir)

    if json_files:
        process_json(json_files, indir, outdir, args.run)
        
    if partial_json_files:
        process_json(partial_json_files, indir, outdir, args.run, partial = True)
        
#-----------------------------------------------------------------------------
