import os
from subprocess import PIPE, run

def process_allele(allele,n,keep_alpha = True):
    temp = []
    alpha = False
    fields = allele.split(':')
    
    if fields[-1][-1].isalpha():
        if keep_alpha and len(fields) < 3:
            alpha = fields[-1][-1]
        fields[-1] = fields[-1][:-1]

    out = fields[:n]
    out = ':'.join(out)
    
    if alpha: out += alpha
    return out

def sort_set(s):
    return sorted(list(s))

def check_path(path):
    if not os.path.isdir(path):
        run_command(['mkdir',path],'',False)
        
    if path[-1] != '/': path += '/'
        
    return path

def remove_files(file_list,keep_files):
    if keep_files:
        return

    for file in file_list:
        run_command(['rm',file],'', False)
        
def run_command(command, message, verbose):
    if type(command) == list:
        command = ' '.join(command)
        
    if verbose:
        print(message,command)
        
    output = run(command, shell=True, stdout=PIPE, stderr=PIPE)
        
    return output