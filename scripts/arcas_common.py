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