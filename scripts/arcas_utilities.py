#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------
#   arcas_utilities.py: functions common to multiple arcasHLA scripts.
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
import re
import logging as log
import uuid
from subprocess import PIPE, STDOUT, run

__version__     = '0.1'
__date__        = '2019-03-14'

#-------------------------------------------------------------------------------

def process_allele(allele, n, keep_alpha = True):
    '''Lowers allele resolution to n-fields.'''
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

def get_gene(allele):
    '''Returns gene of an allele.'''
    return allele.split('*')[0]

def check_path(path):
    '''Check if path exists and if it is terminated by "/".'''
    if not os.path.isdir(path):
        run_command(['mkdir',path])
        
    if path[-1] != '/': path += '/'
        
    return path

def remove_files(files, keep_files):
    '''Removes intermediate files.'''
    if keep_files:
        return

    if type(files) == list:
        for file in files:
            run_command(['rm',file])
    else:
        run_command(['rm',files])
        
def run_command(command, message = ''):
    '''Outputs message and command to log, runs command and returns output.'''
    if type(command) == list:
        command = ' '.join(command)

    if message: 
        log.info(''.join([message,'\n\n\t', command,'\n']))
    
    output = run(command, shell=True, stdout=PIPE, stderr=PIPE)
    
    if message:
        stderr = '\t' + output.stderr.decode('utf-8')
        stderr = re.sub('\n','\n\t',stderr)
        if len(stderr) > 1:
            log.info(stderr)
        
    return output

def create_temp(temp):
    temp = check_path(temp)
    temp_folder = ''.join([temp,'arcas_' + str(uuid.uuid4())])
    return check_path(temp_folder)

def hline():
    log.info('-'*80)

#-------------------------------------------------------------------------------
