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
from subprocess import PIPE, run

#-------------------------------------------------------------------------------

def process_allele(allele,n,keep_alpha = True):
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

def remove_files(file_list,keep_files):
    '''Removes intermediate files.'''
    if keep_files:
        return

    for file in file_list:
        run_command(['rm',file])
        
def run_command(command, message = ''):
    '''Outputs message and command to log, runs command and returns output.'''
    if type(command) == list:
        command = ' '.join(command)

    if message: log.info(message + command)
        
    output = run(command, shell=True, stdout=PIPE, stderr=PIPE)
        
    return output

#-------------------------------------------------------------------------------
