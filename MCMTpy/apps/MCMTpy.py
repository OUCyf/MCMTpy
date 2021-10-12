#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:00:28 2021

@author: Fu Yin (yinfu@mail.ustc.edu.cn) at USTC

This script:
    1) The coding style refers to the beat (https://hvasbath.github.io/beat/);
    2) main app function of MCMTpy;
    3) include user manual and help guide;
    4) include subcommands.

Modify history:
    1) Mar 17 16:00:28 2021    ||    Fu Yin at USTC    ||    The initial release.
    2) ...
    
"""


import sys
from MCMTpy.info import version
from MCMTpy import build_GFs
from MCMTpy import syn
from MCMTpy import sample
from MCMTpy import plot



#%%##########################################################################
#                         -----------------------
#                          1. MCMTpy UESR MANUAL
#                         -----------------------
#############################################################################
### 1.1. program_name
program_name = 'MCMTpy'

### 1.2. subcommand usages
subcommand_usages = {
    'build_GFs':      'build_GFs [pyfk|sem] -c <build_GFs.json>',
    'syn':            'syn [pyfk|sem] -c <syn.json>',
    'sample':         'sample [Grid|MH|HMC] -c <sample.json>',
    'plot':           'plot [pyfk|sem] -c <plot.json>',
}
subcommands = list(subcommand_usages.keys())

### 1.3. user manual
usage = program_name + ''' <subcommand> [options] ...-c  <file.json> 

Version:'''+version+'''
Author: Fu Yin
Email: yinfu@mail.ustc.edu.cn 


Subcommands:

    build_GFs       % create GFs datebase
    syn             % synthesize the waveform
    sample          % sample the parameters
    plot            % visualize the results

To further help:

    MCMTpy <subcommand> --help

''' 





#%%##########################################################################
#                         -------------------------
#                          2. Subcommands Function
#                         -------------------------
#############################################################################
### 2.1. build_GFs function
def command_build_GFs(args):
    
    if args[0]=='--help':
        sys.exit('build_GFs usage: '+subcommand_usages['build_GFs']+'\n'+
                 'for example: MCMTpy build_GFs [pyfk|sem] -c build_GFs.json ')
    
    else:
        if args[1]=='-c' and len(args)==3:
            method = args[0]                                                    # pyfk or sem
            filename = args[2]                                                  # file.json
            build_GFs(filename,method)
        else:
            sys.exit('MCMTpy: command error')



### 2.2. syn function
def command_syn(args):
    
    if args[0]=='--help':
        sys.exit('syn usage: '+subcommand_usages['syn']+'\n'+
                 'for example: MCMTpy syn [pyfk|sem] -c syn.json ')
    
    else:
        if args[1]=='-c' and len(args)==3:
            method = args[0]                                                    # pyfk or sem
            filename = args[2]                                                  # file.json
            syn(filename,method)
        else:
            sys.exit('MCMTpy: command error')



### 2.3. sample function
def command_sample(args):
    
    if args[0]=='--help':
        sys.exit('sample usage: '+subcommand_usages['sample']+'\n'+
                 'for example: MCMTpy sample [grid|MH|HMC] -c sample.json ')
    
    else:
        if args[1]=='-c' and len(args)==3:
            method = args[0]                                                    # pyfk or sem
            filename = args[2]                                                  # file.json
            sample(filename,method)
        else:
            sys.exit('MCMTpy: command error')



### 2.4. plot function
def command_plot(args):
    
    if args[0]=='--help':
        sys.exit('syn usage: '+subcommand_usages['plot']+'\n'+
                 'for example: MCMTpy plot [pyfk|sem] -c plot.json ')
    
    else:
        if args[1]=='-c' and len(args)==3:
            method = args[0]                                                    # pyfk or sem
            filename = args[2]                                                  # file.json
            plot(filename,method)
        else:
            sys.exit('MCMTpy: command error')





#%%##########################################################################
#                         ------------------
#                          3. Main Function
#                         ------------------
#############################################################################

def main():
    args = list(sys.argv)                                                       # funtion: <sys.argv> get parameters from command line. For example: mpirun -n 4 MCMTpy build_GFs pyfk  -c ./build_GFs.json, args will be = ['MCMTpy', 'build_GFs', 'pyfk', '-c', './build_GFs.json']

    if len(args) == 1:                                                          # when command line == MCMTpy
        sys.exit('Usage: %s' % usage)
    else:
        command = args[1]


    if (command in subcommands) and (args[-1] not in ('--help', '-h', 'help')): # execute subcommand
        globals()['command_' + command](args[2:])
    
    elif args[-1] in ('--help', '-h', 'help'):                                  # execute help command
        if command in subcommands:
            globals()['command_' + command](['--help'])
        sys.exit('Usage: %s' % usage)
    
    else:
        sys.exit('MCMTpy: error: no such subcommand: %s' % command)



if __name__ == '__main__':
    main()









