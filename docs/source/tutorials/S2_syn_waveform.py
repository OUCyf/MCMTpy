#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json5
import shutil

# The root directory of the project
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'

#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
# change path to Green_Funcs
notebook_path = sys.path[0]
os.chdir(notebook_path)
os.chdir('../syn')
print(os.listdir())

# show build_GFs.json 
filename = 'syn.json'
with open(filename, 'r',encoding='utf-8') as f1:
    syn_json=json5.load(f1)

# change parameters with absolute path
syn_json["GFs_json_file"] = os.path.join(example_path,"Green_Funcs/build_GFs_new.json")
syn_json["Stf_file"] = os.path.join(example_path,"syn/Stf_file/Stf_file.sac")
syn_json["Output_path"] = os.path.join(example_path,"syn/Synthetic")       
syn_json_new = os.path.join(example_path,'syn/syn_new.json')

with open(syn_json_new,'w') as f2:
    json5.dump(syn_json, f2, indent=2)
f2.close()

shutil.rmtree('./Synthetic')
# !MCMTpy syn pyfk  -c ./syn_new.json  > syn.log