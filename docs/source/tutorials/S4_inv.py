#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json5


# The root directory of the project
example_path = '/Users/yf/3.Project/8.MCMTpy/MCMTpy-master/data/example_yunnan'


#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
# 1. get notebook path
notebook_path = sys.path[0]
os.chdir(notebook_path)
# change path to Green_Funcs
os.chdir('../YN.202105212148_Inv/dc_inv')

# 2.show build_GFs.json 
filename = 'sample_dc.json'
with open(filename, 'r',encoding='utf-8') as f1:
    sample_dc_json=json5.load(f1)

# 3.change parameters with absolute path
sample_dc_json["GFs_json_file"] = os.path.join(example_path,"Green_Funcs/build_GFs_new.json")
sample_dc_json["GFs_file"] =os.path.join(example_path,"Green_Funcs/GFs/GFs_splits" )          
sample_dc_json["Raw_data_file"] = os.path.join(example_path,"YN.202105212148_Inv/YN.202105212148_MCMTpy/YN.202105212148.h5")              
sample_dc_json["Output_path"] = os.path.join(example_path,"YN.202105212148_Inv/dc_inv/Output_YN.202105212148_dc")         
sample_dc_json["Raw_data_inv_info"] = os.path.join(example_path,"YN.202105212148_Inv/dc_inv/Raw_data_inv_info.txt")     
sample_dc_json["Raw_data_inv_info_writed"] = os.path.join(example_path,"YN.202105212148_Inv/dc_inv/Raw_data_inv_info_writed.txt") 
sample_dc_json_new = os.path.join(example_path,'YN.202105212148_Inv/dc_inv/sample_dc_new.json')

sample_dc_json["MPI_n"] = 4                         # CPU num
sample_dc_json["Chains_n"] = 4                      # Each MK-chains num
sample_dc_json["N"] = 1e1                           # each chain‘s search number

sample_dc_json["alpha_max"] = 100                   # The initial value of alpha
sample_dc_json["alpha_min"] = 0                     # The final value of alpha
sample_dc_json["a"] = 10                     
sample_dc_json["b"] = 10               

with open(sample_dc_json_new,'w') as f2:
    json5.dump(sample_dc_json, f2, indent=2)
f2.close()

# 4.run MCMTpy in shell
os.chdir(notebook_path)
os.chdir('../YN.202105212148_Inv/dc_inv')
# !MCMTpy  sample MH  -c ./sample_dc_new.json #> sample.log