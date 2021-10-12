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

# 2. change path to Green_Funcs
os.chdir('../Green_Funcs')
print(os.listdir())

# 3. show build_GFs.json 
filename = 'build_GFs.json'
with open(filename, 'r',encoding='utf-8') as f1:
    gfs_json=json5.load(f1)
f1.close()

# 4. change parameters with absolute path
gfs_json['Source_Station_info'] = os.path.join(example_path,"Green_Funcs/Source_Station_info.txt")
gfs_json["DATADIR"] = os.path.join(example_path,"Green_Funcs/GFs1")
gfs_json["DATADIR_splits"] = os.path.join(example_path,"Green_Funcs1/GFs1/GFs_splits")
gfs_json["Velocity_model"] =  os.path.join(example_path,"v_model/v_model_yunnan.txt")        
gfs_json_new = os.path.join(example_path,'Green_Funcs/build_GFs_new.json')

with open(gfs_json_new,'w') as f2:
    json5.dump(gfs_json, f2, indent=2)
f2.close()

# !mpirun -n 4 MCMTpy build_GFs pyfk -c build_GFs_new.json  > gfs.log