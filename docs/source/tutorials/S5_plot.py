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

# 2.change path to Green_Funcs
os.chdir(notebook_path)
os.chdir('../YN.202105212148_Inv/dc_inv')
print(os.listdir())

# 3.read plot_dc.json 
filename = 'plot_dc.json'
with open(filename, 'r',encoding='utf-8') as f1:
    plot_dc_json=json5.load(f1)

# 4.change parameters with absolute path
plot_dc_json["plot_output_path"] = os.path.join(example_path,"YN.202105212148_Inv/dc_inv/figure_dc")
plot_dc_json["Inv_json_file"] = os.path.join(example_path,"YN.202105212148_Inv/dc_inv/sample_dc_new.json") 
plot_dc_json_new = os.path.join(example_path,'YN.202105212148_Inv/dc_inv/plot_dc_new.json')

# 5.hist
plot_dc_json["plot_hist"] =True
plot_dc_json["N_start"] = 0 
plot_dc_json["N_start_accept"] = 1
plot_dc_json["num_bins"] =  50     
plot_dc_json["num_std"] =5      
plot_dc_json["labels_name"] = ['Mw','strike/°','dip/°','rake/°','z/km']

# 6. misfit
plot_dc_json["plot_misfit"] = True
plot_dc_json["MPI_n_st"] = 0
plot_dc_json["Chains_n_st"] = 0

# 7. waveform
plot_dc_json["plot_waveform"] = True
plot_dc_json["FM_best"] = [6.68,  135.1,  87.0,  -168.1,    25.58,  99.86,  5.95,  -0.57 ]     
plot_dc_json["line_n_sta"] = 3                                     # Draw the data of several stations in a row
plot_dc_json["max_p_ylim"] = 1                                     # max amp y-axis of p wave
plot_dc_json["max_s_ylim"] = 1.0
plot_dc_json["max_surf_ylim"] = 1
plot_dc_json["plot_comp"] = [[1,1,0],[0,0,0],[1,1,1]]              # What components do you want to draw? P、S、Surf's Z/R/T


with open(plot_dc_json_new,'w') as f2:
    json5.dump(plot_dc_json, f2, indent=2)
f2.close()

# 8.run MCMTpy in shell
os.chdir(notebook_path)
os.chdir('../YN.202105212148_Inv/dc_inv')
# !MCMTpy plot pyfk -c plot_dc_new.json > plot.log