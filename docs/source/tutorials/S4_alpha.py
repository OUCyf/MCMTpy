#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np
import matplotlib.pyplot as plt

# parameters
N=5e3
a=1000                  # 0.2*N
b=100                   # 0.04*N
alpha_max=100
alpha_min=0

#------------------------------------------------------------#
# we expect no parameters need to be changed below
#------------------------------------------------------------#
tt=np.linspace(0, round(N), num=round(N))
yy=np.zeros(len(tt))
for i in range(0,len(tt)):
    yy[i]=(alpha_max-alpha_min)*(a+math.exp((tt[0]-a)/b))/(a+math.exp((tt[i]-a)/b))+alpha_min

fig, axs = plt.subplots(1,1)
axs.plot(tt,yy, linestyle='-', color='red', lw=3, alpha=0.6)
axs.set_xlabel("Sample Number",fontsize=17)
axs.set_ylabel("Alpha Value",fontsize=17)
fig.show()