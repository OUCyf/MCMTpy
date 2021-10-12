#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:34:20 2021

@author: yf
"""

import numpy as np
import math,os,sys
import random
import shutil
from mpi4py import MPI
from tqdm import tqdm

from MCMTpy.gfs.pyfk_GFs.read_json import read_GFs_json
from MCMTpy.sampler.pyfk_MH.read_sampler_json import read_Inv_json
from MCMTpy.sampler.pyfk_MH.sampler_module import config_fk,get_GFs_sta_info,get_Sta_rand_data,write_inv_para,get_sigma,get_fixed_fm
from MCMTpy.sampler.pyfk_MH.sampler_module import get_MISFIT_1,get_MISFIT_2,inv_output_file








#%%########################################################################
#                           -------------------
#                            1. main-function
#                           -------------------
###########################################################################w

def sample_MH(filename):
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #########################################################################
    #                         -----------------
    #                           1.1. preprocess
    #                         -----------------
    #########################################################################

    if rank == 0:
        try:
            # read the inversion superparameter dict, station parameter dict, original data dict, source time function (obspy. Trace)
            Inv_para,Sta_inv_info,Sta_raw_data,Stf = read_Inv_json(filename)
            print("Setting: read inversion.json successful in path",filename,flush=True)
        except Exception as inst:
            print(inst,flush=True)
            raise ValueError('read_Inv_json error')

        try:
            # read gfs database info
            Prepro_para,_,Source_Station_info_MPI = read_GFs_json(Inv_para["GFs_json_file"])
            print("Setting: read built_GFs.json successful in path",Inv_para["GFs_json_file"],flush=True)
        except Exception as inst:
            print(inst,flush=True)
            raise ValueError('read_GFs_json error')


        MPI_n = Inv_para['MPI_n']
        Chains_n = Inv_para['Chains_n']


        #------------ a. make directory ------------#
        if os.path.exists(Inv_para["Output_path"]):
            shutil.rmtree(Inv_para["Output_path"])
        if os.path.exists(Inv_para["Output_path"]) == False:
            os.makedirs(Inv_para["Output_path"])
            print("Setting: mkdir Output_path successful in path ",Inv_para["Output_path"],flush=True)


        #------------ b. output parameter info ------------#
        file_path = os.path.join(Inv_para["Output_path"],'Inv_para_info.txt') 
        write_inv_para(file_path,Inv_para,Sta_inv_info,Sta_raw_data,Stf)
        print("Setting: write Inv_para_info.txt successful in path ",Inv_para["Output_path"],flush=True)
        print("\n*****************************************************************",flush=True)
        print(f"Now begin sampling {Chains_n} Markov Chains with {MPI_n} cores ...\n",flush=True)


    else:
        Inv_para,Sta_inv_info,Sta_raw_data,Stf,Prepro_para,Source_Station_info_MPI = [None for _ in range(6)]



    #########################################################################
    #                         ------------------------------
    #                          1.2. broadcast the variables
    #                         ------------------------------
    #########################################################################
    Inv_para = comm.bcast(Inv_para,root=0)
    Sta_inv_info = comm.bcast(Sta_inv_info,root=0)
    Sta_raw_data = comm.bcast(Sta_raw_data,root=0)
    Stf = comm.bcast(Stf,root=0)
    Prepro_para = comm.bcast(Prepro_para,root=0)
    Source_Station_info_MPI = comm.bcast(Source_Station_info_MPI,root=0)



    #########################################################################
    #                     -----------------------------------------
    #                     1.3. MPI MCMC_inv loop through each chunk
    #                     -----------------------------------------
    #########################################################################
    try:
        MCMC_inv(Inv_para,
                 Sta_inv_info,
                 Sta_raw_data,
                 Stf,
                 Prepro_para,
                 Source_Station_info_MPI)
    except Exception as inst:
        print(inst,flush=True)
        raise ValueError('MCMC_inv error')



    #########################################################################
    #                         -----------------------
    #                          1.4. MPI comm.barrier
    #                         -----------------------
    #########################################################################
    comm.barrier()
    if rank == 0:
        print("\n\n*****************************************************************",flush=True)
        print("Successful !\nThe inversion process is completed.\n\n",flush=True)
        sys.exit()







#%%########################################################################
#                           -------------------
#                            2. sub-function
#                           -------------------
###########################################################################

def MCMC_inv(Inv_para,Sta_inv_info,Sta_raw_data,Stf,Prepro_para,Source_Station_info_MPI):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    #########################################################################
    #                         -----------------
    #                           2.1. preprocess
    #                         -----------------
    #########################################################################
    Dt = Inv_para["Dt"]                                                             # sampling rate of the raw data
    dt = Prepro_para["dt"]                                                          # sampling rate of the gfs database

    if Dt!=dt:
        raise ValueError('Error!\nGFs sample rate is not consistent with Raw_data!')

    Chains_n = Inv_para["Chains_n"]
    N=Inv_para["N"]
    com_noise_wave=1/(Inv_para["Noise_level_waveform"]**2)
    com_noise_time=1/(Inv_para["Noise_level_phasetime"]**2)
    InvType=Inv_para["InvType"]
    FM0=Inv_para["FM0"]
    N_k=Inv_para["N_k"]
    N_mag=Inv_para["N_mag"]
    FM_boundary = Inv_para["FM_boundary"]

    source_prem,_,config_prem=config_fk(InvType)
    source_time_function = Stf





    #########################################################################
    #                         -----------------
    #                           2.2. inv
    #                         -----------------
    #########################################################################
    
    for chain in range(0,Chains_n,1):
        #%% 2.2.1 process
        FM_1_all=[];FM_1_accept_all=[]
        FM_2_all=[];FM_2_accept_all=[]
        MISFIT_1_all=[];MISFIT_1_accept_all=[];
        MISFIT_2_all=[];MISFIT_2_accept_all=[];
        accpet_ratio_all=[]
        
        FM=FM0
        MISFIT_1 = 0
        MISFIT_2 = 0
        accept_num=0                                                            # The number of accepted solutions
        generate_num=0                                                          # The number of solutions generated


        if rank==0:
            bar = tqdm(range(0,round(N),1),position=rank,file=sys.stdout)
            print("\tprint rank 0 as example:",flush=True)
        else:
            bar = range(0,round(N),1)

        #%% 2.2.2 inv process
        for jj in bar:
            #------------ a. genrate new FM ------------#
            sigma = get_sigma(Inv_para,jj)                                      # The standard deviation of randomly generating new FM
            Fixed_FM_used = get_fixed_fm(Inv_para,jj)
            FM_new=[]                                                           # [m0, str, dip rake, lat, lon, depth, t0] in dc
            
            for i in range(0,len(FM),1):
                if Fixed_FM_used[i]=='constant':
                    fm_new = FM[i]
                else:
                    fm_new = np.random.normal(loc=FM[i],scale=sigma[i],size=1)[0] 
                
                if Fixed_FM_used[0]=='variable' and i==0 and (jj>=N_k) and jj<(N_k+N_mag):
                    if ('Mw_mean' in locals().keys()):
                        fm_new = Mw_mean+FM[0]
                    else:
                        fm_new = FM[0]
                FM_new.append(fm_new)
            
            FM_new_used = FM_new.copy()
            if InvType=='dc':
                if FM_new[1]>360 and FM_new[1]<720:
                    FM_new_used[1]=FM_new_used[1]-360
                if FM_new[1]<0 and FM_new[1]>(-360):
                    FM_new_used[1]=FM_new_used[1]+360

                if FM_new[3]>180 and FM_new[3]<540:
                    FM_new_used[3]=FM_new_used[3]-360
                if FM_new[3]<(-180) and FM_new[3]>(-540):
                    FM_new_used[3]=FM_new_used[3]+360

            #------------ b. the flag is used to determine whether the new solution is beyond the boundary ------------#
            Is_Beyond = False
            for i in range(0,len(FM),1):
                if (FM_new_used[i]<FM_boundary[i][0])or(FM_new_used[i]>FM_boundary[i][1]):
                    Is_Beyond = True
                    break


            #------------ c.forward model ------------#
            if Is_Beyond==True:
                print(jj+1,'th: Beyond\n')
                continue
                # jj-=1
            else:
                # the number of accepted solutions
                generate_num = generate_num+1
                #------------ % c.1.forward model || gets info of all stations for a given single source ------------#
                GFs_sta_info = get_GFs_sta_info(FM_new_used[-4:-1],Inv_para["NET_STA"],Prepro_para,Source_Station_info_MPI)


                #------------ % c.2.forward model || give FM and stf ------------#
                if InvType=='mt':
                    # a):
                    #   The seismic moment given by the JSON file is MW magnitude.
                    #   pyfk uses dyne-cm units for the synthesis of the moment tensor MT solution, which is converted to M0
                    source_mechanism = np.array(FM_new_used[0:7])
                    mw = source_mechanism[0]
                    source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
                elif InvType=='dc':
                    # b):
                    #   when == dc, pyfk uses moment magnitude, which does not need to be changed
                    source_mechanism = np.array(FM_new_used[0:4])
                elif InvType=='sf':
                    # c):
                    #   The seismic moment given by the JSON file is MW magnitude.
                    #   pyfk uses dyne units when synthesing the moment tensor sf solution, which is converted to M0
                    source_mechanism = np.array(FM_new_used[0:3])
                    mw = source_mechanism[0]
                    source_mechanism[0] = np.power(10., 1.5 * mw + 16.1)
                
                source_prem.update_source_mechanism(source_mechanism)


                #------------ % c.3.forward model syn，and stores all the stations data in a dictionary, and reads the TP Ts ------------#
                Sta_rand_data = get_Sta_rand_data(GFs_sta_info,config_prem,source_time_function)


                #------------ % c.4.MISFIT 1 ------------#
                if jj<N_k:
                    misfit_1 = get_MISFIT_1(FM_new_used,Sta_raw_data,Sta_rand_data)
                    MISFIT_1 = com_noise_time*misfit_1
                    MISFIT_rand = com_noise_time*misfit_1
                    MISFIT_1_all.append(MISFIT_rand)
                    FM_1_all.append(FM_new)
                else:
                #------------ % c.5.MISFIT 2 ------------#
                    misfit_2,Mw_mean = get_MISFIT_2(GFs_sta_info,Sta_raw_data,Sta_rand_data,Sta_inv_info,Inv_para,Prepro_para,FM_new_used[-1])
                    MISFIT_2 = com_noise_wave*misfit_2
                    MISFIT_rand = com_noise_wave*misfit_2
                    MISFIT_2_all.append(MISFIT_rand)
                    FM_2_all.append(FM_new)
                    
                #------------ % c.6. MH algorithm ------------#
                if not ('MISFIT' in locals().keys()):
                    MISFIT = MISFIT_rand

                prob_hr=(- (MISFIT_rand - MISFIT)/2)
                if  prob_hr >= 709:
                    hr=1
                elif jj>=N_k and jj<N_k+N_mag:
                    hr=1
                else:
                    hr=min(1,math.exp(prob_hr))

                r=random.random() 
                if  r<hr:
                    FM=FM_new
                    MISFIT=MISFIT_rand
                    if jj<N_k:
                        MISFIT_1_accept_all.append(MISFIT_rand)
                        FM_1_accept_all.append(FM_new)
                    else:
                        MISFIT_2_accept_all.append(MISFIT_rand)
                        FM_2_accept_all.append(FM_new)
                    accept_num = accept_num+1
                    accept_or_reject = 'accept!'
                else:
                    FM=FM
                    MISFIT=MISFIT
                    accept_or_reject = 'reject!'


                #------------ c.7 rank 0 display ------------#
                # get: MISFIT_1_min FM_1_min_used  || MISFIT_2_min FM_2_min_used
                if jj<N_k:
                    MISFIT_1_min=min(MISFIT_1_all)
                    index_min=MISFIT_1_all.index(min(MISFIT_1_all))
                    FM_1_min=FM_1_all[index_min]
                    FM_1_min_used = FM_1_min.copy()
                else:
                    MISFIT_2_min=min(MISFIT_2_all)
                    index_min=MISFIT_2_all.index(min(MISFIT_2_all))
                    FM_2_min=FM_2_all[index_min]
                    FM_2_min_used = FM_2_min.copy()
                    
                    if InvType=='dc':
                        if FM_2_min[1]>360 and FM_2_min[1]<720:
                            FM_2_min_used[1]=FM_2_min_used[1]-360
                        if FM_2_min[1]<0 and FM_2_min[1]>(-360):
                            FM_2_min_used[1]=FM_2_min_used[1]+360

                        if FM_2_min[3]>180 and FM_2_min[3]<540:
                            FM_2_min_used[3]=FM_2_min_used[3]-360
                        if FM_2_min[3]<(-180) and FM_2_min[3]>(-540):
                            FM_2_min_used[3]=FM_2_min_used[3]+360

                # get accpet_ratio
                accpet_ratio=round(accept_num/generate_num,5)
                accpet_ratio_all.append(accpet_ratio)
                
                
                # get recommand Noise_level_phasetime and Noise_level_waveform
                if rank == 0:
                    # stage 1: get diff_MISFIT_1 diff_MISFIT_1_mean
                    # get difference between 2 sample of wg_misfit
                    if len(MISFIT_1_all) > 1:
                        diff_MISFIT_1 = MISFIT_1 - MISFIT_1_all[-2]
                    else:
                        diff_MISFIT_1 = 0
                    # get mean difference between 2 sample of wg_misfit
                    if len(MISFIT_1_all) < 100:
                        diff_MISFIT_1_mean = np.mean(abs(np.diff(MISFIT_1_all)))/com_noise_time
                    else:
                        diff_MISFIT_1_mean = np.mean(abs(np.diff(np.array(MISFIT_1_all)[-100:])))/com_noise_time

                    # stage 2: get diff_MISFIT_2 diff_MISFIT_2_mean
                    # get difference between 2 sample of wg_misfit
                    if len(MISFIT_2_all) > 1:
                        diff_MISFIT_2 = MISFIT_2 - MISFIT_2_all[-2]
                    else:
                        diff_MISFIT_2 = 0
                    # get mean difference between 2 sample of wg_misfit
                    if len(MISFIT_2_all) < 100:
                        diff_MISFIT_2_mean = np.mean(abs(np.diff(MISFIT_2_all)))/com_noise_wave
                    else:
                        diff_MISFIT_2_mean = np.mean(abs(np.diff(np.array(MISFIT_2_all)[-100:])))/com_noise_wave
                        
                    
                    remin_Noise_level_phasetime = round(math.sqrt(abs(diff_MISFIT_1_mean/math.log(0.2))/2),5)
                    remax_Noise_level_phasetime = round(math.sqrt(abs(diff_MISFIT_1_mean/math.log(0.5))/2),5)
                    remin_Noise_level_waveform = round(math.sqrt(abs(diff_MISFIT_2_mean/math.log(0.2))/2),5)
                    remax_Noise_level_waveform = round(math.sqrt(abs(diff_MISFIT_2_mean/math.log(0.5))/2),5)
                    
                    
                    bar.write(f"{com_noise_time},{com_noise_wave}")
                    
                    bar.write(f"Iter {jj+1} of {chain}th Markov-Chain: ")
                    # MISFIT1 MISFIT2
                    bar.write(f'MISFIT1:{round(MISFIT_1,5)} || min MISFIT1:{round(min(MISFIT_1_all),5)}')
                    if jj<N_k:
                        bar.write('MISFIT2:no || min MISFIT2:no')
                    else:
                        bar.write(f'MISFIT2:{round(MISFIT_2,5)} || min MISFIT2:{round(min(MISFIT_2_all),5)}')
                    
                    # mean MISFIT1 MISFIT2
                    bar.write(f"mean diff_MISFIT1: {round(diff_MISFIT_1_mean,5)} \
|| Recommended Noise_level_phasetime: [min|max] \
{remin_Noise_level_phasetime}|{remax_Noise_level_phasetime}")
                    if jj<N_k:
                        bar.write(f"mean diff_MISFIT2: {round(diff_MISFIT_2_mean,5)} \
|| Recommended Noise_level_waveform: [min|max] \
{remin_Noise_level_waveform}|{remax_Noise_level_waveform}")
                    else:
                        bar.write(f"mean diff_MISFIT2: {round(diff_MISFIT_2_mean,5)} \
|| Recommended Noise_level_waveform: min|max: \
{remin_Noise_level_waveform}|{remax_Noise_level_waveform}")
                    # 
                    bar.write(f"MISFIT_rand:{round(MISFIT_rand,5)}")
                    bar.write(f"accept_ratio:{accpet_ratio}")
                    bar.write(f"prob_hr:{round(prob_hr,3)}")                    # debug:
                    bar.write(f"hr:{round(hr,3)}")
                    
                    fm=''
                    for rr in range(0,len(FM_new_used)):
                        fm+=str(round(FM_new_used[rr],4))+'  '
                    bar.write('FM_rand: '+str(fm))
                    
                    fm1=''
                    for rr in range(0,len(FM_1_min_used)):
                        fm1+=str(round(FM_1_min_used[rr],4))+'  '
                    bar.write('FM_1_min: '+str(fm1))
                    
                    if jj<N_k:
                        bar.write('FM_2_min: no')
                    else:
                        fm2=''
                        for rr in range(0,len(FM_2_min_used)):
                            fm2+=str(round(FM_2_min_used[rr],4))+'  '
                        bar.write('FM_2_min: '+str(fm2))
                    
                    bar.write('Accept_or_Reject: '+accept_or_reject+'\n')
                    bar.set_description(f"Rank {rank}: now MISFIT_weight={round(MISFIT,5)}") 






        #%% 2.2.3 output
        file_path = os.path.join(Inv_para["Output_path"],'rank'+'_'+str(rank)+'_output')
        if not os.path.isdir(file_path):os.mkdir(file_path)

        inv_output_file(file_path,chain,
                        FM_1_all,FM_1_accept_all,
                        FM_2_all,FM_2_accept_all,
                        MISFIT_1_all,MISFIT_2_all,
                        MISFIT_1_accept_all,MISFIT_2_accept_all,
                        accpet_ratio_all)




