1/6*24*60#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 28.09.2021

Created on Tue Jun 29 15:51:21 2021

Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
Description: compares the relative interaction effect in different conditions to each other using a log2 comparison
"""

#%%
import os
import pickle
import numpy as np
import time
import pandas as pd    
import matplotlib as mpl
# working directory: '/home/cathi/Simulations/Software Supplementary'

# define plot formatting settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.figsize'] = [3.5,3.5]
mpl.rcParams['figure.dpi']= 300
mpl.rcParams['lines.linewidth']  = 1
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.title_fontsize'] = 8
mpl.rcParams['axes.labelsize']=6
mpl.rcParams['axes.titlesize']=6
mpl.rcParams['xtick.labelsize']=6
mpl.rcParams['ytick.labelsize']=6
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 0.5
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 0.5
mpl.rcParams['ytick.minor.width'] = 0.3


#%% 
#first create or load df_all from Create_load_df_all.py
# if df_all_s exists, jumpt to line 215 to directly load that one
#%% load df_all from pickle file
path_all = os.getcwd() + "/Simulations/Simu_PW_v17/"
file = open(path_all + "/df_all_v16&v17.pkl", "rb")
df_all = pickle.load(file)

#%% compare no movement (stepwise, nsteps=0) to movement (stepwise,nsteps=1) using log2 for AUC
#and compare stepwise relative change to random relative change using abs difference

t0 = time.time() # start timer

df_all_s = df_all[df_all["int"]!=False][df_all["move_type"]=="stepwise"]
df_all_s = df_all_s[df_all_s["n_steps"]==0]

#order columns correctly
df_all_s['drug'] = pd.Categorical(df_all_s.drug,categories=[5.0,1.75,4.0,2.0],ordered=True)
df_all_s['int'] = pd.Categorical(df_all_s.int,categories=["MIC-","G-","MIC+","G+"],ordered=True)

#add column for log2 comparison of rel change
col_diff = [np.nan]*len(df_all_s)
df_all_s["LC-AUC"]=col_diff # log2change of RC
df_all_s["AD-AUC"]=col_diff # abs difference of rc auc

for index,row in df_all_s.iterrows():
    r1=row["RC-AUC"]
    
    #rel change of movement counterpart   
    r0=list(df_all[df_all["int"]!=False][df_all["drug"]==row['drug']][df_all["drug_type"]==row['drug_type']][df_all["int"]==row['int']][df_all["int_radius"]==row['int_radius']][df_all["int_strength"]==row['int_strength']][df_all["move_type"]=='stepwise'][df_all["n_steps"]==1][df_all["concentration"]==row["concentration"]]["RC-AUC"])[0] 


    if r0 != 0 and r1/r0>0:
        df_all_s.loc[index,"LC-AUC"]= np.log2(r1/r0)
    else:
        df_all_s.loc[index,"LC-AUC"]= np.nan
    
    df_all_s.loc[index,"AD-AUC"]= abs(r1-r0) # |RC_AUC_SW - RC_AUC_R|

t1=time.time()
print(t1-t0)

#% compare stepwise relative change to random relative change using log2 for k_net
#and compare stepwise relative change to random relative change using abs difference

t0 = time.time() # start timer

#add column for log2 comparison of rel change
col_diff = [np.nan]*len(df_all_s)
df_all_s["LC-NetGrowth"]=col_diff # log2change of RC
df_all_s["AD-NetGrowth"]=col_diff # abs difference of rc auc

for index,row in df_all_s.iterrows():
    r1=row["RC-NetGrowth"]
    
    #rel change of movement counterpart   
    r0=list(df_all[df_all["int"]!=False][df_all["drug"]==row['drug']][df_all["drug_type"]==row['drug_type']][df_all["int"]==row['int']][df_all["int_radius"]==row['int_radius']][df_all["int_strength"]==row['int_strength']][df_all["move_type"]=='stepwise'][df_all["n_steps"]==1][df_all["concentration"]==row["concentration"]]["RC-NetGrowth"])[0] 

       
    if r0 != 0 and r1/r0>0:
        df_all_s.loc[index,"LC-NetGrowth"]= np.log2(r1/r0)
    else:
        df_all_s.loc[index,"LC-NetGrowth"]= np.nan
    
    df_all_s.loc[index,"AD-NetGrowth"]= abs(r1-r0) # |RC_AUC_SW - RC_AUC_R|

t1=time.time()
print(t1-t0)

# save df_all_s in pickle file
path_all = os.getcwd() + "/Simulations/Simu_PW_v17/"
file = open(path_all + "/df_all_s_mvnomv_v16&v17_VI.pkl", "wb")
pickle.dump(df_all_s, file)
file.close()  

#%% compare int_rad 1 to int_rad 0
#no movement (stepwise, nsteps=0) to movement (stepwise,nsteps=1) using log2 for AUC
#and compare stepwise relative change to random relative change using abs difference

t0 = time.time() # start timer

df_all_s = df_all[df_all["int"]!=False]
df_all_s = df_all_s[df_all_s["int_radius"]==1]

#order columns correctly
df_all_s['drug'] = pd.Categorical(df_all_s.drug,categories=[5.0,1.75,4.0,2.0],ordered=True)
df_all_s['int'] = pd.Categorical(df_all_s.int,categories=["MIC-","G-","MIC+","G+"],ordered=True)

#add column for log2 comparison of rel change
col_diff = [np.nan]*len(df_all_s)
df_all_s["LC-AUC"]=col_diff # log2change of RC
df_all_s["AD-AUC"]=col_diff # abs difference of rc auc

for index,row in df_all_s.iterrows():
    r1=row["RC-AUC"]
    
    #rel change of int_rad =0 counterpart   
    r0=list(df_all[df_all["int"]!=False][df_all["drug"]==row['drug']][df_all["drug_type"]==row['drug_type']][df_all["int"]==row['int']][df_all["int_radius"]==0][df_all["int_strength_string"]==row['int_strength_string']][df_all["move_type"]==row["move_type"]][df_all["n_steps"]==row["n_steps"]][df_all["concentration"]==row["concentration"]]["RC-AUC"])[0] 

       
    if r0 != 0 and r1/r0>0:
        df_all_s.loc[index,"LC-AUC"]= np.log2(r1/r0)
    else:
        df_all_s.loc[index,"LC-AUC"]= np.nan
    
    df_all_s.loc[index,"AD-AUC"]= abs(r1-r0) # |RC_AUC_SW - RC_AUC_R|

t1=time.time()
print(t1-t0)

#% compare stepwise relative change to random relative change using log2 for k_net
#and compare stepwise relative change to random relative change using abs difference

t0 = time.time() # start timer

#add column for log2 comparison of rel change
col_diff = [np.nan]*len(df_all_s)
df_all_s["LC-NetGrowth"]=col_diff # log2change of RC
df_all_s["AD-NetGrowth"]=col_diff # abs difference of rc auc

for index,row in df_all_s.iterrows():
    r1=row["RC-NetGrowth"]
    
    #rel change of movement counterpart   
    r0=list(df_all[df_all["int"]!=False][df_all["drug"]==row['drug']][df_all["drug_type"]==row['drug_type']][df_all["int"]==row['int']][df_all["int_radius"]==0][df_all["int_strength_string"]==row['int_strength_string']][df_all["move_type"]==row["move_type"]][df_all["n_steps"]==row["n_steps"]][df_all["concentration"]==row["concentration"]]["RC-NetGrowth"])[0] 

       
    if r0 != 0 and r1/r0>0:
        df_all_s.loc[index,"LC-NetGrowth"]= np.log2(r1/r0)
    else:
        df_all_s.loc[index,"LC-NetGrowth"]= np.nan
    
    df_all_s.loc[index,"AD-NetGrowth"]= abs(r1-r0) # |RC_AUC_SW - RC_AUC_R|

t1=time.time()
print(t1-t0)

#%% save df_all_s in pickle file
path_all = os.getcwd() + "/Simulations/Simu_PW_v17/"
file = open(path_all + "/df_all_s_intrad1|0_v16&v17_VI.pkl", "wb")
pickle.dump(df_all_s, file)
file.close() 