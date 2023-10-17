#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 28.09.2021

Created on Tue Jun 29 15:51:21 2021

Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
Description: extract a subset of single agent data according to the timepoints saved in df_dict_3
"""
#%%
import os
import pickle
import pandas as pd
import time
#%% load db_model to determine k_net_max and Time of k_net_max for every it

folder_names16=['Simulations_v16_PW_Const_MIC-VI_1B',
                'Simulations_v16_PW_Const_MIC+VI_1B',
                'Simulations_v16_PW_Const_G-VI_DrugIntStrengths_1B',
                'Simulations_v16_PW_Const_G+VI_DrugIntStrengths_1B']

folder_names17=['Simulations_v17_PW_Const_MIC-VI_1B',
                'Simulations_v17_PW_Const_MIC+VI_1B',
                'Simulations_v17_PW_Const_G-VI_DrugIntStrengths_1B',
                'Simulations_v17_PW_Const_G+VI_DrugIntStrengths_1B']
folder_name_range=folder_names16+folder_names17

path_core=os.getcwd() + "/Data/"
    
#%% load df_dict and db_agents, keep just the Time data
# save as db_agents_part
# combine db_agents_part from all folders
# calculate measures

#load data from all simulaton runs, save agent data of one timepoint for all iterattions in df_agents_part (dataframe) 
df_agents_part=pd.DataFrame(columns=['Species','Target','timestep']+['drug_type',"drug","n_steps","int","int_radius", "int_strength" ,"move_type","concentration","iteration"])

#load and extract timepoint data from v16 folders
# for folder_name in folder_name_range:
print('df_agents_part')
for folder_name in folder_names16:
    t0 = time.time() # start timer
    path = path_core + folder_name
    file = open(path + "/db_agents.pkl", "rb")
    db_agents = pickle.load(file)
    file = open(path + "/df_dict_3.pkl", "rb")
    df_dict = pickle.load(file)
    key=list(db_agents.keys())[0]
    for key in db_agents:
        tp=int(df_dict[df_dict['drug']==key[0]][df_dict['n_steps']==key[1]][df_dict['int']==key[2]][df_dict['int_radius']==key[3]][df_dict['int_strength']==key[4]][df_dict['move_type']==key[5]][df_dict['concentration']==key[6]][df_dict['iteration']==key[7]]['NetGrowthTimestep'])
        df=db_agents[key].loc[tp].copy()
        df=df[df['Species']=='A']
        df=df.drop(['Growthrate','Killrate','NetGrowthrate','Time','AgentID'],axis=1)
        #rename MIC_I or Max_replicationrate column as target
        if 'MIC_I' in df.columns:
            df.rename(columns={'MIC_I':'Target'},inplace=True)
        elif 'Max_replicationrate' in df.columns:
            df.rename(columns={'Max_replicationrate':'Target'},inplace=True)
        cols_to_add=["drug","n_steps","int","int_radius", "int_strength" ,"move_type","concentration","iteration"]
        for i in range(len(cols_to_add)):
            df[cols_to_add[i]]=key[i]
        df['drug_type']='add'
        df_agents_part=df_agents_part.append(df).reset_index(drop=True)
        del df
    del db_agents
    del df_dict
    print(folder_name)
    t1=time.time()
    print(t1-t0)   
    
# load and extract timepoint data from v17 folders
# for folder_name in folder_name_range:
for folder_name in folder_names17:
    t0 = time.time() # start timer
    path = path_core + folder_name
    file = open(path + "/db_agents.pkl", "rb")
    db_agents = pickle.load(file)
    file = open(path + "/df_dict_3.pkl", "rb")
    df_dict = pickle.load(file)
    key=list(db_agents.keys())[0]
    for key in db_agents:
        t=int(df_dict[df_dict['drug']==key[0]][df_dict['n_steps']==key[1]][df_dict['int']==key[2]][df_dict['int_radius']==key[3]][df_dict['int_strength']==key[4]][df_dict['move_type']==key[5]][df_dict['concentration']==key[6]][df_dict['iteration']==key[7]]['NetGrowthTimestep'])
        df=db_agents[key][db_agents[key]['Time']==t][db_agents[key]['Species']=='A'].copy()
        df=df.drop(['Growthrate','Killrate','NetGrowthrate','timestep','AgentID'],axis=1)
        cols_to_add=["drug","n_steps","int","int_radius", "int_strength" ,"move_type","concentration","iteration"]
        for i in range(len(cols_to_add)):
            df[cols_to_add[i]]=key[i]
        df['drug_type']='prop'
        df_agents_part=df_agents_part.append(df).reset_index(drop=True)
        del df
    del db_agents
    del df_dict 
    print(folder_name)
    t1=time.time()
    print(t1-t0) 
    
#save df_agents_part 
file = open(path_core + "/df_agents_part_3.pkl", "wb")
pickle.dump(df_agents_part, file)
file.close()
# del df_agents_part
print('metrics')
t0 = time.time() # start timer
# calculate mean, std and cv for each iteration and save in df_agents_part_metrics
df_agents_part['Target']=df_agents_part['Target'].astype(float)
df_agents_part['concentration']=df_agents_part['concentration'].astype(float)

df_metrics=df_agents_part.groupby(['Species', 'Time', 'drug_type', 'drug', 'n_steps', 'int', 'int_radius', 'int_strength', 'move_type', 'concentration','iteration'])[['Target']].describe()
df_metrics.reset_index(inplace=True)

df_metrics[('Target','CV')]=df_metrics[('Target','std')]/df_metrics[('Target','mean')]
file = open(path_core + "/df_metrics_3.pkl", "wb")
pickle.dump(df_metrics, file)
file.close()