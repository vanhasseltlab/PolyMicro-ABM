#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 28.09.2021

Created on Tue Jun 29 15:51:21 2021

Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
Description: calculates the maximal net growth rate only including datapoints were at least 10 agents are alive, saved in df_dict_3.pkl for each simulation set
"""
#%%
import os
import pickle
import numpy as np
import pandas as pd
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
folder_name='Simulations_v16_PW_Const_MIC-VI_1B'

path_core=os.getcwd() + "/Data/"
#load data from multiple simulaton runs, calculate max k_net Time, save in df_dict
print('df_dict')
for folder_name in folder_name_range:
    path = path_core + folder_name
    file = open(path + "/db_model.pkl", "rb")
    db_model = pickle.load(file)
    file.close()
    
    # make dataframe with max k_nets and auc from db_model, save as dataframe
    data=[]
    columns=["drug","n_steps","int","int_radius", "int_strength" ,"move_type","concentration","iteration","NetGrowth","AUC","NetGrowthTimestep"]
    for key in db_model:
        i=np.nanargmax(db_model[key][db_model[key]["#A"]>9]['Max_Av_k_net_A'].abs())
        k_nets=db_model[key]['Max_Av_k_net_A'].loc[i]
        t=db_model[key]['timestep'].loc[i]
        auc=np.nansum(db_model[key]['#A'])
        d = [key[0],key[1],key[2],key[3],key[4],key[5],key[6],key[7],k_nets,auc,t]
        data.append(d)
    del db_model
    df_dict = pd.DataFrame(data, columns=columns)
    del data

    #save df_dicts in pickle file to be able to retrieve timepoints for calculation of k_nets
    path = path_core + folder_name
    file = open(path + "/df_dict_3.pkl", "wb")
    pickle.dump(df_dict, file)
    file.close()
    del df_dict
    print(folder_name)

