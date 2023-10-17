#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Date: 28.09.2021

Created on Tue Jun 29 15:51:21 2021

Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
Description: creates combined dataframe of db_model data for all v17 simulation sets & combines combined dataframes for v16 and v17 into one, saves in df_all_v16&v17.pkl

"""
#%%
#%reset -f
#%matplotlib qt
#%%
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import time
# working directory: '/home/cathi/Simulations/Software Supplementary'

# define plot formatting settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.figsize'] = [3.5,3.5]
mpl.rcParams['figure.dpi']= 300
mpl.rcParams['lines.linewidth']  = 1
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['legend.title_fontsize'] = 8
mpl.rcParams['axes.labelsize']=8
mpl.rcParams['axes.titlesize']=8
mpl.rcParams['xtick.labelsize']=8 
mpl.rcParams['ytick.labelsize']=8
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

folder_name0="Simulations_v17_PW_Const_VI"
folder_name1="Simulations_v17_PW_Const_MIC-VI_1B"
folder_name2="Simulations_v17_PW_Const_G-VI_DrugIntStrengths_1B"
folder_name3="Simulations_v17_PW_Const_MIC+VI_1B"
folder_name4="Simulations_v17_PW_Const_G+VI_DrugIntStrengths_1B"

folder_name_range=[folder_name0,folder_name1,folder_name2,folder_name3,folder_name4]

df_all=pd.DataFrame(columns=["drug","int","int_radius", "int_strength" ,"move_type","n_steps","concentration", "Av-NetGrowth","Av-AUC","Std-NetGrowth","Std-AUC"])
#load data from multiple simulaton runs, calculate max k_net and auc and save in one dataframe df_all
for folder_name in folder_name_range:
    path = os.getcwd() + "/Data" + folder_name    
    file = open(path + "/db_model.pkl", "rb")
    db_model = pickle.load(file)
    

    data=[]
    columns=["drug","int","int_radius", "int_strength" ,"move_type","n_steps","concentration","iteration", "NetGrowthMean","AUC","NetGrowthSum",'TimestepMean','TimestepSum']
    for key in db_model:
        i=np.nanargmax(db_model[key]['Max_Av_k_net_A'].abs())
        k_nets = db_model[key]['Max_Av_k_net_A'].loc[i]
        tp_mean=db_model[key]['timestep'].loc[i]
        
        i=np.nanargmax(db_model[key]['Max_Av_k_net_A'].abs()*db_model[key]['#A'])
        k_nets_sum=db_model[key]['Max_Av_k_net_A'].loc[i]*db_model[key]['#A'].loc[i]
        tp_sum=db_model[key]['timestep'].loc[i]
        
        auc=np.nansum(db_model[key]['#A'])
        
        d = [key[0],key[2],key[3],key[4],key[5],key[1],key[6],key[7],k_nets,auc,k_nets_sum,tp_mean,tp_sum]
        data.append(d)
    del db_model
    df_dict = pd.DataFrame(data, columns=columns) 
    del data
    
    #%
    #calculate mean and std across all iterations, save as dataframe
    #
    rows=df_dict.drop_duplicates(subset=["drug","int","int_radius", "int_strength" ,"move_type","n_steps","concentration"],keep='first')[["drug","int","int_radius", "int_strength" ,"move_type","n_steps","concentration"]]
    columns=["drug","int","int_radius", "int_strength" ,"move_type","n_steps","concentration", "Av-NetGrowthMean","Av-AUC","Av-NetGrowthSum",'Av-TpMean','Av-TpSum',"Std-NetGrowthMean","Std-AUC","Std-NetGrowthSum",'Std-TpMean','Std-TpSum']
    data=[]
    for index,r in rows.iterrows():
        #calculate mean and std, save in d
        df=df_dict[df_dict['drug']==r[0]][df_dict['int']==r[1]][df_dict['int_radius']==r[2]][df_dict['int_strength']==r[3]][df_dict['move_type']==r[4]][df_dict['n_steps']==r[5]][df_dict['concentration']==r[6]]
        d = list(r) + [np.nanmean(df['NetGrowthMean']),np.nanmean(df['AUC']),np.nanmean(df['NetGrowthSum']),np.nanmean(df['TimestepMean']),np.nanmean(df['TimestepSum']),np.nanstd(df['NetGrowthMean']),np.nanstd(df['AUC']),np.nanstd(df['NetGrowthSum']),np.nanstd(df['TimestepMean']),np.nanstd(df['TimestepSum'])]
        data.append(d)
    del df_dict
    df_dict_summary = pd.DataFrame(data, columns=columns) 
    del data
    
    df_all=df_all.append(df_dict_summary) #dataframe with data from all ints together
    del df_dict_summary
    print(folder_name)
df_all=df_all.reset_index(drop=True)

#correct Int names | only necessary for no_limit v16 simulations
for index,row in df_all.iterrows():
    if row['int']=='MIC-' and row['int_strength']>0:
        df_all.at[index,'int']= 'MIC+'
    if row['int']=='G-' and row['int_strength']>0:
        df_all.at[index,'int']= 'G+'


#replace conc value with new names in df_all
drugs=list(dict.fromkeys(df_all['drug']))
#for rows of each drug, replace names separately  in extra column
for drug in drugs:
    df = df_all[df_all['drug']==drug]
    old_names=sorted(list(dict.fromkeys(df['concentration'])))
    new_names=list(range(len(old_names)))
    dict_names=dict(zip(old_names, new_names))
    for index,row in df.iterrows():
            df_all.at[index,'concentration_simple']=dict_names[df_all.loc[index]['concentration']]

#replace int strengths with new names in df_all 
drugs=list(dict.fromkeys(df_all['drug']))
int_types=list(dict.fromkeys(df_all['int']))
int_rad=list(dict.fromkeys(df_all['int_radius']))

#for rows of each drug, replace names separately in extra column
for drug in drugs:
    for int_t in int_types:
        for int_r in int_rad:
            df = df_all[df_all['drug']==drug]
            df = df[df['int']==int_t]
            df = df[df['int_radius']==int_r]
            old_names=sorted(list(dict.fromkeys(df['int_strength'])),key=abs)
            new_names=['weak','medium','strong']
            dict_names=dict(zip(old_names, new_names))
            for index,row in df.iterrows():
                    df_all.at[index,'int_strength_string']=dict_names[df_all.loc[index]['int_strength']]
            del df

#replace random with r and stepwise with s
for index,row in df_all.iterrows():
    if df_all.at[index,'move_type']=='stepwise':
        df_all.at[index,'move_type_simple']='s'
    if df_all.at[index,'move_type']=='random':
        df_all.at[index,'move_type_simple']='r'
        
# define categoricals
df_all['drug'] = pd.Categorical(df_all.drug,categories=[5.0,1.75,4.0,2.0],ordered=True)
df_all['int_strength_string'] = pd.Categorical(df_all.int_strength_string,categories=['weak','medium','strong'],ordered=True)

#%% calculate relative change for AUC

col_RC = []

t0 = time.time() # start timer

for index,row in df_all.iterrows():
    
    row_c_x_AUC = df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==row["concentration_simple"]][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-AUC"]
    row_c_last_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==8][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-AUC"])[0]
    row_c_first_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==0][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-AUC"])[0]
                                                                 
    RC = list((row["Av-AUC"]-row_c_x_AUC)/abs(row_c_last_AUC-row_c_first_AUC))[0]
    # (v2-c2)/abs(c4-c0)
    df_all.at[index,'RC-AUC']=RC     
    

t1=time.time()
print(t1-t0)

#%% calculate relative change for k_net

col_RC = []

t0 = time.time() # start timer

for index,row in df_all.iterrows():
    
    row_c_x_AUC = df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==row["concentration_simple"]][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthMean"]
    row_c_last_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==8][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthMean"])[0]
    row_c_first_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==0][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthMean"])[0]
                                                                 
    RC = list((row["Av-NetGrowthMean"]-row_c_x_AUC)/abs(row_c_last_AUC-row_c_first_AUC))[0]
    # (v2-c2)/abs(c4-c0)
    df_all.at[index,'RC-NetGrowthMean']=RC     
    
t1=time.time()
print(t1-t0)

#%% calculate relative change for k_net sum

col_RC = []

t0 = time.time() # start timer

for index,row in df_all.iterrows():
    
    row_c_x_AUC = df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==row["concentration_simple"]][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthSum"]
    row_c_last_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==8][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthSum"])[0]
    row_c_first_AUC = list(df_all[df_all["int"]==False][df_all["drug"]==row["drug"]][df_all["concentration_simple"]==0][df_all["n_steps"]==row["n_steps"]][df_all["move_type"]==row["move_type"]]["Av-NetGrowthSum"])[0]
                                                                 
    RC = list((row["Av-NetGrowthSum"]-row_c_x_AUC)/abs(row_c_last_AUC-row_c_first_AUC))[0]
    # (v2-c2)/abs(c4-c0)
    df_all.at[index,'RC-NetGrowthSum']=RC     
    
t1=time.time()
print(t1-t0)
#%% save df_all in pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_DrugIntStrengths_v17_2_VI.pkl", "wb")
pickle.dump(df_all, file)
file.close()  


#%% load df_all from two sets of simulations v16 and v17 and combine
#combine df_all from Simv16 and Simv17 into one
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_DrugIntStrengths_2_VI.pkl", "rb")
df_all1 = pickle.load(file)

# add column drug type
df_all1['drug_type']='prop'

path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_DrugIntStrengths_v17_2_VI.pkl", "rb")
df_all2 = pickle.load(file)

# add column drug type
df_all2['drug_type']='add'
    
df_all=df_all1.append(df_all2) #dataframe with data from all ints together
df_all=df_all.reset_index(drop=True)
df_all['drug_type'] = pd.Categorical(df_all.drug_type,categories=['prop','add'],ordered=True)


for index,row in df_all.iterrows():
    if row['drug']==5 or row['drug']==1.75:
        df_all.at[index,'kill_type']='cidal'
    elif row['drug']==4 or row['drug']==2:
        df_all.at[index,'kill_type']='static'
    if row['drug']==5 or row['drug']==4:
        df_all.at[index,'Cdep']='ci'
    elif row['drug']==1.75 or row['drug']==2:
        df_all.at[index,'Cdep']='cd'
del df_all1
del df_all2

# save df_all in pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_v16&v17.pkl", "wb")
pickle.dump(df_all, file)
file.close()  
