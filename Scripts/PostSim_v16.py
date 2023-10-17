#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 16:07:24 2021
An Agent based modeling framework of an interactive, multispecies community during antimicrobial treatment
Postsimulation processing script
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
"""

# import os
import pickle
import numpy as np
import pandas as pd
import time

# ------------------------------------------------------------------------------------------------------------------                              
# ------------------------------------------------------------------------------------------------------------------               
# POST-SIMULATION PROCESSING OF RAW DATA
# ------------------------------------------------------------------------------------------------------------------                              
# ------------------------------------------------------------------------------------------------------------------               
def main(data_interval,path,db_agents,db_model,all_species,all_interactions,all_interactions_env,ab_conc_range,capacity,n_init,its,runtime,gridrange,move_type,stepsize,data_collection,variable_params,fixed_params,
         calculate_mean=[]):
#%%
    ##Determine species names for all simulated scenarios
    species_names=[x[0] for x in all_species]
  
    #Remove None columns from db_agents and db_model
    col_model = [list(db_model[key].columns) for key in db_model]
    col_agents = [list(db_agents[key].columns) for key in db_agents]

    
    if 'AB_fix_I' not in col_model[0]: # if it has not been processed already
        # drop all columns for which NONE was recorded
        for key, columns in zip(db_agents,col_agents):
            for col in columns:
                if db_agents[key][col][0,0]==None:
                    db_agents[key]=db_agents[key].drop([col],axis=1)
                    
        for key,columns in zip(db_model,col_model):
            for col in columns:
                if db_model[key][col][0]==None:
                    db_model[key]=db_model[key].drop([col],axis=1)
                    
        # add columns, do calculations before reformatting
        for key in db_agents:
            # split Number_species column in columns with specific #species_names
            if 'Number_Species' in db_model[key].columns:
                # split Number_Species column into species types
                for y in range(len(species_names)): # for every drug in the simulation
                    db_model[key] = pd.concat([db_model[key],
                               pd.DataFrame({'#'+species_names[y]: [x[y] for x in db_model[key]['Number_Species']]},
                                        index=db_model[key].index)],
                           axis=1, join='inner')
            #remove old column
            db_model[key]=db_model[key].drop(['Number_Species'],axis=1)
        
        for key in db_agents:
            # split Number_species column in columns with specific #species_names
            if 's_species' in db_model[key].columns:
                # split Number_Species column into species types
                for y in range(len(species_names)): # for every drug in the simulation
                    db_model[key] = pd.concat([db_model[key],
                               pd.DataFrame({'s_'+species_names[y]: [x[y] for x in db_model[key]['s_species']]},
                                        index=db_model[key].index)],
                           axis=1, join='inner')
            
                #remove old column
                db_model[key]=db_model[key].drop(['s_species'],axis=1)
            
            
            # split k_net column in columns with specific #species_names
            if 'Max_Av_k_net' in db_model[key].columns:
                # split Number_Species column into species types
                for y in range(len(species_names)): # for every drug in the simulation
                    db_model[key] = pd.concat([db_model[key],
                               pd.DataFrame({'Max_Av_k_net_'+species_names[y]: [x[y] for x in db_model[key]['Max_Av_k_net']]},
                                        index=db_model[key].index)],
                           axis=1, join='inner')
            
                #remove old column
                db_model[key]=db_model[key].drop(['Max_Av_k_net'],axis=1)
            
            
            if 'MIC' in db_agents[key].columns:
                # split MIC column for each drug present in simulation
                for y in range(len(ab_conc_range[0])): # for every drug in the simulation
                    db_agents[key] = pd.concat([db_agents[key],
                               pd.DataFrame({'MIC_'+"I"*(y+1): [x[y] for x in db_agents[key]["MIC"]]},
                                        index=db_agents[key].index)],
                           axis=1, join='inner') 
                               
                #remove old column
                db_agents[key]=db_agents[key].drop(['MIC'],axis=1)
            
            if 'AB' in db_model[key].columns:
            # split AB column, multiple drugs
                for y in range(len(ab_conc_range[0])): # for every drug in the simulation
                    db_model[key] = pd.concat([db_model[key],
                               pd.DataFrame({'AB_'+"I"*(y+1): [x[y] for x in db_model[key]["AB"]]},
                                        index=db_model[key].index)],
                           axis=1, join='inner')
                             
                    #remove old column
                db_model[key]=db_model[key].drop(['AB'],axis=1)
            
            if 'AB_fix' in db_model[key].columns:
    
                # split AB_fix column, multiple drugs
                for y in range(len(ab_conc_range[0])): # for every drug in the simulation
                    db_model[key] = pd.concat([db_model[key],
                               pd.DataFrame({'AB_fix_'+"I"*(y+1): [x[y] for x in db_model[key]["AB_fix"]]},
                                        index=db_model[key].index)],
                           axis=1, join='inner')
                             
                #remove old columns
                db_model[key]=db_model[key].drop(['AB_fix'],axis=1)
                db_model[key]=db_model[key].drop(['Int_current'],axis=1)
                db_model[key]=db_model[key].drop(['Int_failed'],axis=1)
                db_model[key]=db_model[key].drop(['Int_give'],axis=1)
                db_model[key]=db_model[key].drop(['Int_receive'],axis=1)
                db_model[key]=db_model[key].drop(['affected'],axis=1)
                
    
    if fixed_params['data_agents']: # only if data for agents was saved    
        # add timestep, agentID and Time column to each dataframe in agent dictionary
        for key in db_agents:
            if 'Time' not in db_agents[key].columns: # if it has not been processed already
                timestep = [x[0] for x in db_agents[key].index]
                time_list = list(db_model[key]['Time'])
                # generate TIME column for agent dataframes
                time_a = []
                j=0
                timestep_0=timestep[0]
                for i in range(len(timestep)):
                    if timestep[i]!=timestep_0:
                        timestep_0=timestep[i]
                        j+=1
                    time_a.append(time_list[j])
                
                
                # add time and timepoint columns
                db_agents[key] = pd.concat([db_agents[key],
                           pd.DataFrame({'timestep': timestep,
                                     'AgentID':[x[1] for x in db_agents[key].index],
                                     "Time": time_a} ,
                                    index=db_agents[key].index)],
                       axis=1, join='inner')
        
        # save processed version of raw data, all timepoints               
        file = open(path + "/db_agents.pkl", "wb")
        pickle.dump(db_agents, file)
        file.close()    
    
    file = open(path + "/db_model.pkl", "wb")
    pickle.dump(db_model, file)
    file.close()  
    #
    # ------------------------------------------------------------------------------------------------------------------                              
    # ------------------------------------------------------------------------------------------------------------------               
    # select part of the dataframes for calculations and count number of each species
    #%%               
    #DEFINE INTERVAL TO COLLECT SUBSET OF DATA FOR CALCULATIONS
    t0 = time.time() # start timer
    interval=data_interval# calculate species count etc every x mins, choose so that this intervals >= datacollection interval
    t_endpoints_time =  list(range(0,round(runtime),interval)) + [runtime]
    
    
    # take a subset of the dataframe at the defined andpoints
    db_m_part = db_model.copy()
    db_a_part = db_agents.copy()
    del [db_model, db_agents]
    for key in db_m_part:
    
        time_list = list(db_m_part[key]["Time"])
        Number=list(db_m_part[key]["Number"])
        # remove the last timepoint from model datacollector, because it doesn't exist in agent datacollcetion
        if Number[-1]==0:
            time_list=time_list[:-1]
        t_endpoints_a = [x for x in t_endpoints_time if x<=time_list[-1]+interval] # remove endpoints that haven't been simulated
        mins = [[abs(x-y) for y in time_list] for x in t_endpoints_a]
        min_index = [i for x in mins for i in range(len(x)) if min(x)==x[i]] # index of which timestep is closest to endpoint
        t_endpoints_a = [time_list[i] for i in min_index] # generate endpoints from timsteps closest
          
        db_m_part[key] = db_m_part[key][db_m_part[key]['Time'].isin(t_endpoints_a)]
        if fixed_params['data_agents']: # only if data for agents was saved    
            db_a_part[key] = db_a_part[key][db_a_part[key]['Time'].isin(t_endpoints_a)]
   
    # add iteration number to row
    for key in db_m_part:
        l = len(db_m_part[key])
        db_m_part[key]=pd.concat([db_m_part[key],
                                         pd.DataFrame({'iteration' :[key[-1]]*l
                                                       },index=db_m_part[key].index)]
                                 , axis=1)
    #%%
    if fixed_params['data_agents']: # only if data for agents was saved    
        #caclulate mean and std for some columns from agent dictionary and add to db_m_part
        # this is only interesting if variable in column are heterogeneous across population
        col_agents = [list(db_a_part[key].columns) for key in db_a_part]
        col_model = [list(db_m_part[key].columns) for key in db_m_part]
        column_names = [[y for y in x if (y not in ["Species","Time","AgentID","timestep"])] for x in col_agents]
        
        for key,column in zip(db_m_part,column_names):
            tp=list(dict.fromkeys([x[0] for x in db_a_part[key].index])) # timepoints in agent dataframe
            indexes=[]
            for t in tp: # go through all time points
                # make lists with indexes per species for each time point
                indexes_species = [[tuple([t,index]) for index,row in db_a_part[key].loc[t].iterrows() if row["Species"]==spec_name] for spec_name in species_names]
                indexes.append(indexes_species)
            indexes = [[x[i] for x in indexes] for i in range(len(indexes[0]))] # resort, switch 1st and 2nd axis
                
            for col in column:
                for s in range(len(species_names)):
                    mn = [np.nanmean([db_a_part[key].loc[x_t][col]]) for x_t in indexes[s]]
                    stds = [np.nanstd([db_a_part[key].loc[x_t][col]]) for x_t in indexes[s]]
                    # calculate mean and std for each column in list
                        # add to db_m_part
                    db_m_part[key]=pd.concat([db_m_part[key],
                                             pd.DataFrame({'m-'+col + '_'+ str(species_names[s]) :mn ,
                                                           'std-'+col + '_'+ str(species_names[s]) :stds
                                                           },index=db_m_part[key].index)],
                             axis=1, join='inner')
    else:
        column_names = [[] for i in range(len(db_m_part))]
    #%%
                               
    # Calculations on part of the dataframe
    # mean and std of ab distribution
    # this is only interesting if variable AB is hetergogeneous on grid
    
    for key in db_m_part:
        if 'AB_I' in db_m_part[key].columns:
            for y in range(len(ab_conc_range[0])): # for every drug in the simulation
                db_m_part[key] = pd.concat([db_m_part[key],
                           pd.DataFrame({'m-AB_'+"I"*(y+1): [np.mean(x) for x in db_m_part[key]['AB_'+"I"*(y+1)]],
                                         'std-AB_'+"I"*(y+1): [np.std(x) for x in db_m_part[key]['AB_'+"I"*(y+1)]]},
                                    index=db_m_part[key].index)],
                       axis=1, join='inner') 
                           
    for key in db_m_part: 
        if 'AB_I' in db_m_part[key].columns:
        # drop columns
            db_m_part[key]=db_m_part[key].drop(['AB_'+'I'*(y+1) for y in range(len(ab_conc_range[0]))],axis=1)

        
    # save processed data of only select timepoints
    file = open(path + "/db_m_part.pkl", "wb")
    pickle.dump(db_m_part, file)
    file.close()
    
    if fixed_params['data_agents']: # only if data for agents was saved    
        file = open(path + "/db_a_part.pkl", "wb")
        pickle.dump(db_a_part, file)
        file.close()   
    
    t1=time.time()
    print(t1-t0)   
    

    #%%
    # ------------------------------------------------------------------------------------------------------------------                              
    # ------------------------------------------------------------------------------------------------------------------                   
    # Dictionary of Dataframe with info about all iterations together, average values etc. per combination: dict_model
    # Calculate mean and std of species number between iterations at each time point
    t0 = time.time() # start timer   
    # make a copy only with select columns for which to fill up dbm and then calculate the population mean and stds
    column_s = [['iteration','Time','Number'] 
                +['#'+str(species_names[i]) for i in range(len(species_names))]
                +['Max_Av_k_net_'+str(species_names[i]) for i in range(len(species_names))]
                +['m-'+column_names[j][c_i]+'_'+str(species_names[i]) for i in range(len(species_names)) for c_i in range(len(column_names[j]))]
                for j in range(len(column_names))]
    db_cs = dict(zip(db_a_part.keys(),column_s))
    
    db_m = db_m_part.copy()
    for key in db_m:
        column_s_j = db_cs[key]
        db_m[key]=db_m[key][column_s_j]
    
    # make a dataframe with the number and time from all iterations
    all_keys = list(db_m.keys())
    keys = list(dict.fromkeys([x[:-1] for x in all_keys])) # same keys but without the iteration is the key list for the new dict
    
    bool_lengths = []
    bool_times = []
    dfs=[] # list of dataframes
    for key in keys:
        key_i=[key+(i,) for i in range(its)]
        # determine maximal length of simulation for one comb
        max_len = max([len(db_m[k_i]) for k_i in key_i])       
        last_time = [list(db_m[k_i]['Time'])[-1] for k_i in key_i]
           
        same_lengths = all([len(db_m[k_i])==max_len for k_i in key_i]) # this is set to false when not all iteration have the same length
        same_times = all([list(db_m[k_i]['Time'])[-1]==last_time[0] for k_i in key_i]) # this is set to false when not all last time points are the same
        bool_lengths.append(same_lengths)
        bool_times.append(same_times)
        
        if not same_lengths:
            for k_i in key_i:
                while len(db_m[k_i]) is not max_len:
                    db_m[k_i]=db_m[k_i].append(pd.DataFrame([[list(k_i)[-1],np.nan]+[0]*(len(species_names)+1)+[np.nan]*(len(species_names))*1+[np.nan]*len(column_names[0])*(len(species_names))],columns=db_cs[k_i]),ignore_index=True) 
       
        same_lengths = all([len(db_m[k_i])==max_len for k_i in key_i]) # recalculate
        
        if same_lengths:
            mean_time_vector = [int(np.nanmean([x[i] for x in [list(db_m[k_i]['Time']) for k_i in key_i]])) for i in range(max_len)]
            means = [[np.nanmean([x[i] for x in [list(db_m[k_i][col]) for k_i in key_i]]) for i in range(max_len)] for col in db_cs[key_i[0]][2::]]
            stds = [[np.nanstd([x[i] for x in [list(db_m[k_i][col]) for k_i in key_i]]) for i in range(max_len)] for col in db_cs[key_i[0]][2::]]
            x=[mean_time_vector]+means+stds
            x=[[y[i] for y in x ] for i in range(len(x[0]))]
            columns_new=[db_cs[key_i[0]][1]]+['m-'+ x for x in db_cs[key_i[0]][2::]]+['std-'+ x for x in db_cs[key_i[0]][2::]]

        dm_means=pd.DataFrame(x,columns=columns_new)
        dfs.append(dm_means)
        
    del [db_m]
    dict_model = dict(zip(keys,dfs))
    #%%
    # add AB_fix to dict_model
    for key in dict_model:
        for y in range(len(ab_conc_range[0])): # for every drug in the simulation
            df=pd.DataFrame(db_m_part[key+(0,)]['AB_fix_'+'I'*(y+1)])
            df.reset_index(drop=True, inplace=True)
            dict_model[key]=pd.concat([dict_model[key],df],axis=1)
               
    # save processed data with averages over all iterations
    file = open(path + "/dict_model.pkl", "wb")
    pickle.dump(dict_model, file)
    file.close()  
       
    t1=time.time()
    print(t1-t0) 
