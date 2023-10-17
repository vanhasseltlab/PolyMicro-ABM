# -*- coding: utf-8 -*-
"""
Agent based model of an interactive, multispecies community during antimicrobial treatment
Script to run set of simulations
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
"""
#%reset -f
#%%
import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+"/..")
from Model_v19 import MyModel
from PostSim_v16 import main
sys.path.append(dir_path)
from Functions import scale_intradius
import pickle
import mesa
import numpy as np
import pandas as pd
import time
import numpy
import copy

#%%
#-----------------------------------------------------------------------------------------------------------------------
# SET SPECIES CHARACTERISTICS
##-----------------------------------------------------------------------------------------------------------------------
maxgrowthrate = 0.015 # 1/min , approximately 66 min generation time
species_MIC = [6] # list of MIC values of species A to every drug in the simulation, in this example there is only one drug, so it's a list with one element
k_d0 = 0.2/60# natural death rate
species_A = ["A", species_MIC, maxgrowthrate, k_d0] # characteristics of one species in the format [species name, species_MIC, maximum growth rate, maximum kill rate]
species_B = ["B", species_MIC, maxgrowthrate, k_d0] # characteristics of one species in the format [species name, species_MIC, maximum growth rate, maximum kill rate]
all_species = [species_A,species_B] # all_species is a list of the characteristics of all species

# interactions
all_interactions = [
        [[False, False],
         ["MIC-", False]]
        ,
        [[False,False],
         [[1], False]]
        ,
        [[False,False], # third list in all_iteractions specifies the interaction radius for each interaction, 
         [False, False]]    # if int_radius=0: only cells on the own patch, int_radius=1: cell on the own patch and the 8 cells around, int_radius =2: like int_radius=1+ 1 row/column of cells around
        ]

int_radius_range = [1,0]

all_interactions_range = [copy.deepcopy(all_interactions) for r in range(len(int_radius_range))]
for r in range(len(int_radius_range)):
    all_interactions_range[r][2][1][0] = int_radius_range[r]

all_interactions_env = [[False,False],
                        [False,False],
                        [False,False]]
all_interactions_env=[all_interactions_env]

#drugs
k = [5,1.75,4,2] # hypothetical (original) & Tobramycin & Ceftazidime , hill parameters for each drug # Mouton& Vink Table 1 # independent of units
G_min = [-3/60, -0.0654, None, None] # for bactericidal antibiotics only
killtype = ['C','C','S','S'] # ADD static 'S' or for cidal 'C' to the drug list
drugtype=['prop','prop','prop','prop']

all_drugs= [[[G_min[i]],[k[i]],[killtype[i]],[drugtype[i]]] for i in range(len(k))]

ab_conc_range_CI = [[x] for x in [x*species_MIC[0] for x in [1,1.15,1.3,1.45,1.65,2.1,3,4,6]]]
ab_conc_range_CD = [[x] for x in [x*species_MIC[0] for x in [1,1.5,2,2.6,3.25,3.75,4.25,4.75,6]]]
ab_conc_range_SI = [[x] for x in [x*species_MIC[0] for x in [1,1.075,1.15,1.26,1.4,1.65,2.5,4,6]]]
ab_conc_range_SD = [[x] for x in [x*species_MIC[0] for x in [1,1.125,1.25,1.5,2,2.8,3.5,4.5,6]]]
ab_conc_range_alldrugs = [ab_conc_range_CI,ab_conc_range_CD,ab_conc_range_SI,ab_conc_range_SD]
# ab_conc_range_alldrugs = [ab_conc_range_CI] #just for test
#interaction parameters
# DEFAULT PARAMETERS FOR INT_RADIUS=1
# MIC-
MIC_change_range_original = [-0.01,-0.025,-0.05] #y for all drugs

#zaro is the natural lower limit for MIC, k_repmax and ab conc
# y_range = [-0.9999, np.nan] # range for MIC changes, TRY TO AVOID UPPER LIMIT
                    
#-----------------------------------------------------------------------------------------------------------------------
# SET SIMULATION PARAMETERS
#-----------------------------------------------------------------------------------------------------------------------
runtime = 1/6*24*60 # runtime in mins
# runtime = 2/6*60 #just for test
its=50 # number of iterations
# its=2 #just for test
data_collection = 10 # interval for datacollection in mins
# data_agents = False # set to False for test simulations, requires less storage
data_agents = True
capacity = 10**3# maximal number of species in simulation at the same time
N_c = 6.8 # number of cells per patch at capacity, with assumed realistic capacity=10**9 CFU/ml
N_p = capacity/N_c
gridrange = round(numpy.sqrt(N_p)) # this is 12 for a capacity of 10**3
gridrange_range=[gridrange]
# n_steps_move_type=[[0,0.1,0.25,0.5,0.75,1],[1]]
n_steps_move_type=[[0],[1]] #just for test
n_init = round(1* N_c*gridrange*gridrange/2) # number of bacteria at the start of the simulation, 100% of total capacity
stepsize = [10] # length of one timestep in simulation, in mins
max_steps = round(runtime/stepsize[0]) # maximal number of steps to be performed
# move_type=['stepwise','random']
move_type=['stepwise'] #just for test

#%%
# ------------------------------------------------------------------------------------------------------------------               
#-----------------------------------------------------------------------------------------------------------------------
# BATCH RUN SIMULATIONs
#-----------------------------------------------------------------------------------------------------------------------    
# ------------------------------------------------------------------------------------------------------------------               

time_start=time.time()
db_agents={}
db_model={}

# run simulations in batch 
# for one species at different drug concentrations
for drugs,ab_conc_range in zip(all_drugs,ab_conc_range_alldrugs):
    ab_conc_time_range = [[list(np.repeat(x, max_steps+1)) for x in ab_conc] for ab_conc in ab_conc_range] # create vectors with concentration over time at every timestep for each defined concentration
    drugs=[drugs]
    for m_type,n_steps_r in zip(move_type,n_steps_move_type):
        for n_steps in n_steps_r:
            for all_interactions in all_interactions_range:
                int_radius = all_interactions[2][1][0]
                all_interactions=[all_interactions]
                MIC_change_range = [scale_intradius(v,1,int_radius) for v in MIC_change_range_original]
                for MIC_change in MIC_change_range:
                    for ab_conc,ab_conc_time in zip(ab_conc_range,ab_conc_time_range):
                        all_species = [[["A", [6], 0.015, 0.2/60],["B", ab_conc, 0.015, 0.2/60]]]
                        ab_conc_time = [ab_conc_time]
                        params = {"capacity": capacity,
                                  "all_species": all_species,
                                  "all_interactions": all_interactions,
                                  "all_interactions_env": all_interactions_env,
                                  "ab_conc_time": ab_conc_time,
                                  "drugs": drugs,
                                  "runtime": runtime,
                                  "gridrange": gridrange_range[0],
                                  "stepsize": stepsize[0],
                                  "data_collection": data_collection,
                                  "MIC_change": MIC_change,
                                  "data_agents": data_agents,
                                  "n_steps": n_steps,
                                  "move_type": m_type,
                                  "n": range(n_init,n_init+1)}
                        
                        fixed_params = dict(list(params.items())[::-1])
                        variable_params = {"n": range(n_init,n_init+1)}
    		            
                        results = mesa.batch_run(MyModel,parameters=params,iterations=its,max_steps=round(runtime/stepsize[0])-1,number_processes=None,data_collection_period=int(data_collection/stepsize[0]),display_progress=True)
                        
                        results_df = pd.DataFrame(results)
                        
                        #recreate data_batch_agents and data_batch_model dataframes from results_df to be able to use old PostSim file Postsim_v16 and older
                        
                        for i in range(its):
                            data_model = results_df.copy()
                            #select iteration
                            data_model=data_model[data_model['iteration']==i]
                            data_model=data_model.drop(list(data_model.index)[0]) #drop first row because it is duplicate
                            data_model=data_model.reset_index(drop=True)
                            #select columns
                            data_model=data_model[['timestep', 'Time', 'Number', 'Number_Species', 'AB_fix', 'Int_current','Int_failed', 'Int_give', 'Int_receive', 'affected', 'Max_Av_k_net']]
                            db_model[(results_df['drugs'].loc[0][1][0],
                                results_df['n_steps'].loc[0],
                                results_df['all_interactions'].loc[0][0][1][0],
                                results_df['all_interactions'].loc[0][2][1][0],
                                round(results_df['MIC_change'].loc[0],ndigits=6),
                                results_df['move_type'].loc[0],
                                round(results_df['ab_conc_time'].loc[0][0][0],ndigits=2),
                                i)] = data_model.copy()
                            if data_agents==False:
                                db_agents[(results_df['drugs'].loc[0][1][0],
                                results_df['n_steps'].loc[0],
                                results_df['all_interactions'].loc[0][0][1][0],
                                results_df['all_interactions'].loc[0][2][1][0],
                                round(results_df['MIC_change'].loc[0],ndigits=6),
                                results_df['move_type'].loc[0],
                                round(results_df['ab_conc_time'].loc[0][0][0],ndigits=2),
                                i)] = pd.DataFrame()
                            elif data_agents==True:
                                data_batch_agents = results_df.copy()
                                #select iteration
                                data_batch_agents=data_batch_agents[data_batch_agents['iteration']==i]
                                data_batch_agents=data_batch_agents.set_index(['Step','AgentID'])
                                #select columns
                                data_batch_agents=data_batch_agents[['MIC', 'Max_replicationrate', 'Species', 'Growthrate', 'Killrate', 'NetGrowthrate']]
                                db_agents[(results_df['drugs'].loc[0][1][0],
                                results_df['n_steps'].loc[0],
                                results_df['all_interactions'].loc[0][0][1][0],
                                results_df['all_interactions'].loc[0][2][1][0],
                                round(results_df['MIC_change'].loc[0],ndigits=6),
                                results_df['move_type'].loc[0],
                                round(results_df['ab_conc_time'].loc[0][0][0],ndigits=2),
                                i)] = data_batch_agents.copy()

time_end=time.time()
print(time_end-time_start)
##%% save data to files  

folder_name="Simulations_v16_PW_Const_MIC-VI"
path = dir_path + "/../../Data/" + folder_name
os.mkdir(path)

##%%
file = open(path + "/species.pkl", "wb")
pickle.dump([all_species[0],all_interactions_range[0],all_interactions_env[0],ab_conc_range], file)
file.close()  

file = open(path + "/settings.pkl", "wb")
pickle.dump([capacity,n_init,its,runtime,gridrange,move_type,stepsize,data_collection,variable_params,fixed_params], file)
file.close()  

file = open(path + "/db_agents.pkl", "wb")
pickle.dump(db_agents, file)
file.close()    

file = open(path + "/db_model.pkl", "wb")
pickle.dump(db_model, file)
file.close() 

#%%
data_interval = 30
main(data_interval,path,db_agents,db_model,all_species[0],all_interactions_range,all_interactions_env[0],ab_conc_range,capacity,n_init,its,runtime,gridrange,move_type,stepsize,data_collection,variable_params,fixed_params)