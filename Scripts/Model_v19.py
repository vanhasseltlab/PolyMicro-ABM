# -*- coding: utf-8 -*-
"""
An Agent based modeling framework of an interactive, polymicrobial community during antimicrobial treatment
Main model script
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
"""
# load packages
from mesa import Agent, Model
from mesa.time import RandomActivation
# from mesa.time import StagedActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
import numpy as np
import copy

# Agent object class, each agent represents a bacterium
class MyMicrobe(Agent):
    # initialization of a new agent
    def __init__(self, model, name, species, interactions):
        super().__init__(name, model)
        self.alive = True
        # Species specific parameters, these parameters do not change during the simulation
        self.species = species # variable containing most important species characteristics
        self.species_name = self.species[0] # name of species
        self.MIC_original = self.species[1] # original species-specific MIC of each drug present in simulation in list format: [MIC of drug I, MIC of drug II]
        self.k_repmax_original = self.species[2] # original maximal replication rate for this species
        self.k_d0 = self.species[3] # maximal replication rate for this species
        self.interactions = interactions # agent-agent interactions of this species; in list format [interaction partner&type, interaction target, interaction radius]
        # Parameters specific to each agent, these parameters do not change during the simulation
        self.name = name # unique agent identifier, an integer assigned when agent is created, reused by first daughter, automatically collected by datacollector
        # Parameters specific to each agent, these parameters can change during the simulation as a result of interactions
        self.MIC = copy.deepcopy(self.MIC_original) # current MIC of this agent, initialized to be species specific original MIC  
        self.k_repmax = self.k_repmax_original # current k_repmax for this agent, initialized to be species specific original replication rate 
                
    # STEP FUNCTION FOR AGENTS
    def step(self):
        self.count() # running counts for data collection
        self.divide_die()
        self.interact()
        self.move_agent()
        self.update_interaction_effects()

    # FUNCTIONS FOR ACTIONS DIVIDE OR DIE
    
    # evaluates net growth probability to decide if the agent dies or divides
    def divide_die(self):
        if self.alive:
        # UPDATE K_NET
            self.set_netgrowth([x[self.pos[0]][self.pos[1]] for x in self.model.ab_conc_grid_alldrugs])           
        #  decide whether DIVIDE OR DIE based on if K_NET is positive or negative
            if self.K_NET < 0:
                if self.model.random.random() < abs(self.model.stepsize*self.K_NET):
                    self.die()
            elif self.K_NET > 0:
                if self.model.random.random() < (self.model.stepsize*self.K_NET):
                    self.divide()
            self.model.k_nets[[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]].append(self.K_NET) #collects k_net of every agent at one timestep
    
    # calculates net growth of the agent at current timestep taking interaction effects into account
    def set_netgrowth(self, conc): #calculates net growth for simulations with one drug only
        k_rep = self.set_k_rep(self.k_repmax)
        self.set_drugeffect(conc)
        if self.model.drugtype[0]=='prop':
            self.replicationrate = k_rep* self.E_S # actual replication rate including drug and interaction effect at this time point
            self.deathrate = self.k_d0* self.E_C # actual death rate including drug and interaction effect at this time point
        elif self.model.drugtype[0]=='add':
            self.replicationrate = k_rep-self.E_S # actual replication rate including drug and interaction effect at this time point
            self.deathrate = self.k_d0+self.E_C # actual death rate including drug and interaction effect at this time point
        self.K_NET = self.replicationrate - self.deathrate # overall net growth rate including drug effect and interaction effects
    
    # applies capacity limitation to replication rate
    # whether or not interaction effects are taking into account depends on the k_repmax input
    def set_k_rep(self, k_repmax):
        k_d0 = self.k_d0
        B = self.model.number_microbes
        K = self.model.capacity
        k_rep = k_repmax *(1- B/K*(1- k_d0/k_repmax))
        return k_rep

    # calculates drug effect from the local drug concentration conc taking into account interaction effects; for simulations with one drug only
    def set_drugeffect(self, conc):
        k_rep = self.set_k_rep(self.k_repmax_original) # replication rate disregarding interaction effects or drug effects
        all_E_S=[]
        all_E_C=[]
        
        for drugnumber in range(len(conc)):
            if self.model.killtype[drugnumber] == 'C' and self.model.drugtype[drugnumber]=='prop': # for bactericidal drugs
                k=self.model.k[drugnumber]
                G_min = self.model.G_min[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                E_C = 1 + (((k_rep-G_min)/k_d0 - 1)*(c/MIC)**k)/((c/MIC)**k + G_min/(k_d0-k_rep))
                E_S = 1 # no effect
                
            elif self.model.killtype[drugnumber] == 'S'and self.model.drugtype[drugnumber]=='prop':  # for bacteriostatic drugs
                k=self.model.k[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                E_S = 1 - ((c/MIC)**k)/(((c/MIC)**k) + k_d0/(k_rep-k_d0))
                E_C = 1 # no effect
                
            elif self.model.killtype[drugnumber] == 'C' and self.model.drugtype[drugnumber]=='add': # for bactericidal drugs
                k=self.model.k[drugnumber]
                G_min = self.model.G_min[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                E_C = (k_rep-k_d0-G_min)*(((c/MIC)**k)/((c/MIC)**k - (G_min/(k_rep-k_d0))))
                E_S = 0 # no effect
                
            elif self.model.killtype[drugnumber] == 'S' and self.model.drugtype[drugnumber]=='add':  # for bacteriostatic drugs
                k=self.model.k[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                E_S = k_rep* (((c/MIC)**k)/(((c/MIC)**k) + k_d0/(k_rep-k_d0)))
                E_C = 0 # no effect
                
            all_E_S.append(E_S)
            all_E_C.append(E_C)
        
        self.E_S = all_E_S[0]
        self.E_C = all_E_C[0]

    # "kills" the agent by setting the alive parameter of the agent to False and removing it from the counters, does not remove the agent object yet
    def die(self):
        self.model.kill_microbes_list.append(self) # add to kill list
        self.model.number_microbes =  self.model.number_microbes - 1
        self.model.number_microbes_grid[self.pos[0]][self.pos[1]] -=1
        self.model.number_per_species[[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]] -= 1
        self.model.number_per_species_grid[self.pos[0]][self.pos[1]][[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]] -=1
        self.alive = False
    
    # duplicates the agent, the first daughter is created by resetting the mother agent, for the second data a new agent object is created which copies the characteristics of the first one        
    def divide(self):
        # reuse the set agent characteristics back to imitate newborn agent (with or) without inheriting charactersistics, 1st daughter
        self.reset_microbe() # reset and move to random patch, reusing the already existing agent saves computation time
        # create one more new agent of the same species as copy of 1st daughter on model level
        self.model.create_copy(self)
        self.model.ab_conc # WHY IS THIS HERE?
    
    # resets characteristics of mother agent to create first daughter, unique agent identifier "name" stays the same    
    def reset_microbe(self):
        #reset interaction effects
        self.MIC = copy.deepcopy(self.MIC_original)# current MIC of this agent, initialized to be species specific MIC
        self.k_repmax = self.k_repmax_original
        if self.model.move_type == "random":
            self.move_rand()

    # FUNCTIONS FOR ACTION INTERACT
    # initiates interactions if the agent is alive and has interaction abilities
    def interact(self):
        if self.alive:
            # INTERACT with agents
            # interactions with other agents in interaction radius
            #CHANGE TO INTERACTION RADIUS AND NOT ONLY THIS PATCH 
            if any(self.interactions[0]): #and len(self.model.grid.get_cell_list_contents([self.pos]))>1: #if this agent interacts with any other species and there is at least one more agent in the same patch
                self.interact_agents()
    
    # executes interactions of 4 different types: positive resistance interaction (MIC+), negative resistance interaction (MIC-), positive growth interaction (G+) and negative growth interaction (G-)
    def interact_agents(self): # interact with agents from the target species in the interaction radius
        inter_types = [x for x in self.interactions[0] if x is not False] # types of interactions that this agent can initiate
        inter_species = [self.model.species_types[i] for i in range(len(self.interactions[0])) if self.interactions[0][i] is not False] # target species affected by any interaction
        inter_drug = [self.interactions[1][i] for i in range(len(self.interactions[0])) if self.interactions[0][i] is not False] # drug for which the resistance interaction effects are executed, irrelevant for simulations with one drug 
        inter_radius = [self.interactions[2][i] for i in range(len(self.interactions[0])) if self.interactions[0][i] is not False] # interaction radius for each outgoing int edge
        cellmates_done = []
        for r in range(max(inter_radius)+1): # starting with the patch on which the agent resides, go through layers of patches around the active agent (self) until the most outer layer within the interaction radius is reached
            if r==0:
                cellmates = self.model.grid.get_cell_list_contents([self.pos]) # all cells in the same patch
            else:
                cellmates = self.model.grid.get_neighbors(self.pos, moore=True,include_center=False,radius=r)# True:Moore neighborhood (incl diagonals), do not include center
                if r>1:
                    cellmates = [x for x in cellmates if x not in cellmates_done]# remove cellmates from inner layers, only necessary for r>1, if r=1 center is not included anyhow
            for other_agent in cellmates: # go through all cellmates
                if self.random.random()<self.model.interaction_prob_agents: # interact with interaction probability, default 100%
                    for i in range(len(inter_species)): # for every possible interaction partner
                        if r<= inter_radius[i]:# check if the interaction that would happen with other_agent based on its species is in the int radius based on the type of interaction
                            if other_agent.species_name == inter_species[i]: # check if its a species to interacts with
                                
                            # positive resistance interaction    
                                for drug_number in range(len(self.model.ab_conc)):
                                    if inter_types[i] == "MIC+": # check which type of interaction
                                        if inter_drug[i][drug_number]>0: # resistance interactions only happen when a drug is present
                                            if other_agent.MIC[drug_number] < (self.model.y_max+1)*other_agent.MIC_original[drug_number]: #interact if the agent has not reached the maximal MIC
                                                # APPLY INTERACTION EFFECT TO MIC
                                                other_agent.MIC[drug_number] = min(other_agent.MIC[drug_number] + self.model.MIC_change*other_agent.MIC_original[drug_number]*inter_drug[i][drug_number], (self.model.y_max+1)*other_agent.MIC_original[drug_number])
                                                self.model.int_current[0]+=1
                                                self.model.int_give[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]]+=1
                                                self.model.int_receive[[i for i in range(len(self.model.species_types)) if other_agent.species_name==self.model.species_types[i]][0]]+=1
                                            else:
                                                self.model.int_failed[0]+=1
                                                
                                    # negative resistance interaction            
                                    elif inter_types[i] == "MIC-":
                                        if inter_drug[i][drug_number]>0:
                                            if other_agent.MIC[drug_number] > (self.model.y_min+1)*other_agent.MIC_original[drug_number]: #interact if the agent has not reached the minimal MIC
                                                # APPLY INTERACTION EFFECT TO MIC
                                                other_agent.MIC[drug_number] = max(other_agent.MIC[drug_number] + self.model.MIC_change*other_agent.MIC_original[drug_number]*inter_drug[i][drug_number], (self.model.y_min+1)*other_agent.MIC_original[drug_number])
                                                self.model.int_current[1]+=1
                                                self.model.int_give[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]]+=1
                                                self.model.int_receive[[i for i in range(len(self.model.species_types)) if other_agent.species_name==self.model.species_types[i]][0]]+=1        
                                            else:
                                                self.model.int_failed[1]+=1                                  
                                    
                                # negative growthrate interaction 
                                if inter_types[i] == "G-":
                                    if (other_agent.k_repmax) > (self.model.z_min+1)*other_agent.k_repmax_original: #interact only if minimal k_repmax has not been reached yet
                                        # APPLY INTERACTION EFFECT TO K_REPMAX
                                        other_agent.k_repmax = max((other_agent.k_repmax+self.model.G_change*other_agent.k_repmax_original),(self.model.z_min+1)*other_agent.k_repmax_original)
                                        self.model.int_current[2]+=1
                                        self.model.int_give[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]]+=1
                                        self.model.int_receive[[i for i in range(len(self.model.species_types)) if other_agent.species_name==self.model.species_types[i]][0]]+=1
                                    else: 
                                        self.model.int_failed[2]+=1
                                        
                                # positive growthrate interaction         
                                elif inter_types[i] == "G+":
                                    if (other_agent.k_repmax) < (self.model.z_max+1)*other_agent.k_repmax_original: #interact only if maximal k_repmax has not been reached yet
                                        # APPLY INTERACTION EFFECT TO K_REPMAX
                                        other_agent.k_repmax = min((other_agent.k_repmax+self.model.G_change*other_agent.k_repmax_original),(self.model.z_max+1)*other_agent.k_repmax_original)
                                        self.model.int_current[3]+=1
                                        self.model.int_give[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]]+=1
                                        self.model.int_receive[[i for i in range(len(self.model.species_types)) if other_agent.species_name==self.model.species_types[i]][0]]+=1
                                    else:
                                        self.model.int_failed[3]+=1
                                    
            if r>0:
                cellmates_done = cellmates_done + cellmates # adds all cells that have been looked at in radius of r into this list, only save cells not in center, those are automatcally exluded in get_neighbours with include_center=False
        
    # count if this agent has been affected by any interactions at time of activation
    def count(self):
        # order of ints in affected: ["MIC+","MIC-","G-","G+"]
        # assuming each species is only affected by one interaction type at a time and assuming one drug present at a time
        if self.MIC[0]>self.MIC_original[0]:
            self.model.affected[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]][0]+=1
        if self.MIC[0]<self.MIC_original[0]:
            self.model.affected[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]][1]+=1
        if self.k_repmax < self.k_repmax_original:
            self.model.affected[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]][2]+=1
        if self.k_repmax > self.k_repmax_original:
            self.model.affected[[i for i in range(len(self.model.species_types)) if self.species_name==self.model.species_types[i]][0]][3]+=1


    # FUNCTIONS FOR ACTION MOVE
    
    # moves agent n_steps according to random or stepwise settings
    def move_agent(self):
        if self.alive:
            self.model.number_microbes_grid[self.pos[0]][self.pos[1]] -=1 # remove counter from old position
            self.model.number_per_species_grid[self.pos[0]][self.pos[1]][[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]] -=1
            
            p=self.model.n_steps%1
            
            if self.model.move_type == "stepwise":
                if self.random.random()<p: # if p is 0.3 and n_steps is 1.3, then in 33 percent of cases the agent takes 2 steps and in 66% of cases 1 step
                    for n in range(int(np.ceil(self.model.n_steps))): # n_steps for each timestep
                        # 1) MOVE on grid
                        self.move()
                else:
                    for n in range(int(np.floor(self.model.n_steps))): # n_steps for each timestep
                        # 1) MOVE on grid
                        self.move()
                        
            elif self.model.move_type == "random": # random always moves every agent to a random position
                self.move_rand()

            self.model.number_microbes_grid[self.pos[0]][self.pos[1]] +=1 # add counter at new position
            self.model.number_per_species_grid[self.pos[0]][self.pos[1]][[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]] +=1
    
    # move agent to one of 8 neighbouring patches or stay on the same patch, all with equal chances
    def move(self):
        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            moore=True,
            include_center=True)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)
            
        
    # UPDATE INTERACTION EFFECTS IN PREPARATION FOR NEXT TIMESTEP
    # at the end of an activation of each agent the interaction effects are set back to original
    def update_interaction_effects(self):
        # this means: interactions that are happending before the agent is activated take effect in the current timestep and interactions that happen after an agent is activated count for the next timestep
        if self.alive:
            for drug_number in range(len(self.model.ab_conc)):
                self.MIC[drug_number] = copy.deepcopy(self.MIC_original[drug_number])
            self.k_repmax = self.k_repmax_original
        

    # FUNCTIONS (MAINLY) CALLED FROM MODEL LEVEL       
    # move agent to position pos        
    def move_daughter(self, pos):
        self.model.number_microbes_grid[self.pos[0]][self.pos[1]] -=1 # remove counter from old position
        self.model.grid.move_agent(self,(pos[0],pos[1]))
        self.model.number_microbes_grid[self.pos[0]][self.pos[1]] +=1 # add counter at new position
    
    # move agent to a random position on the grid
    def move_rand(self):
        coords = (self.model.random.randrange(0, self.model.gridrange), self.model.random.randrange(0, self.model.gridrange))
        self.model.grid.move_agent(self, coords)

# main model object class, for each iteration of the simulation, one model object is created which then creates agent objects        
class MyModel(Model):
    # Model initialization
    def __init__(self, n, capacity, all_species, all_interactions, all_interactions_env, ab_conc_time, drugs, runtime, gridrange, move_type="stepwise",
                 n_steps=60, data_collection = 10, data_agents = False, stepsize = 1,MIC_change = 0.1, G_change = 0.1, y_range = [-0.9999, 10], 
                 z_range=[-0.9999, 2],interaction_prob_agents  = 1):
        # varaiable all_interactions_env unused: eliminate from input scripts
        super().__init__()
        self.running = True # necessary for batch run, is set to False when all bacteria have died, This stops the simulation run
        
        # simulation parameters
        self.runtime = runtime
        self.stepsize = stepsize # in mins
        self.starting_number = n
        self.time = 0
        self.step_id = 0
        self.data_collection = data_collection
        self.data_agents = data_agents
        self.schedule = RandomActivation(self)
        
        # model parameters
        self.capacity = capacity # maximal number of agents
        self.y_min = y_range[0] # one value from which the minimum MIC can be calcuated: MIC_min = (y_min+1)*MIC_original for each species
        self.y_max = y_range[1] 
        self.z_min = z_range[0] # one value 
        self.z_max = z_range[1] # one value
         
        # growth dynamics and antibiotic treatment parameters, simulation specific
        self.MIC_change = MIC_change # unit change of MIC for one interaction # in fraction of original value
        self.G_change = G_change # unit change of k_repmax forone interaction # in fraction of original value
        self.interaction_prob_agents = interaction_prob_agents
      
        # spatial grid and movement settings
        self.gridrange = gridrange
        self.grid = MultiGrid(self.gridrange, self.gridrange, torus=True)
        self.move_type = move_type
        self.n_steps = n_steps
        
        #create grid for initial distribution of antibiotic concentration
        self.ab_conc_time = ab_conc_time # applied oncentration over time as predefined, not on grid yet
        self.set_ab_conc() # apply concentration to each patch on grid
        
        # drug parameter and characteristics
        self.G_min=drugs[0]
        self.k = drugs[1]
        self.killtype = drugs[2]
        self.drugtype = drugs[3]
        
        
        # counter initialization
        self.last_name = -1 # to assign unique name, start at 0 for first agent.name, so initialize last_name to -1
        self.number_microbes = 0 # counts alive agents
        self.number_microbes_grid = [[0 for j in range(self.gridrange)] for i in range(self.gridrange)]
        self.number_species = len(all_species)
        self.species_types=[x[0] for x in all_species] # make a list of all the species
        self.number_per_species = [0 for x in range(self.number_species)]
        self.number_per_species_grid = [[[0 for x in range(self.number_species)] for j in range(self.gridrange)] for i in range(self.gridrange)]
        self.nm=copy.deepcopy(self.number_per_species)
        self.nmg=copy.deepcopy(self.number_per_species_grid)
        self.int_current=[] # creates variable that stores interactions that happened in the time interval between two data collections
        self.int_give=[]
        self.int_receive=[]
        self.affected = [] # creates variable that collects number of bacteria of each species affected by each type of interaction # what if they are G+ and G- interactions at the same time and they level out
        self.int_failed=[]
        self.k_nets = [[np.nan] for i in range(self.number_species)] #list that collects k_net of every agent at every time step in datacollection interval, from that an average is calculated # one list for each species
        self.k_net_mean = [[np.nan] for i in range(self.number_species)] #counter that collects average (in population) k_net of every timestep
        self.k_net_max = [np.nan for i in range(self.number_species)]
        self.kill_microbes_list = list() # initialize kill list
        
        # create a random starting ditribution of species with approximately equal numbers of each species
        x = self.starting_number/self.number_species
        ns = [self.random.gauss(x,x*0.1) for i in range(self.number_species)]
        ns = [ i/sum(ns) for i in ns ] #divide each element by sum
        ns = [int(i*self.starting_number) for i in ns] # multiply each element by desired sum, this leads to approximately n_init bacteria at the start
        ns = [x+1 if x==0 else x for x in ns]
            
        # create agents
        for i in range(self.number_species):
            self.create_microbes(ns[i], all_species[i], [x[i] for x in all_interactions])# create chosen number of each type of species
        
        # set up data collection during runtime 
        if self.data_agents: #only collect agent reporters if data_agents is True
            self.datacollector = DataCollector(
                model_reporters={"timestep": "step_id", "Time": "time", "Number": "number_microbes",
                                 "Number_Species": "nm", 
                                 "AB_fix":"ab_conc", 
                                  "Int_current": lambda m: copy.deepcopy(m.int_current) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_failed": lambda m: copy.deepcopy(m.int_failed) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_give": lambda m: copy.deepcopy(m.int_give) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_receive": lambda m: copy.deepcopy(m.int_receive) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "affected": lambda m: copy.deepcopy(m.affected) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                 "Max_Av_k_net" : lambda m: copy.deepcopy(m.k_net_max)
                                 }, 

                agent_reporters={"MIC": lambda a: copy.deepcopy(a.MIC) if (any([x=='MIC-' for y in all_interactions[0] for x in y]) or any([x=='MIC+' for y in all_interactions[0] for x in y])) else None, 
                                 "Max_replicationrate":lambda a: a.k_repmax if (any([x=='G-' for y in all_interactions[0] for x in y]) or any([x=='G+' for y in all_interactions[0] for x in y])) else None,
                                 "Species":"species_name", 
                                 "Growthrate": "replicationrate",
                                 "Killrate": "deathrate",
                                 "NetGrowthrate": "K_NET"}
                )
        else: 
            self.datacollector = DataCollector(
                model_reporters={"timestep": "step_id", "Time": "time", "Number": "number_microbes",
                                 "Number_Species": "nm", 
                                 "AB_fix":"ab_conc", 
                                  "Int_current": lambda m: copy.deepcopy(m.int_current) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_failed": lambda m: copy.deepcopy(m.int_failed) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_give": lambda m: copy.deepcopy(m.int_give) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "Int_receive": lambda m: copy.deepcopy(m.int_receive) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                  "affected": lambda m: copy.deepcopy(m.affected) if any([x!=False for y in all_interactions[0] for x in y]) or any([x!=False for x in all_interactions[0]]) else None,
                                 "Max_Av_k_net" : lambda m: copy.deepcopy(m.k_net_max)
                                 })
            
        self.collect_data() # collect initial data before simulation starts
        
   
    def step(self):

        # UPDATE time and conc for next step
        self.step_id += 1
        self.time += self.stepsize
        self.set_ab_conc() # set antibiotic concentration to next concentration
        
        #reset counters or interactions
        self.int_current = [0]*5 # counts number of each type of interaction happening in this TIME INTERVAL with this order ["MIC+","MIC-","G-","G+","AB-"]
        self.int_give=[0]*self.number_species
        self.int_receive=[0]*self.number_species
        self.affected = [[0]*4]*self.number_species # counts number of affected bacteria from each type of interspecies interaction and for each species
        self.int_failed = [0]*5 #interactions that did not happen because the maximum change was reached
        self.k_nets = [[np.nan] for i in range(self.number_species)] #reset to np.nan for each timestep
        
        # 1) STEP AGENTS
        self.schedule.step() # activate agents in random order according to schedule and then run the agent's step function
        
        # 2) EXECUTE KILLING
        self.kill_microbes() # remove agent object of agents that have been placed on the kill list
        # this has to be done on the model level and cannot be done on the agent level
        
        # 4) STOP?
        # stop simulation run when all bacteria have died
        if self.number_microbes == 0:
            self.running = False
            
        # 3) COLLECT DATA 
        #calculate average k_net of this timestep to save for calculation of exact max k_net
        for i in range(self.number_species):
            self.k_net_mean[i].append(np.nanmean(self.k_nets[i]))
        
        # save variables for each time step for data collection
        # collect data approximately every x minutes as defined in data_collection
        # collect also for the last time step, both when stopped because no cells are left or when max_steps is reached
        if self.time%self.data_collection < self.stepsize or self.running==False or round(self.runtime/self.stepsize)==self.step_id: 
            self.k_net_max = [max(np.nanmax(self.k_net_mean[i]),np.nanmin(self.k_net_mean[i]),key=abs) for i in range(self.number_species)] #calculate the max of average k_nets every datacollection interval to save in datacollector
            self.k_net_mean =[[np.nan] for i in range(self.number_species)] #reset k_net_mean
            self.collect_data()


    # model functions in order of appearance in step function
    
    # set antibiotic concentration at every patch
    def set_ab_conc(self):
        self.ab_conc = [x[self.step_id] for x in self.ab_conc_time]
        self.ab_conc_grid_alldrugs=[[[self.ab_conc[drug_number] for j in range(self.gridrange)] for i in range(self.gridrange)] for drug_number in range(len(self.ab_conc))]
        
    # creates new agent objects
    def create_microbes(self, n, species, interactions): # use in initialization and during runtime through create_copy to set agent characteristics correctly
        for i in range(n):
            a = MyMicrobe(self, self.last_name + 1, species, interactions)# create new agent with unique name
            self.schedule.add(a) # schedule agent to be avtivated
            coords = (self.random.randrange(0, self.gridrange), self.random.randrange(0, self.gridrange))
            self.grid.place_agent(a, coords) # place agent on grid
            self.last_name += 1 # assign unique identifier
            # count agent in running counts
            self.number_microbes += 1
            self.number_microbes_grid[coords[0]][coords[1]] +=1
            self.number_per_species[[i for i in range(self.number_species) if species[0]==self.species_types[i]][0]] += 1
            self.number_per_species_grid[coords[0]][coords[1]][[i for i in range(self.number_species) if species[0]==self.species_types[i]][0]] +=1
        return a
  
    def create_copy(self, a): # creates the 2nd daugther as a copy of 1st daughter
        a_copy = self.create_microbes(1, a.species, a.interactions) # position is random
        # Species specific parameters, these parameters do not change during the simulation
        a_copy.species = a.species # variable containing most important species characteristics
        a_copy.species_name = a.species_name # name of species
        a_copy.MIC_original = copy.deepcopy(a.MIC_original) # original MIC for this species
        a_copy.k_repmax_original = a.k_repmax_original # maximal growth rate for this species
        a_copy.interactions = a.interactions # agent-agent interactions of this species
        # Parameters specific to each agent
        # MIC and k_repmax are the same as the 1st daughters
        a_copy.MIC = copy.deepcopy(a.MIC)
        a_copy.k_repmax = a_copy.k_repmax_original
        
       # move newly created agent into position according to settings
        if self.move_type != "random":
            a_copy.move_daughter(a.pos) # place at the same position as first daughter
        elif self.move_type == "random":
            a_copy.move_rand() # place randonly

    # remove agent object of agents that have been placed on the kill list
    def kill_microbes(self):
        for x in self.kill_microbes_list:
            self.grid.remove_agent(x)
            self.schedule.remove(x)
        self.kill_microbes_list = list() # reset list for next step

    # collects data
    def collect_data(self):
        self.nmg=copy.deepcopy(self.number_per_species_grid)
        self.nm=copy.deepcopy(self.number_per_species)
        self.datacollector.collect(self)
