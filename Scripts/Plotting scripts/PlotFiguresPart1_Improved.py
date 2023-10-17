#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:27:12 2022
Plotting script
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
"""

#%%
import os
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# working directory: '/home/cathi/Simulations/Software Supplementary'
# dir_path = os.path.dirname(os.path.abspath("__file__"))
# os.chdir(dir_path+'/../..')

#%% plot settings
plt.style.use('default')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.dpi']= 600
mpl.rcParams['lines.linewidth']  = 0.5
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.title_fontsize'] = 8
mpl.rcParams['axes.labelsize']=8
mpl.rcParams['axes.titlesize']=8
mpl.rcParams['font.size']=8
mpl.rcParams['xtick.labelsize']=6 
mpl.rcParams['ytick.labelsize']=6
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.size'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 0.5
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 0.5
mpl.rcParams['ytick.minor.width'] = 0.3
mpl.rcParams['grid.linewidth'] = 0.3
mpl.rcParams['lines.markersize'] = 0.5
mpl.rcParams['lines.marker'] = 'o'

#colors
colors_range=[['deepskyblue','mediumslateblue','deeppink'],['limegreen','gold','orangered']]
base_color='blue'

#%% 
# load combined df_all in pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_v16&v17.pkl", "rb")
df_all = pickle.load(file)
file.close() 

#%% plot CE curves - Figure "Figure6_highdpi_v2"
mpl.rcParams['figure.figsize']=[6.52,7.2]

m_type='stepwise'
int_rad=1
n_step=1
for int_rad in [0]:
    for n_step in [1]:            
        #take subset to plot
        df_plot = df_all.copy()
        df_plot=df_plot[df_plot['Cdep']=='ci']
        df_plot = df_plot[df_plot['move_type']==m_type] 
        df_plot = df_plot[df_plot['n_steps']==n_step]
        df_plot=df_plot.reset_index(drop=True)
        
        drugs_key=list(dict.fromkeys(df_plot['drug']))
        drugs_type_key=list(dict.fromkeys(df_plot['drug_type']))
        kill_type_key=list(dict.fromkeys(df_plot['kill_type']))
        cd_type_key=list(dict.fromkeys(df_plot['Cdep']))
        drugs_type_key_range=[drugs_type_key,drugs_type_key]
        ints_key=list(dict.fromkeys(df_plot['int']))[1::]  #remove False  
        
        fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'wspace':0.3,'hspace':0.05}, sharex=True,sharey=False)

        #plot no interaction CE curves
        j=0
        for k_type,d_type_key in zip(kill_type_key,drugs_type_key_range):
            for d_type in d_type_key:
                for i in range(2):
        #        for ints_key in ints_key_pairs:
                    for cd_type in cd_type_key:
                        df = df_plot[df_plot['int']==False] #select interaction
                        df = df[df['kill_type']==k_type]
                        df = df[df['drug_type']==d_type]
                        df = df[df['Cdep']==cd_type]
                        ax=axes.flat[j]
                        #select int_radius setting
                        y = df["Av-NetGrowth"][:-1]
                        x = df["concentration"][:-1]/6
                        ax.plot(x, y,linestyle='solid',label=list(df['n_steps'])[0],color=base_color)
                        #add std
                        y_s = df["Std-NetGrowth"][:-1]
                        ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=base_color)
                        ax.grid(which='major')
                        j+=1
        
        #plot interaction curves
        ints_key_pairs=[['MIC-', 'MIC+'], ['G-', 'G+']]    
        j=0
        for k_type,d_type_key in zip(kill_type_key,drugs_type_key_range):
            for d_type in d_type_key:
                for ints_key_pair in ints_key_pairs:
                    ymin=0
                    ymax=0
                    for cd_type in cd_type_key:
                        for ints_key in ints_key_pair:
                            if ints_key=='MIC-' or ints_key=='G-':
                                colors=colors_range[1]
                            elif ints_key=='MIC+' or ints_key=='G+':
                                colors=colors_range[0]
                            df = df_plot[df_plot['int']==ints_key] #select interaction
                            df = df[df['kill_type']==k_type]
                            df = df[df['drug_type']==d_type]
                            df = df[df['Cdep']==cd_type]
                            df = df[df['int_radius']==int_rad] #always zero for no interaction (False)
                            strengths_key=sorted(list(dict.fromkeys(df['int_strength'])),key=abs)
                            ax=axes.flat[j]
                            c_i=0
                            for s in strengths_key:
                                df_s=df[df['int_strength']==s]
                                y = df_s["Av-NetGrowth"][:-1]
                                x = df_s["concentration"][:-1]/6
                                ax.plot(x, y,label=list(df['n_steps'])[0],color=colors[c_i])
                                #add std
                                y_s = df_s["Std-NetGrowth"][:-1]
                                ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=colors[c_i])
                                ax.grid(which='major')
                                ymin=min(ymin,min(y))
                                ymax=max(ymax,max(y))
                                c_i+=1                        
                        j+=1

        #add subfigure letters   
        labels=['a','b','c','d','e','f','g','h']
        for i in range(8):
            ax=axes.flat[i]
            ax.text(0.9, 0.95, labels[i], fontsize=10,fontweight='bold',transform=ax.transAxes, verticalalignment='top')
                       
        #add x labels
        axes_i=[6,7]
        for j in axes_i:
            ax=axes.flat[j]
            ax.set_xlabel('Drug concentration $[mg$ $ml^{-1}$ $MIC]$')  
        
        # add y labels
        axes_i=[0,1,2,3,4,5,6,7]
        for j in axes_i:
            ax=axes.flat[j]
            ax.set_ylabel(r'$Max(k_{net})$')  
            
        plt.subplots_adjust(bottom=0.05, left=0.1, right=0.80, top=0.95)
        
        # add legend
        legend_elements = [
            Line2D([0], [0], color=colors_range[0][2], label='strong +'),
            Line2D([0], [0], color=colors_range[0][1],  label='medium +'), 
            Line2D([0], [0], color=colors_range[0][0],  label='weak +'), 
            Line2D([0], [0], color=base_color,  label='none (control)'), 
            Line2D([0], [0], color=colors_range[1][0],  label='weak -'), 
            Line2D([0], [0], color=colors_range[1][1],  label='medium -'), 
            Line2D([0], [0], color=colors_range[1][2],  label='strong -')]
        ax=axes.flat[3]
        ax.legend(handles=legend_elements,title='Interaction',ncol=1,loc='lower right', bbox_to_anchor=(1.65, -0.5))
        
        # add text boxes
        props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
        ax=axes.flat[0]
        ax.text(0.3, 1.17, 'MIC interaction',transform=ax.transAxes, verticalalignment='top', bbox=props)
        ax=axes.flat[1]
        ax.text(0.2, 1.17, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
        x=1.32
        y=0.75
        ax=axes.flat[1]
        ax.text(x,y , 'Bactericidal \n drug with \n proportional \n effect', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
        ax=axes.flat[3]
        ax.text(x, 0.94, 'Bactericidal \n drug with \n additive \n effect',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
        ax=axes.flat[5]
        ax.text(x, 0.5, 'Bacteriostatic \n drug with \n proportional \n effect', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
        ax=axes.flat[7]
        ax.text(x, y, 'Bacteriostatic \n drug with \n additive \n effect', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
        
        #save figure
        saveloc = os.getcwd() + "/Figures"+"/Figure6_highdpi_v2.png"
        fig.savefig(saveloc, format="png",dpi=600)
        
#%% plot ci and cd for one drug - Figure "Figure7_highdpi_v2"

mpl.rcParams['figure.figsize']=[6.52,3]

m_type='stepwise'
int_rad=1
n_step=1
for int_rad in [0]:
    for n_step in [1]:            
        #take subset to plot
        df_plot = df_all.copy()
        df_plot = df_plot[df_plot['move_type']==m_type] 
        df_plot = df_plot[df_plot['n_steps']==n_step]
        df_plot=df_plot.reset_index(drop=True)
        
        drugs_key=list(dict.fromkeys(df_plot['drug']))
        drugs_type_key=list(dict.fromkeys(df_plot['drug_type']))
        kill_type_key=[list(dict.fromkeys(df_plot['kill_type']))[0]]
        cd_type_key=list(dict.fromkeys(df_plot['Cdep']))
        drugs_type_key_range=[drugs_type_key[0]]
        ints_key=list(dict.fromkeys(df_plot['int']))[1::]  #remove False  
        
        fig, axes = plt.subplots(nrows=2,ncols=2,gridspec_kw={'wspace':0.02,'hspace':0.05}, sharex=True,sharey=True)
        #plot no interaction CE curves
        j=0
        for cd_type in cd_type_key:
            for k_type,d_type in zip(kill_type_key,drugs_type_key_range):
                for i in range(2):
                    df = df_plot[df_plot['int']==False] #select interaction
                    df = df[df['kill_type']==k_type]
                    df = df[df['drug_type']==d_type]
                    df = df[df['Cdep']==cd_type]
                    ax=axes.flat[j]
                    #select int_radius setting
                    y = df["Av-NetGrowth"]
                    x = df["concentration"]/6
                    ax.plot(x, y,label=list(df['n_steps'])[0],color=base_color)
                    #add std
                    y_s = df["Std-NetGrowth"]
                    ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=base_color)
                    ax.grid(which='major', linewidth=0.3)
                    j+=1
        
        #plot interaction curves
        ints_key_pairs=[['MIC-', 'MIC+'], ['G-', 'G+']]    
        
        j=0
        for cd_type in cd_type_key:
            for k_type,d_type in zip(kill_type_key,drugs_type_key_range):
                for ints_key_pair in ints_key_pairs:
                    if ints_key_pair[0]=='MIC-':
                        colors=colors_range[0]
                    elif ints_key_pair[0]=='G-':
                        colors=colors_range[1]
                    for ints_key in ints_key_pair:
                        if ints_key=='MIC-' or ints_key=='G-':
                            colors=colors_range[1]
                        elif ints_key=='MIC+' or ints_key=='G+':
                            colors=colors_range[0]
                                
                        df = df_plot[df_plot['int']==ints_key] #select interaction
                        df = df[df['kill_type']==k_type]
                        df = df[df['drug_type']==d_type]
                        df = df[df['Cdep']==cd_type]
                        df = df[df['int_radius']==int_rad] #always zero for no interaction (False)
                        strengths_key=sorted(list(dict.fromkeys(df['int_strength'])),key=abs)
                        ax=axes.flat[j]
                        c_i=0
                        for s in strengths_key:
                            df_s=df[df['int_strength']==s]
                            y = df_s["Av-NetGrowth"]
                            x = df_s["concentration"]/6
                            ax.plot(x, y,label=list(df['n_steps'])[0],color=colors[c_i])
                            #add std
                            y_s = df_s["Std-NetGrowth"]
                            ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=colors[c_i])
                            ax.grid(which='major')
                            c_i+=1                        
                    j+=1
        
        # add subfigure labels
        labels=['a','b','c','d','e','f','g','h']
        for i in range(4):
            ax=axes.flat[i]
            ax.text(0.9, 0.95, labels[i],fontweight='bold',transform=ax.transAxes, verticalalignment='top')
        
        #add x labels
        axes_i=[2,3]
        for j in axes_i:
            ax=axes.flat[j]
            ax.set_xlabel('Drug concentration $[mg$ $ml^{-1}$ $MIC]$')  
        
        #add y labels
        axes_i=[0,2]
        for j in axes_i:
            ax=axes.flat[j]
            ax.set_ylabel(r'$Max(k_{net})$')  
        
        plt.subplots_adjust(bottom=0.11, left=0.1, right=0.77, top=0.9)        
        
        # add legend
        legend_elements = [
            Line2D([0], [0], color=colors_range[0][2], label='strong +'),
            Line2D([0], [0], color=colors_range[0][1],  label='medium +'), 
            Line2D([0], [0], color=colors_range[0][0],  label='weak +'), 
            Line2D([0], [0], color=base_color,  label='none (control)'), 
            Line2D([0], [0], color=colors_range[1][0],  label='weak -'), 
            Line2D([0], [0], color=colors_range[1][1],  label='medium -'), 
            Line2D([0], [0], color=colors_range[1][2],  label='strong -')]
        ax=axes.flat[3]
        ax.legend(handles=legend_elements,title='Interaction',ncol=1,loc='lower right', bbox_to_anchor=(1.65,0.4))
        
        # place text boxes
        props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
        ax=axes.flat[0]
        ax.text(0.3, 1.18, 'MIC interaction',transform=ax.transAxes, verticalalignment='top', bbox=props)
        ax=axes.flat[1]
        ax.text(1.35, 1, 'Concentration- \n independent drug',transform=ax.transAxes, horizontalalignment='center', verticalalignment='top', bbox=props)
        ax=axes.flat[1]
        ax.text(0.25, 1.18, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
        ax=axes.flat[3]
        ax.text(1.35, 0.4, 'Concentration- \n dependent drug',transform=ax.transAxes, horizontalalignment='center',verticalalignment='top', bbox=props)
        
        saveloc = os.getcwd() + "/Figures"+"/Figure7_highdpi_v2.png"        
        fig.savefig(saveloc, format="png",dpi=600)
        
