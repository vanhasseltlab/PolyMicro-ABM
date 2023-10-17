#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 18:14:06 2022
Plotting script
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
"""
import os
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#%%

path_core=os.getcwd() + "/Data"
file = open(path_core + "/df_metrics_3.pkl", "rb")
df_metrics = pickle.load(file)
file.close()

#%% plot formatting
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = ['Arial']
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

#%% for one single cidal drug : mean and std in separate plots
#plot Mean_its(Mean_pop(Target))=:Targetmeanmean over concentrations with Mean_its(Std_pop(Target))=:Targetstdmean

# select subset and calculate metrics for plot
df_plot=df_metrics.copy()
# df_plot=df_plot.sample(n=100)
#flatten columns
indexes=[x[0]+x[1] for x in df_plot.columns.to_flat_index()]
df_plot.columns=indexes
#calculate Mean_its(Targetmean) and Mean_its(Targetstd)
df_plot_metrics=df_plot.groupby(['Species',
  'drug_type',
  'n_steps',
  'int',
  'int_radius',
  'int_strength',
  'move_type',
  'concentration',
  'Cdep',
  'kill_type',
  'concentration_simple',
  'int_strength_string',
  'move_type_simple'])['Targetmean','Targetstd'].describe()
df_plot_metrics.reset_index(inplace=True)
indexes=[x[0]+x[1] for x in df_plot_metrics.columns.to_flat_index()]
df_plot_metrics.columns=indexes
df_plot_metrics=df_plot_metrics.drop(columns=['Targetmeancount', 'Targetmeanmin',
'Targetmean25%', 'Targetmean50%', 'Targetmean75%', 'Targetmeanmax',
'Targetstdcount', 'Targetstdmin',
'Targetstd25%', 'Targetstd50%', 'Targetstd75%', 'Targetstdmax'])

#%% plot only means portrait: mean plot for other settings: n_steps and int_rad
#plot

mpl.rcParams['figure.figsize'] = [3.2,3.2]

m_type='stepwise'
int_rad=0
n_step=1
# drug_type='add'
drug_type='prop'
cd_type='ci'
kill_type='cidal'
# kill_type='static'

#take subset to plot
df_plot = df_plot_metrics.copy()
#select subset of data
df_plot = df_plot[df_plot['move_type']==m_type] 
df_plot = df_plot[df_plot['n_steps']==n_step]
df_plot = df_plot[df_plot['int_radius']==int_rad]
df_plot = df_plot[df_plot['drug_type']==drug_type]
df_plot = df_plot[df_plot['kill_type']==kill_type]
df_plot = df_plot[df_plot['Cdep']==cd_type]

df_plot=df_plot.reset_index(drop=True)

drugs_type_key=list(dict.fromkeys(df_plot['drug_type']))
kill_type_key=list(dict.fromkeys(df_plot['kill_type']))
cd_type_key=list(dict.fromkeys(df_plot['Cdep']))

colors_range=[['deepskyblue','mediumslateblue','deeppink'],['limegreen','gold','orangered']]
base_color='blue'

mpl.rcParams['figure.figsize']=[3.2,3.2]

fig, axes = plt.subplots(nrows=2,ncols=1,gridspec_kw={'wspace':0,'hspace':0.3}, sharex=True,sharey=False)

ints_key_pairs=[['MIC-', 'MIC+'], ['G-', 'G+']]    

j=0
for k_type in kill_type_key:
    df_p = df_plot[df_plot['kill_type']==k_type]
    for i in range(1):
        for ints_key_pair in ints_key_pairs:
            ax=axes.flat[j]
            for ints_key in ints_key_pair:
                df = df_p[df_p['int']==ints_key] #select interaction
                
                strengths_key=['weak','medium','strong']
                
                if ints_key=='MIC-' or ints_key=='G-':
                    colors=colors_range[1]
                elif ints_key=='MIC+' or ints_key=='G+':
                    colors=colors_range[0]
                if i==0:
                    if ints_key=='MIC-' or ints_key=='MIC+':
                        y=[6]*9
                        x=df[df['int_strength_string']=='weak']["concentration"]/6
                        
                        ax.plot(x, y,marker='o',linestyle='solid',color=base_color)
                    elif ints_key=='G-' or ints_key=='G+':
                        y=[0.015]*9
                        x=df[df['int_strength_string']=='weak']["concentration"]/6
                        ax.plot(x, y,color=base_color)
    
                c_i=0
                for s in strengths_key:
                    df_s=df[df['int_strength_string']==s]
                    if i==0:
                        y = df_s["Targetmeanmean"]
                        x = df_s["concentration"]/6
                        ax.plot(x, y,color=colors[c_i])
                    if i==1:
                        #add std
                        y = df_s["Targetstdmean"]
                        x = df_s["concentration"]/6
                        ax.plot(x, y,color=colors[c_i])
                        y=[0]*9
                        x=df_s["concentration"]/6
                        ax.plot(x, y,color=base_color)
                    ax.grid(which='major', linewidth=0.3)
                    c_i+=1
            j+=1

#add x labels
axes_i=[1]
for j in axes_i:
    ax=axes.flat[j]
    ax.set_xlabel('Drug concentration $[mg$ $ml^{-1}$ $MIC]$')  

#add y labels
axes_i=[0,1]

labels= [r'Average $\overline{MIC}$',r'Average $\overline{k_{rep,max}}$']
for j in axes_i:
    ax=axes.flat[j]
    ax.set_ylabel(labels[j])  

props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
#adjust yaxes
axes.flat[0].set_ylim(3.5,8.5)
axes.flat[1].set_ylim(0,0.035)

labels=['d','e','f','g','h']
for i in range(2):
    ax=axes.flat[i]
    ax.text(0.9, 0.9, labels[i],fontweight='bold',transform=ax.transAxes, verticalalignment='top')
        
# add legend
plt.subplots_adjust(bottom=0.1, left=0.15, right=0.62, top=0.92)        
legend_elements = [
    Line2D([0], [0], color=colors_range[0][2], label='strong +'),
    Line2D([0], [0], color=colors_range[0][1],  label='medium +'), 
    Line2D([0], [0], color=colors_range[0][0],  label='weak +'), 
    Line2D([0], [0], color=base_color,  label='none (control)'), 
    Line2D([0], [0], color=colors_range[1][0],  label='weak -'), 
    Line2D([0], [0], color=colors_range[1][1],  label='medium -'), 
    Line2D([0], [0], color=colors_range[1][2],  label='strong -')]
ax=axes.flat[1]
ax.legend(handles=legend_elements,title='Interaction',ncol=1,loc='lower right', bbox_to_anchor=(1.84,0.5))
props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
ax=axes.flat[0]
# place a text box in upper left in axes coords
y=1.15
x=0.2
ax.text(x, y, 'MIC interaction',transform=ax.transAxes, verticalalignment='top', bbox=props)
ax=axes.flat[1]
ax.text(x, y, 'Growth interaction',transform=ax.transAxes, verticalalignment='top', bbox=props)
   
saveloc = os.getcwd() + "/Figures/"+"/Figure4_highdpi.png"

fig.savefig(saveloc, format="png",dpi=600)

#%% plot histogram of target heterogeneity for one interaction for each concentration - part 1
    
file = open(path_core + "/df_agents_part_3.pkl", "rb")
df = pickle.load(file)
file.close()    
df=df[df['int_radius']==0]
df=df[df['n_steps']==1]
df=df[df['drug_type']=='prop']
df=df[df['kill_type']=='cidal']
df=df[df['Cdep']=='ci']
df=df[df['int']=='MIC+']

mpl.rcParams['figure.figsize']=[3.2,1.2]

int_strengths=list(dict.fromkeys(df['int_strength_string']))
int_strengths=['medium']

concs=list(dict.fromkeys(df['concentration_simple']))
concs=[0,3,6]

# for c in concs:
fig, axes = plt.subplots(nrows=len(int_strengths),ncols=len(concs),gridspec_kw={'wspace':0.1,'hspace':0.2}, sharex=True,sharey=True)
j=0
for i in int_strengths:
    df_plot=df.copy()
    df_plot=df_plot[df_plot['int_strength_string']==i]
    for conc in concs:
        df_p=df_plot.copy()
        df_p=df_p[df_p['concentration_simple']==conc]
        its=list(dict.fromkeys(df_p['iteration']))
        df_list=list(df_p['Target']) #all iterations
        
        ax=axes.flat[j]

        maxdf=max(df_plot['Target'])
        mindf=min(df_plot['Target'])

        binwidth=0.05*(maxdf-mindf)
        binwidth=0.02*(maxdf-mindf)
        ax.hist(df_list, bins=np.arange(mindf, maxdf + binwidth, binwidth),color='deepskyblue')#,log=True)
        
        # ax.hist(df_list,color='deepskyblue',log=True)
        ax.axvline(np.mean(df_list), color='r', linestyle='dashed', linewidth=0.5)

        ax.set_axisbelow(True)
        ax.grid(which='major', linewidth=0.2)
        

        j+=1

props = dict(boxstyle='round', facecolor='grey', alpha=0.1)

#add x labels
axes_i=[0,1,2]
for j in axes_i:
    ax=axes.flat[j]
    ax.set_xlabel('MIC')
    ax.set_xticks(np.arange(6, 18, step=4))  # Set label locations.

axes_i=[0]
for j in axes_i:
    ax=axes.flat[j]
    ax.set_ylabel('Number of \n bacteria')

labels=['a','b','c','d']
for i in range(3):
    ax=axes.flat[i]
    ax.text(0.85, 0.95, labels[i], fontsize=8,fontweight='bold',transform=ax.transAxes, verticalalignment='top')

labels=['c=1','c=1.45','c=4']
# [6, 6.9, 7.8, 8.7, 9.9, 12.6, 18, 24, 36]
for i in range(3):
    ax=axes.flat[i]
    ax.text(0.5, 1.25, labels[i], fontsize=8,transform=ax.transAxes, verticalalignment='top',horizontalalignment='center')

plt.subplots_adjust(bottom=0.25, left=0.2, right=0.95, top=0.8)  

saveloc = os.getcwd() + "/Figures/"+"/Figure3B_highdpi_allitsTEST.png"

fig.savefig(saveloc, format="png",dpi=600)    

#%% plot histogram of target heterogeneity
    
file = open(path_core + "/df_agents_part_3.pkl", "rb")
df = pickle.load(file)
file.close()    
df=df[df['int_radius']==0]
df=df[df['n_steps']==1]
df=df[df['drug_type']=='prop']
df=df[df['kill_type']=='cidal']
df=df[df['Cdep']=='ci']
df=df[df['concentration_simple']==4]

mpl.rcParams['figure.figsize']=[6.52,3.2]


min(df[df['int']=='MIC-']['Target'])
max(df[df['int']=='MIC-']['Target'])
min(df[df['int']=='MIC+']['Target'])
max(df[df['int']=='MIC+']['Target'])
min(df[df['int']=='G-']['Target'])
max(df[df['int']=='G-']['Target'])
min(df[df['int']=='G+']['Target'])
max(df[df['int']=='G+']['Target'])

#%
int_strengths=list(dict.fromkeys(df['int_strength_string']))
concs=list(dict.fromkeys(df['concentration_simple']))
for c in concs:
    fig, axes = plt.subplots(nrows=len(int_strengths),ncols=4,gridspec_kw={'wspace':0.1,'hspace':0.2}, sharex=False,sharey=True)
    j=0
    for i in int_strengths:
        df_plot=df.copy()
        df_plot=df_plot[df_plot['int_strength_string']==i]
        for int_key in ['MIC-','MIC+','G-','G+']:
            df_p=df_plot.copy()
            df_p=df_p[df_p['int']==int_key]
            its=list(dict.fromkeys(df_p['iteration']))
            df_list=list(df_p[df_p['iteration']==1][df_p['concentration_simple']==c]['Target'])
            df_list=list(df_p[df_p['concentration_simple']==c]['Target']) #all iterations
            
            ax=axes.flat[j]
            if int_key=='MIC-':
                maxdf=max(df[df['int']=='MIC-']['Target'])
                mindf=min(df[df['int']=='MIC-']['Target'])
                col='deepskyblue'
            elif int_key=='MIC+':
                maxdf=max(df[df['int']=='MIC+']['Target'])
                mindf=min(df[df['int']=='MIC+']['Target'])
                col='deepskyblue'
            elif int_key=='G-':
                maxdf=max(df[df['int']=='G-']['Target'])
                mindf=min(df[df['int']=='G-']['Target'])
                col= 'deeppink'
            elif int_key=='G+':
                maxdf=max(df[df['int']=='G+']['Target'])
                mindf=min(df[df['int']=='G+']['Target'])
                col='deeppink'
            binwidth=0.05*(maxdf-mindf)
            binwidth=0.01*(maxdf-mindf)
            ax.hist(df_list, bins=np.arange(mindf, maxdf + binwidth, binwidth),color=col,log=True)
            ax.set_axisbelow(True)
            ax.grid(which='major', linewidth=0.2)
            

            j+=1
    props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
    y=1.5
    x=0.5
    ax=axes.flat[0]
    ax.text(x,y, 'Negative MIC \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[1]
    ax.text(x, y, 'Positive MIC \n interaction',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[2]
    ax.text(x,y, 'Negative Growth \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[3]
    ax.text(x, y, 'Positive Growth \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    y=0.6
    x=-0.72
    ax=axes.flat[0]
    ax.text(x,y, 'weak', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[4]
    ax.text(x, y, 'medium',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[8]
    ax.text(x,y, 'strong', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    
    #add x labels
    axes_i=[8,9]
    for j in axes_i:
        ax=axes.flat[j]
        ax.set_xlabel('MIC')
        ax.xaxis.label.set_color('deepskyblue')
    axes_i=[10,11]
    for j in axes_i:
        ax=axes.flat[j]
        ax.set_xlabel(r'$k_{rep,max}$')
        ax.xaxis.label.set_color('deeppink')
    
    axes_i=[0,4,8]
    for j in axes_i:
        ax=axes.flat[j]
        ax.set_ylabel('Number of \n bacteria')
    
    # set ylim 
    ax=[0,4,8]
    for x in ax:
        axes.flat[x].set_ylim([pow(10,0),pow(10,4)])
    ax=[1,5,9]
    for x in ax:
        axes.flat[x].set_ylim([pow(10,0),pow(10,4)])

    ax=[2,6,10]
    for x in ax:
        axes.flat[x].set_ylim([pow(10,0),pow(10,4)])

    ax=[3,7,11]
    for x in ax:
        axes.flat[x].set_ylim([pow(10,0),pow(10,4)])

    ax=range(8)
    for x in ax:
        ax=axes.flat[x]
        ax.get_xaxis().set_ticklabels([])
    
    labels=['a','b','c','d','e','f','g','h','i','j','k','l']
    for i in range(12):
        ax=axes.flat[i]
        ax.text(0.9, 0.95, labels[i], fontsize=8,fontweight='bold',transform=ax.transAxes, verticalalignment='top')
    
    plt.subplots_adjust(bottom=0.1, left=0.18, right=0.99, top=0.85)  
    
    saveloc = os.getcwd() + "/Figures/"+"/Figure3_highdpi_v4.png"

    fig.savefig(saveloc, format="png",dpi=600)
#%% plot avergae CV as histogram

file = open(path_core + "/df_metrics_3.pkl", "rb")
df_metrics = pickle.load(file)
file.close()

indexes=[x[0]+x[1] for x in df_metrics.columns.to_flat_index()]
df_metrics.columns=indexes


mpl.rcParams['figure.figsize'] = [3.2,3.2]


df_meanCV=df_metrics.groupby(['int','int_strength_string'])['TargetCV'].describe()
df_meanCV.reset_index(inplace=True)
df_meanCV=df_meanCV[['int','int_strength_string','mean','std']]
ints=['MIC+','G+','MIC-','G-']

fig, axes = plt.subplots(nrows=2,ncols=2,gridspec_kw={'wspace':0,'hspace':0.06},sharex=True,sharey=True)
j=0
colors_range=[['deepskyblue','mediumslateblue','deeppink'],['limegreen','gold','orangered']]

for i in ints:
    ax=axes.flat[j]
    x=df_meanCV[df_meanCV['int']==i]['int_strength_string']
    mn=df_meanCV[df_meanCV['int']==i]['mean']
    std=df_meanCV[df_meanCV['int']==i]['std']
    ax.grid(which='major', linewidth=0.1)
    if i=='MIC-' or i=='G-':
        colors=colors_range[1]
    else: 
        colors=colors_range[0]
    
    ax.bar(x,mn,color='deepskyblue')
    ax.set_ylim(0,1.6)
    j+=1
    
plt.subplots_adjust(bottom=0.05, left=0.13, right=0.99, top=0.98)        

axes.flat[0].set_ylabel(r'$\overline{CV}$')
axes.flat[2].set_ylabel(r'$\overline{CV}$')
props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
y=0.95
x=0.5
ax=axes.flat[0]
ax.text(x, y, 'Positive MIC \n interaction',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
ax=axes.flat[1]
ax.text(x, y, 'Positive Growth \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
ax=axes.flat[2]
ax.text(x,y, 'Negative MIC \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
ax=axes.flat[3]
ax.text(x,y, 'Negative Growth \n interaction', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)

saveloc = os.getcwd() + "/Figures/"+"/Figure5_highdpi_v3TEST.png"
plt.savefig(saveloc, format="png",dpi=600)