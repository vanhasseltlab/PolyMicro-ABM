#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 13:27:12 2022
Plotting script
Author: Catharina Herzberg, LACDR, Leiden University, The Netherlands
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
import matplotlib.pyplot as plt
import seaborn as sns

#%% PART 3
# load df_all_s from pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_s_intrad1|0_v16&v17_VI.pkl", "rb")
df_all_s = pickle.load(file)
file.close()

#%%
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

mpl.rcParams['figure.figsize'] = [3.2,4]

n_steps_range=[0,1] # n_steps=0 or structured lifestyle and n_steps=1 for unstructured lifestyle
for n_steps in n_steps_range:
    # plot k_net Logchange as heatmap Interaction radius comparison
    df_dict_CE_rad=df_all_s.copy()
    df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['move_type']=='stepwise']
    df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['n_steps']==n_steps]
    df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['Cdep']=='ci']
    df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['drug_type']=='prop']
    
    
    fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'wspace':0.05,'hspace':0.8}, sharex=False,sharey=False)
    cbar_ax = fig.add_axes([.27, .08, 0.7, 0.03]) #add axis for scale bar 
     
    f=0
    maxv=0
    minv=0
    for int_type_pair in [['MIC+', 'G+'],['MIC-','G-']]:
        for kill_type in list(dict.fromkeys(df_dict_CE_rad['kill_type'])):
            for drug_type in list(dict.fromkeys(df_dict_CE_rad['drug_type'])):
                for int_type in int_type_pair:
                    for c_dep in list(dict.fromkeys(df_dict_CE_rad['Cdep'])):
    
                        # pandas pivot
                        ax=axes.flat[f]
                        heatmap1_data = pd.pivot_table(df_dict_CE_rad[df_dict_CE_rad["int"]==int_type][df_dict_CE_rad["Cdep"]==c_dep][df_dict_CE_rad["kill_type"]==kill_type][df_dict_CE_rad["drug_type"]==drug_type], values='LC-NetGrowth',
                                                       index=["int_strength_string" ],
                                                       columns=["concentration_simple" ],dropna=False)
                        #add sign to int_strength_string in heatmap
                        heatmap1_data=heatmap1_data.rename({'weak':'weak'+int_type[-1],'medium':'medium'+int_type[-1],'strong':'strong'+int_type[-1]})
                        
                        maxv=max(maxv,np.nanmax(np.array(heatmap1_data)))
                        minv=max(minv,np.nanmin(np.array(heatmap1_data)),key=abs)
                        colormap='seismic_r'
                        vmin=-1.85
                        vmax=1.85
                        if f==0:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=True,cbar_ax=cbar_ax,cbar_kws={'label': 'Log fold change','orientation':'horizontal'},vmin=vmin, vmax=vmax,xticklabels=False)
                            axes.flat[f].set_xlabel('')
                            axes.flat[f].set_ylabel('')
                        elif f==2:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
                            axes.flat[f].set_xlabel('')
                            axes.flat[f].set_ylabel('')
                        elif f==4:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
                            axes.flat[f].set_xlabel('')
                            axes.flat[f].set_ylabel('                           Interaction strength \n                          negative (-)                               positive (+)')
                        elif f==6:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
                            axes.flat[f].set_ylabel('')
                            axes.flat[f].set_xlabel('Drug concentration')
                            axes.flat[f].annotate('', xy=(-0.02, -0.05), xycoords='axes fraction', xytext=(1.025, -0.05), arrowprops=dict(arrowstyle="<-", color='black',linewidth=0.6))
    
                        elif f==7:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,yticklabels=False,xticklabels=False)
                            axes.flat[f].set_ylabel('')
                            axes.flat[f].set_xlabel('Drug concentration')
                            axes.flat[f].annotate('', xy=(-0.02, -0.05), xycoords='axes fraction', xytext=(1.025, -0.05), arrowprops=dict(arrowstyle="<-", color='black',linewidth=0.6))
    
    
                        else:
                            g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,yticklabels=False,xticklabels=False)
                            axes.flat[f].set_ylabel('')
                            axes.flat[f].set_xlabel('')
                        
                        f+=1
                              
    plt.subplots_adjust(bottom=0.2, left=0.24, right=0.99, top=0.88)        
    
    # place text boxes 
    props = dict(boxstyle='round', facecolor='grey', alpha=0.1)
    ax=axes.flat[0]
    y=1.95
    ax.text(0.2, y, 'MIC interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
    ax=axes.flat[1]
    ax.text(0, y, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    x=-0.05
    y=1.4
    ax=axes.flat[1]
    ax.text(x, y, 'Bactericidal drug with proportional effect',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[3]
    ax.text(x, y, 'Bacteriostatic drug  with proportional effect',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[5]
    ax.text(x, y, 'Bactericidal drug  with proportional effect',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    ax=axes.flat[7]
    ax.text(x, y, 'Bacteriostatic drug  with proportional effect', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
    
    ax=axes.flat[0]
    if n_steps==0:
        ax.text(-0.6,-6.5, '\n Structured \n   lifestyle \n',transform=ax.transAxes, verticalalignment='center', bbox=props)
        saveloc = os.getcwd() +"/Figures"+"/Figure10-1_highdpi_v2.png"
    else:
        ax.text(-0.6,-6.5, '\nUnstructured\n    lifestyle \n',transform=ax.transAxes, verticalalignment='center', bbox=props)
        saveloc = os.getcwd() +"/Figures"+"/Figure10-2_highdpi_v2.png"
    
    plt.savefig(saveloc, format="png",dpi=600)

#%%
# #alternative heatmap for two columns
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['ps.fonttype'] = 42
# mpl.rcParams['font.family'] = ['Arial']
# mpl.rcParams['figure.dpi']= 600
# mpl.rcParams['figure.figsize']=[6.52,4]

# mpl.rcParams['lines.linewidth']  = 0.5
# mpl.rcParams['legend.fontsize'] = 8
# mpl.rcParams['legend.title_fontsize'] = 8
# mpl.rcParams['axes.labelsize']=8
# mpl.rcParams['axes.titlesize']=8
# mpl.rcParams['font.size']=8
# mpl.rcParams['xtick.labelsize']=6 
# mpl.rcParams['ytick.labelsize']=6
# mpl.rcParams['axes.linewidth'] = 0.5
# mpl.rcParams['xtick.major.size'] = 0.5
# mpl.rcParams['xtick.major.width'] = 0.5
# mpl.rcParams['xtick.minor.size'] = 0.5
# mpl.rcParams['xtick.minor.width'] = 0.3
# mpl.rcParams['ytick.major.size'] = 1
# mpl.rcParams['ytick.major.width'] = 0.5
# mpl.rcParams['ytick.minor.size'] = 0.5
# mpl.rcParams['ytick.minor.width'] = 0.3
# mpl.rcParams['grid.linewidth'] = 0.3
# mpl.rcParams['lines.markersize'] = 0.5
# mpl.rcParams['lines.marker'] = 'o' 

# # plot k_net Logchange as heatmap Interaction radius comparison
# df_dict_CE_rad=df_all_s.copy()
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['move_type']=='stepwise']
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['n_steps']==0]
# # df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['n_steps']==1]
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['Cdep']=='ci']
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['drug_type']=='prop']


# # df_dict_CE_rad['Cdep'] = pd.Categorical(df_dict_CE_rad.Cdep,categories=['ci','cd'],ordered=True)

# #df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad["conc"]==7.0]
# # mpl.rcParams['ytick.labelsize']=10
# # mpl.rcParams['ytick.labelsize']=10
# # mpl.rcParams['axes.labelsize']=10
# # mpl.rcParams['axes.titlesize']=14

# fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'wspace':0.03,'hspace':0.1}, sharex=False,sharey=False)
# # cbar_ax = fig.add_axes([.71, .3, .03, .4]) #add axis for scale bar
# cbar_ax = fig.add_axes([.76, .3, .03, .4]) #add axis for scale bar
 
# f=0
# maxv=0
# minv=0
# for int_type_pair in [['MIC+', 'G+'],['MIC-','G-']]:
#     for kill_type in list(dict.fromkeys(df_dict_CE_rad['kill_type'])):
#         for drug_type in list(dict.fromkeys(df_dict_CE_rad['drug_type'])):
#             for int_type in int_type_pair:
#                 for c_dep in list(dict.fromkeys(df_dict_CE_rad['Cdep'])):

#                     # pandas pivot
#                     ax=axes.flat[f]
#                     heatmap1_data = pd.pivot_table(df_dict_CE_rad[df_dict_CE_rad["int"]==int_type][df_dict_CE_rad["Cdep"]==c_dep][df_dict_CE_rad["kill_type"]==kill_type][df_dict_CE_rad["drug_type"]==drug_type], values='LC-NetGrowth',
#                                                    index=["int_strength_string" ],
#                                                    columns=["concentration_simple" ],dropna=False)
#                     #add sign to int_strength_string in heatmap
#                     heatmap1_data=heatmap1_data.rename({'weak':'weak'+int_type[-1],'medium':'medium'+int_type[-1],'strong':'strong'+int_type[-1]})
                    
#                     maxv=max(maxv,np.nanmax(np.array(heatmap1_data)))
#                     minv=max(minv,np.nanmin(np.array(heatmap1_data)),key=abs)
#                     colormap='coolwarm_r'
#                     colormap='seismic_r'
#                     vmin=-1.85
#                     vmax=1.85
#                     if f==0:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=True,cbar_ax=cbar_ax,cbar_kws={'label': 'Log fold change'},vmin=vmin, vmax=vmax,xticklabels=False)
#                         # axes.flat[f].set_ylabel('Positive interaction                    \n \n \n Cidal \n')
#                         axes.flat[f].set_xlabel('')
#                         axes.flat[f].set_ylabel('')
#                     elif f==2:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
#                         # axes.flat[f].set_ylabel('                                        \n \n \n Static \n')
#                         axes.flat[f].set_xlabel('')
#                         axes.flat[f].set_ylabel('')
#                     elif f==4:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
#                         # axes.flat[f].set_ylabel('Negative interaction                    \n \n \n Cidal \n')
#                         axes.flat[f].set_xlabel('')
#                         axes.flat[f].set_ylabel('                   Interaction strength')
#                     elif f==6:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,xticklabels=False)
#                         axes.flat[f].set_ylabel('')
#                         # axes.flat[f].set_ylabel('                                        \n \n \n Static \n')
#                         axes.flat[f].set_xlabel('Drug concentration')
#                     elif f==7:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,yticklabels=False,xticklabels=False)
#                         axes.flat[f].set_ylabel('')
#                         # axes.flat[f].set_ylabel('                                        \n \n \n Static \n')
#                         axes.flat[f].set_xlabel('Drug concentration')

#                     else:
#                         g=sns.heatmap(heatmap1_data, cmap=colormap,ax=axes.flat[f],cbar=False,vmin=vmin, vmax=vmax,yticklabels=False,xticklabels=False)
#                         axes.flat[f].set_ylabel('')
#                         axes.flat[f].set_xlabel('')
                    
#                     f+=1
                    
# # plt.subplots_adjust(bottom=0.05, left=0.26, right=0.70, top=0.85)  
# plt.subplots_adjust(bottom=0.05, left=0.1, right=0.75, top=0.87)        

# # axes.flat[0].set_title('Resistance interaction \n')
# # axes.flat[1].set_title('Growth interaction \n')
# # fig.suptitle('            Resistance interactions                                                 Growth interactions')
# # fig.tight_layout(rect=[0, 0, .9, 1])
# props = dict(boxstyle='round', facecolor='grey', alpha=0.1)

# ax=axes.flat[0]
# # place a text box in upper left in axes coords
# ax.text(0.25, 1.25, 'MIC interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
# ax.text(-0.25, 1.58, '                                                                             No movement                                                                               ',transform=ax.transAxes, verticalalignment='top', bbox=props)
# # ax.text(-0.25, 1.58, '                                                                           Fast movement                                                                               ',transform=ax.transAxes, verticalalignment='top', bbox=props)

# ax=axes.flat[1]
# ax.text(0.17, 1.25, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)

# # x=-0.65
# # y=2
# # ax=axes.flat[2]
# # ax.text(x,y, '\n \n \n \n \n Positive \n interactions \n \n \n \n \n ',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)

# # ax=axes.flat[6]
# # ax.text(x,y, '\n \n \n \n \n Negative \n interactions \n \n \n \n \n ', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)


# x=1.55
# y=0.75
# ax=axes.flat[1]
# ax.text(x, y, 'Bactericidal \n drug',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
# ax=axes.flat[3]
# ax.text(x, y, 'Bacteriostatic \n drug',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
# ax=axes.flat[5]
# ax.text(x, y, 'Bactericidal \n drug',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)
# ax=axes.flat[7]
# ax.text(x, y, 'Bacteriostatic \n drug', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)




# saveloc = os.getcwd() +"/Simulations/Software Supplementary/Figures"+"/Figure10-1_highdpi.png"
# # saveloc = os.getcwd() +"/Simulations/Software Supplementary/Figures"+"/Figure10-2_highdpi.png"
# saveloc = os.getcwd() + "/Figures"+"/Figure10-1_highdpi.png"

# plt.savefig(saveloc, format="png",dpi=600)