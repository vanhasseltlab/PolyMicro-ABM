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
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
# working directory: '/home/cathi/Simulations/Software Supplementary'
# dir_path = os.path.dirname(os.path.abspath("__file__"))
# os.chdir(dir_path+'/../..')

#%%
# define plot formatting settings
plt.style.use('default')
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



#%%
# load combined df_all in pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_v16&v17.pkl", "rb")
df_all = pickle.load(file)
file.close() 

#%% plot Figure8_highdpi 
#Plot an example CE curve with movement and no movement in separate plots

mpl.rcParams['figure.figsize']=[6.52,3]
    
#take subset to plot
df_plot = df_all.copy()
#select movement subset of data
df_plot = df_plot[df_plot['move_type']=='stepwise'] 
df_plot = df_plot[df_plot['drug_type']=='add']
df_plot = df_plot[df_plot['kill_type']=='cidal']
df_plot = df_plot[df_plot['Cdep']=='ci']
int_rad=0
df_plot = df_plot[df_plot['int_radius']==int_rad]
df_plot=df_plot.reset_index(drop=True)

base_color='blue'
colors_range=[['deepskyblue','mediumslateblue','deeppink'],['limegreen','gold','orangered']]

fig, axes = plt.subplots(nrows=2,ncols=2,gridspec_kw={'wspace':0.02,'hspace':0.05}, sharex=True,sharey=True)
#plot no interaction CE curves

j=0
   
df = df_plot[df_plot['int']==False] #select interaction
df = df[df['n_steps']==0]

for i in range(4):
    ax=axes.flat[j]
    #select int_radius setting
    y = df["Av-NetGrowth"]
    x = df["concentration"]/6
    
    ax.plot(x, y,marker='o',markersize=1,linestyle='solid',label=list(df['n_steps'])[0],color=base_color)
    #add std
    y_s = df["Std-NetGrowth"]
    ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=base_color)
    ax.grid(which='major', linewidth=0.3)
    j+=1
    
ints_key_pairs=[['MIC+','MIC-'],['G+','G-']]

 
#plot interaction CE curves movement
j=0
for int_pair in ints_key_pairs:  
    for i in int_pair:
        if i=='MIC-' or i=='G-':
            colors=colors_range[1]
        elif i=='MIC+' or i=='G+':
            colors=colors_range[0]

        df = df_plot[df_plot['n_steps']==1]
        df = df[df['int']==i]
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

#plot interaction CE curves no movement
j=2
for int_pair in ints_key_pairs:
    for i in int_pair:
        if i=='MIC-' or i=='G-':
            colors=colors_range[1]
        elif i=='MIC+' or i=='G+':
            colors=colors_range[0]

        df = df_plot[df_plot['n_steps']==0]
        df = df[df['int']==i]
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

#add row titles
axes_i=[0,2]
labels=['Movement','No movement']
   
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
ax.text(0.3, 1.18, 'MIC interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
ax=axes.flat[1]
ax.text(1.35, 0.95, 'unstructured \n lifestyle', transform=ax.transAxes, horizontalalignment='center', verticalalignment='top', bbox=props)
ax=axes.flat[1]
ax.text(0.25, 1.18, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)
ax=axes.flat[3]
ax.text(1.35, 0.3, 'structured \n lifestyle', transform=ax.transAxes, horizontalalignment='center',verticalalignment='top', bbox=props)
        
#save figure
saveloc = os.getcwd() + "/Figures"+"/Figure8_highdpi.png"
fig.savefig(saveloc, format="png",dpi=600)

#%% 
# load df_all_s from pickle file
path_all = os.getcwd() + "/Data"
file = open(path_all + "/df_all_s_mvnomv_v16&v17_VI.pkl", "rb")
df_all_s = pickle.load(file) 
file.close() 

#%% plot Figure9_highdpi_v5
# Heatmap only ci drugs and only prop 
# plot k_net Logchange as heatmap Movement comparison
mpl.rcParams['figure.figsize'] = [3.2,4]

df_dict_CE_rad=df_all_s.copy()
df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['int_radius']==0]
df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['Cdep']=='ci']
df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['drug_type']=='prop']

fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'wspace':0.05,'hspace':0.8}, sharex=False,sharey=False)
cbar_ax = fig.add_axes([.27, .08, 0.7, 0.03]) #add axis for scale bar 
f=0
for int_type_pair in [['MIC+', 'G+'],['MIC-','G-']]:
    for kill_type in list(dict.fromkeys(df_dict_CE_rad['kill_type'])):
        for drug_type in list(dict.fromkeys(df_dict_CE_rad['drug_type'])):
            for int_type in int_type_pair:
                for c_dep in list(dict.fromkeys(df_dict_CE_rad['Cdep'])):

                    # pandas pivot
                    ax=axes.flat[f]
                    df_dict_CE_rad['concentration_simple']=df_dict_CE_rad['concentration_simple'].astype(int)
                    heatmap1_data = pd.pivot_table(df_dict_CE_rad[df_dict_CE_rad["int"]==int_type][df_dict_CE_rad["Cdep"]==c_dep][df_dict_CE_rad["kill_type"]==kill_type][df_dict_CE_rad["drug_type"]==drug_type], values='LC-NetGrowth',
                                                   index=["int_strength_string" ],
                                                   columns=["concentration_simple"],dropna=False)
                    #add sign to int_strength_string in heatmap
                    heatmap1_data=heatmap1_data.rename({'weak':'weak'+int_type[-1],'medium':'medium'+int_type[-1],'strong':'strong'+int_type[-1]})
                    
                    vmin=-1.85
                    vmax=1.85
                    colormap='seismic_r'

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

# place a text box
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

saveloc = os.getcwd() + "/Figures"+"/Figure9_highdpi_v5.png"
plt.savefig(saveloc, format="png",dpi=600)

#%% alternative heatmap version
# #      # Heatmap only ci drugs and only prop 
# # plot k_net Logchange as heatmap Movement comparison
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


# df_dict_CE_rad=df_all_s.copy()
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['int_radius']==0]
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['Cdep']=='ci']
# df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad['drug_type']=='prop']

# # df_dict_CE_rad['Cdep'] = pd.Categorical(df_dict_CE_rad.Cdep,categories=['ci','cd'],ordered=True)

# #df_dict_CE_rad=df_dict_CE_rad[df_dict_CE_rad["conc"]==7.0]

# fig, axes = plt.subplots(nrows=4,ncols=2,gridspec_kw={'wspace':0.03,'hspace':0.05}, sharex=False,sharey=False)
# cbar_ax = fig.add_axes([.74, .3, .03, .4]) #add axis for scale bar
 
# f=0
# for int_type_pair in [['MIC+', 'G+'],['MIC-','G-']]:
#     for kill_type in list(dict.fromkeys(df_dict_CE_rad['kill_type'])):
#         for drug_type in list(dict.fromkeys(df_dict_CE_rad['drug_type'])):
#             for int_type in int_type_pair:
#                 for c_dep in list(dict.fromkeys(df_dict_CE_rad['Cdep'])):

#                     # pandas pivot
#                     ax=axes.flat[f]
#                     df_dict_CE_rad['concentration_simple']=df_dict_CE_rad['concentration_simple'].astype(int)
#                     heatmap1_data = pd.pivot_table(df_dict_CE_rad[df_dict_CE_rad["int"]==int_type][df_dict_CE_rad["Cdep"]==c_dep][df_dict_CE_rad["kill_type"]==kill_type][df_dict_CE_rad["drug_type"]==drug_type], values='LC-NetGrowth',
#                                                    index=["int_strength_string" ],
#                                                    columns=["concentration_simple"],dropna=False)
#                     #add sign to int_strength_string in heatmap
#                     heatmap1_data=heatmap1_data.rename({'weak':'weak'+int_type[-1],'medium':'medium'+int_type[-1],'strong':'strong'+int_type[-1]})
                    
#                     vmin=-1.85
#                     vmax=1.85
#                     colormap='coolwarm_r'
#                     colormap='seismic_r'

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
#                         axes.flat[f].set_ylabel('                           Interaction strength \n                          negative (-)                               positive (+)')
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

# # labels=['proportional drug effect','proportional drug effect','additive drug effect','additive drug effect','proportional drug effect','proportional drug effect','additive drug effect','additive drug effect']
# # for i in range(8):
# #     ax=axes.flat[i]
# #     ax.set_title(labels[i])

                    
# plt.subplots_adjust(bottom=0.05, left=0.12, right=0.73, top=0.87)        
# # axes.flat[0].set_title('Resistance interaction \n')
# # axes.flat[1].set_title('Growth interaction \n')
# # fig.suptitle('            Resistance interactions                                                 Growth interactions')
# # fig.tight_layout(rect=[0, 0, .9, 1])
# props = dict(boxstyle='round', facecolor='grey', alpha=0.1)

# ax=axes.flat[0]
# # place a text box in upper left in axes coords
# ax.text(0.25, 1.3, 'MIC interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)

# ax=axes.flat[1]
# ax.text(0.17, 1.3, 'Growth interaction', transform=ax.transAxes, verticalalignment='top', bbox=props)

# # x=-0.65
# # y=1.95
# # ax=axes.flat[2]
# # ax.text(x,y, '\n \n \n \n \n \n Positive \n interactions \n \n \n \n \n ',transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)

# # ax=axes.flat[6]
# # ax.text(x,y, '\n \n \n \n \n \n Negative \n interactions \n \n \n \n \n ', transform=ax.transAxes,horizontalalignment='center', verticalalignment='top', bbox=props)


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



# #    #        ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)

# # saveloc = "/home/cathi/Simulations/Simu_PW_v16_grid/Comparing no movement to movement/Heatmap_LC_NetGrowth_VI"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Heatmap_ci_LC_NetGrowth_VI_intrad0.png"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Heatmap_ci_prop_highdpi_LC_NetGrowth_VI_intrad0_Improved.png"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Figure9_highdpi.png"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Figure9_highdpi_v2.png"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Figure9_highdpi_v3.png"
# saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Figure9_highdpi_v4.png"
# saveloc = os.getcwd() + "/Simulations/Software Supplementary/Figures"+"/Figure9_highdpi_v4.png"
# saveloc = os.getcwd() + "/Figures"+"/Figure9_highdpi_v4TEST.png"

# # saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Heatmap_LC_NetGrowth_VI_intrad1"
# # saveloc = "/home/cathi/Simulations/Simu_PW_v17/Part2/Heatmap_LC_NetGrowth_VI_intrad1_scale"

# #
# plt.savefig(saveloc, format="png",dpi=600)