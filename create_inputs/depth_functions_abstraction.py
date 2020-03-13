# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:55:29 2020

@author: engelen
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

import os

def assign_bin_depths(df, depth_col, max_depth):
    df["depth (m)"] = pd.cut(df[depth_col], np.arange(0, max_depth, 20), duplicates="drop")

def sum_Q_per_group(df, group, *Q_cols):
    return(df.groupby(group).sum()[list(Q_cols)])

def count_Q_per_group(df, group, *Q_cols):
    return(df.groupby(group).count()[list(Q_cols)])

def bin_to_mids(index):
    return([(a.left + a.right)/2 for a in index])

def get_Q_per_aqf(df, *Q_cols, func=sum_Q_per_group):
    Q_per_aqf = func(df, "aqf_nr", *Q_cols)
    #Scale
    Q_per_aqf /= Q_per_aqf.values.sum()    
    return(Q_per_aqf)
    
def get_Q_per_depth(df, depth_col, max_depth, *Q_cols, func=sum_Q_per_group, scale=True):
    assign_bin_depths(df, depth_col, max_depth)
    
    Q_per_depth = func(df, "depth (m)", *Q_cols)
    
    mid_depth = bin_to_mids(Q_per_depth.index)

    #Scale
    if scale==True:
        Q_per_depth /= Q_per_depth.values.sum()
    
    #For some reason I had to reindex it twice befire the actual index was assigned
    Q_per_depth = Q_per_depth.reindex(mid_depth)
    Q_per_depth = Q_per_depth.reindex(mid_depth) 
    return(Q_per_depth)

def hbar_plot_df(df, xlabel, path_out):
    df = df.sort_index(ascending=False)
    ax = df.plot(kind="barh", stacked=True)
    ax.set_xlabel(xlabel)
    
    plt.tight_layout()
    plt.savefig(path_out, dpi=300)
    plt.close()    

#%%Path management
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Nile\Nile_wells.csv"
path_Mekong = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mekong\Wells_Mekong.csv"
path_Mississippi = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mississippi\miss_wells.xlsx"
figure_folder   = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions"

#%%Read
abs_Nile = pd.read_csv(path_Nile)
abs_Mekong = pd.read_csv(path_Mekong)
N_Missi= pd.read_excel(path_Mississippi)

#%%Prepare data Nile
abs_Nile["Zmid"] = abs_Nile["Zmid"] * -1
abs_Nile["Q"] = abs_Nile["Q"] * -1

#%%Prepare data Mississippi
miss_index = pd.IntervalIndex.from_arrays(left  = N_Missi["LowDepth"], 
                                          right = N_Missi["HighDepth"])

N_Missi["z"] = miss_index.mid
N_Missi = N_Missi.loc[N_Missi["z"] < 500.]

#%%Process data for analysis
args_Nile = abs_Nile, "Zmid", 261, "Q"
Q_sum_depth_Nile     = get_Q_per_depth(*args_Nile, func=sum_Q_per_group)
Q_count_depth_Nile   = get_Q_per_depth(*args_Nile, func=count_Q_per_group)
#Q_sum_aqf_Nile       = get_Q_per_aqf(abs_Nile, *args_Nile[3:], func=sum_Q_per_group)

Nile_Q = get_Q_per_depth(*args_Nile, func=sum_Q_per_group, scale=False)
print(Nile_Q.sum()/497038.69440162915)

args_Mekong = abs_Mekong, "z", 641, "Q"
Q_sum_depth_Mekong   = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group)
Q_count_depth_Mekong = get_Q_per_depth(*args_Mekong, func=count_Q_per_group)
Q_sum_aqf_Mekong     = get_Q_per_aqf(abs_Mekong, *args_Mekong[3:], func=sum_Q_per_group)

Mekong_Q   = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group, scale=False)
print(Mekong_Q.sum()/879522.3821543201)

args_Mississippi = N_Missi, "z",501, "Count" #Count per 20m bin by summing the count
Q_count_Mississippi = get_Q_per_depth(*args_Mississippi, func=sum_Q_per_group)

#%%Testing
a = abs_Mekong.loc[abs_Mekong["depth (m)"] == abs_Mekong["depth (m)"].cat.categories[6]]

#%%Plot
sns.set()

hbar_plot_df(Q_sum_depth_Nile, "rel GW abstracted (-)", os.path.join(figure_folder, "Nile_sum.png"))
hbar_plot_df(Q_count_depth_Nile, "rel number of wells (-)", os.path.join(figure_folder, "Nile_count.png"))
hbar_plot_df(Q_sum_depth_Mekong, "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum.png"))
hbar_plot_df(Q_count_depth_Mekong, "rel number of wells (-)", os.path.join(figure_folder, "Mekong_count.png"))
hbar_plot_df(Q_sum_aqf_Mekong,  "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum_aqf.png"))
hbar_plot_df(Q_count_Mississippi, "rel number of wells (-)", os.path.join(figure_folder, "Mississippi_count.png"))