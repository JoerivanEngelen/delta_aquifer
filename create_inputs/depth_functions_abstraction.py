# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:55:29 2020

@author: engelen
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from pkg_resources import resource_filename
from glob import glob
import os
import json

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
    Q_per_depth = reindex_with_mids(Q_per_depth, scale=scale)
    return(Q_per_depth)

def reindex_with_mids(df, scale=True):
    mid_depth = bin_to_mids(df.index)
    
    #Scale
    if scale==True:
        df /= df.values.sum()
    
    #For some reason I had to reindex it twice before the actual index was assigned
    df = df.reindex(mid_depth)
    df = df.reindex(mid_depth) 
    return(df)

def hbar_plot_df(df, xlabel, path_out):
    df = df.sort_index(ascending=False)
    ax = df.plot(kind="barh", stacked=True)
    ax.set_xlabel(xlabel)
    
    plt.tight_layout()
    plt.savefig(path_out, dpi=300)
    plt.close()    

def reindex_str_to_IntervalIndex(df):
    """Saving Intervals with pandas is horrible at the moment as there are no parsers yet, so this
    clunky code is needed unfortunately.
    """
    idx = pd.IntervalIndex([pd.Interval(*json.loads(i)) for i in df.index.str.replace("(", "[").to_list()])
    df.index = idx
    return(df)

#%%Path management
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Nile\Nile_wells.csv"
path_Mekong = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mekong\Wells_Mekong.csv"
path_Mississippi = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mississippi\miss_wells.xlsx"
path_NL = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Rhine-Meuse\wel_NL.csv"

figure_folder   = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions"

globpath_validation = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results_30deltas\validation_deltas\*"
paths_validation = glob(globpath_validation)

delta_names = ["Nile", "Mekong", "Rhine-Meuse"]#, "Mississippi"]
datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)
idxs = [par.loc[par["Delta"]==delta].index for delta in delta_names]

paths_val = [[os.path.join(globpath_validation, "..", "AD_i{:03d}".format(i), "well_dist_depth.csv") for i in idx] for idx in idxs]

#%%Read
abs_Nile = pd.read_csv(path_Nile, index_col=0)
abs_Mekong = pd.read_csv(path_Mekong, index_col=0)
N_Missi= pd.read_excel(path_Mississippi)
abs_NL = pd.read_csv(path_NL, index_col=0)

#%%Prepare model depth distributions
dfs_validate = [[pd.read_csv(path).set_index("well_depth_bins").add_suffix(sim_nr) for sim_nr, path in zip(idxs[i], paths_val[i])] for i in range(len(idxs))]
dfs_validate = [pd.concat(df, axis=1) for df in dfs_validate]
dfs_validate = [reindex_str_to_IntervalIndex(df).sort_index() for df in dfs_validate]
dfs_validate = [reindex_with_mids(df, scale=False) for df in dfs_validate]
dfs_validate = [df.set_index(df.index * -1) for df in dfs_validate]
#%%Prepare data Nile
abs_Nile["Zmid"] = abs_Nile["Zmid"] * -1
abs_Nile["Q"] = abs_Nile["Q"] * -1

#%%Prepare data Mississippi
miss_index = pd.IntervalIndex.from_arrays(left  = N_Missi["LowDepth"], 
                                          right = N_Missi["HighDepth"])

N_Missi["z"] = miss_index.mid
N_Missi = N_Missi.loc[N_Missi["z"] < 500.]

#%%Prepare data NL
abs_NL["z"] = abs_NL["z"] * -1
abs_NL["q"] = abs_NL["q"] * -1

#%%Process data for analysis
args_Nile = abs_Nile, "Zmid", 261, "Q"
Q_sum_depth_Nile     = get_Q_per_depth(*args_Nile, func=sum_Q_per_group)
Q_count_depth_Nile   = get_Q_per_depth(*args_Nile, func=count_Q_per_group)
#Q_sum_aqf_Nile       = get_Q_per_aqf(abs_Nile, *args_Nile[3:], func=sum_Q_per_group)

Nile_Q = get_Q_per_depth(*args_Nile, func=sum_Q_per_group, scale=False)
#print(Nile_Q.sum()/497038.69440162915)
print(Nile_Q.sum()/1.8e6)
args_Mekong = abs_Mekong, "z", 641, "Q"
Q_sum_depth_Mekong   = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group)
Q_count_depth_Mekong = get_Q_per_depth(*args_Mekong, func=count_Q_per_group)
Q_sum_aqf_Mekong     = get_Q_per_aqf(abs_Mekong, *args_Mekong[3:], func=sum_Q_per_group)

Mekong_Q   = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group, scale=False)
print(Mekong_Q.sum()/2.2e6)

args_Mississippi = N_Missi, "z",501, "Count" #Count per 20m bin by summing the count
Q_count_Mississippi = get_Q_per_depth(*args_Mississippi, func=sum_Q_per_group)

args_NL = abs_NL, "z", 361, "q"
Q_sum_depth_NL   = get_Q_per_depth(*args_NL, func=sum_Q_per_group)
Q_count_depth_NL = get_Q_per_depth(*args_NL, func=count_Q_per_group)

##%%Testing
#a = abs_Mekong.loc[abs_Mekong["depth (m)"] == abs_Mekong["depth (m)"].cat.categories[6]]

#%%Setup plot
agu_whole = (19/2.54, 23/2.54)
#gridspec = 

#%%Plot
sns.set()

hbar_plot_df(Q_sum_depth_Nile, "rel GW abstracted (-)", os.path.join(figure_folder, "Nile_sum.png"))
hbar_plot_df(Q_count_depth_Nile, "rel number of wells (-)", os.path.join(figure_folder, "Nile_count.png"))
hbar_plot_df(Q_sum_depth_Mekong, "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum.png"))
hbar_plot_df(Q_count_depth_Mekong, "rel number of wells (-)", os.path.join(figure_folder, "Mekong_count.png"))
hbar_plot_df(Q_sum_aqf_Mekong,  "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum_aqf.png"))
hbar_plot_df(Q_count_Mississippi, "rel number of wells (-)", os.path.join(figure_folder, "Mississippi_count.png"))
hbar_plot_df(Q_sum_depth_NL,  "rel GW abstracted (-)", os.path.join(figure_folder, "NL_sum_depth.png"))
hbar_plot_df(Q_count_depth_NL, "rel number of wells (-)", os.path.join(figure_folder, "NL_count_depth.png"))