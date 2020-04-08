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

from itertools import product

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

def save_hbar_df(df, xlabel, path_out):
    ax = hbar_df(df)
    ax.set_xlabel(xlabel)
    
    plt.tight_layout()
    plt.savefig(path_out, dpi=300)
    plt.close()    

def hbar_df(df, ax=None):
    df = df.sort_index(ascending=False)
#    ax = df.plot(kind="barh", stacked=True, ax = ax, width=1, align="center")
    ax.barh(df.index, df.values, align="center", height=20)
    return(ax)

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
dfs_validate = [dfs.loc[dfs.index<260] for dfs in dfs_validate]
#dfs_validate = [df.set_index(df.index * -1) for df in dfs_validate]
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
abs_NL = abs_NL.rename(columns={"q" : "Q"})

#%%Process data for analysis
dfs_obs = {}

args_Nile = abs_Nile, "Zmid", 251, "Q"
dfs_obs["Nile"] = dict(sum = get_Q_per_depth(*args_Nile, func=sum_Q_per_group),
              count =get_Q_per_depth(*args_Nile, func=count_Q_per_group))

#Nile_Q = get_Q_per_depth(*args_Nile, func=sum_Q_per_group, scale=False)
#print(Nile_Q.sum()/1.8e6)

args_Mekong = abs_Mekong, "z", 301, "Q"

dfs_obs["Mekong"] = dict(sum = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group),
                count = get_Q_per_depth(*args_Mekong, func=count_Q_per_group),
                sum_aqf = get_Q_per_aqf(abs_Mekong, *args_Mekong[3:], func=sum_Q_per_group))

#Mekong_Q   = get_Q_per_depth(*args_Mekong, func=sum_Q_per_group, scale=False)
#print(Mekong_Q.sum()/2.2e6)

args_Mississippi = N_Missi, "z",301, "Count" #Count per 20m bin by summing the count
dfs_obs["Mississippi"] = dict(count = get_Q_per_depth(*args_Mississippi, func=sum_Q_per_group))

args_NL = abs_NL, "z", 301, "Q"
dfs_obs["Rhine-Meuse"] = dict(sum = get_Q_per_depth(*args_NL, func=sum_Q_per_group),
                        count = get_Q_per_depth(*args_NL, func=count_Q_per_group))

#%%Setup plot
agu_whole = (19/2.54, 23/2.54)

#%%Prepare figure and gridspec
sns.set()

fig = plt.figure(figsize=agu_whole)
gs_main = fig.add_gridspec(3,2)
gs_obs = [gs_main[i, 0].subgridspec(1,1) for i in range(3)] #Slice mapping not available on gridspec objects, therefore work with these lists
gs_mod = [gs_main[i, 1].subgridspec(3,3) for i in range(3)]

#%%Initiate subplots
ax_obs = []
ax_mod = [[], [], []]
ax_title = []

for r in range(len(gs_obs)):
    ax_obs.append(fig.add_subplot(gs_obs[r][0,0]))
    for i, j in product(range(3), repeat=2):
            ax_mod[r].append(fig.add_subplot(gs_mod[r][i, j], 
                  sharex = ax_obs[r], 
                  sharey = ax_obs[r]))

#%%Plot
for r, delta in enumerate(delta_names):
    hbar_df(dfs_obs[delta]["sum"]["Q"], ax=ax_obs[r])
    ax_obs[r].invert_yaxis()
    for i, (ax, idx) in enumerate(zip(ax_mod[r], idxs[r])):
        hbar_df(dfs_validate[r]["sum{}".format(idx)], ax=ax)
        ax.label_outer()

#%%Fix ticklabels and axis titles
row_title_idx = range(2, 9, 3)
col_title_idx = range(3)

row_titles = ["$K_{h}$ min", "$K_{h}$ max", "$K_{h}$ mean"]
col_titles = ["$K_{v}$ max", "$K_{v}$ min", "$K_{v}$ mean"]

for r in range(len(gs_obs)):
    if r == 0:
        ax_obs[r].set_title("Observed", weight = "bold")
    ax_obs[r].set_ylabel("depth (m)")
    ax_obs[r].set_xlabel("$Q$/$Q_{tot}$ (-)")
    ax_obs[r].annotate(delta_names[r], xy=(-0.25, .5), xycoords="axes fraction", 
          rotation=90, ha="right", va="center", weight = "bold")
    
    for id, title in zip(row_title_idx, row_titles):
        ax_mod[r][id].annotate(title, xy=(1.02, .5), xycoords="axes fraction",
            rotation=270, ha="left", va="center")
    for id, title in zip(col_title_idx, col_titles):
        ax_mod[r][id].set_title(title)

plt.tight_layout()
fig.canvas.draw() #Ensure ticklabels are created so we can only show every second ticklabel
for r, delta in enumerate(delta_names):
    xlabels = [txt_obj._text.rstrip("0") if ((i % 2) == 0) else "" for i, txt_obj in enumerate(ax_obs[r].get_xticklabels())]
    ax_obs[r].set_xticklabels(xlabels)
    ylabels = [txt_obj._text if ((i % 2) == 0) else "" for i, txt_obj in enumerate(ax_obs[r].get_yticklabels())]
    ax_obs[r].set_yticklabels(ylabels)

#%%Save
plt.savefig(os.path.join(figure_folder, "Comparison_Validation.png"), dpi=300)