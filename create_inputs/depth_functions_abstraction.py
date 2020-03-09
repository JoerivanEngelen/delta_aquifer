# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:55:29 2020

@author: engelen
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
import os
import imod

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
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Nile_Delta\Data\Well_Data\Extractions_All.xlsx"
path_Mekong_folder = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mekong\Voor_joeri\VERSION_1"
figure_folder   = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions"

M_ipfs = glob(os.path.join(path_Mekong_folder, "WEL", "*.IPF"))
M_bots = glob(os.path.join(path_Mekong_folder, "BOT", "*.IDF"))
M_tops = glob(os.path.join(path_Mekong_folder, "TOP", "*.IDF"))

#%%Read
abs_Nile = pd.read_excel(path_Nile)

#%%Prepare data Nile
abs_Nile["Screen_Mid"] = (abs_Nile["Screen_Top"] + abs_Nile["Screen_Bottom"])/2

abs_Nile = abs_Nile[abs_Nile["Screen_Mid"] < 1000.]
abs_Nile = abs_Nile.dropna(subset=["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"])

#%%Prepare data Mekong
import xarray as xr

d = {"aqf_name" : ['QH', 'QP3', 'QP2-3', 'QP23', 'QP1', 'N22', 'N21', 'N13'],
     "layer"  : [3, 5, 7, 7, 9, 11, 13, 15],
     "aqf_nr" : [1, 2, 3, 3, 4, 5, 6, 7]}

lookup_df = pd.DataFrame(d).set_index("aqf_name")

#Filter names
filter_ipfs=True

if filter_ipfs:
    M_ipfs = [ipf for ipf in M_ipfs if "2015" not in ipf]
    M_ipfs = [ipf for ipf in M_ipfs if "_HCMC" not in ipf]

names = [os.path.splitext(os.path.basename(ipf))[0] for ipf in M_ipfs]
abs_Mekong = [imod.ipf.read(ipf) for ipf in M_ipfs]
abs_Mekong = pd.concat([df.assign(name=names[i]) for i, df in enumerate(abs_Mekong)])
abs_Mekong = abs_Mekong.reset_index(drop = True)
if not filter_ipfs:
    abs_Mekong = abs_Mekong.drop(columns = ["Aquifer", "Year"])
abs_Mekong["aquifer"] = abs_Mekong["name"].str.replace("2015", "")
abs_Mekong["aquifer"] = abs_Mekong["aquifer"].str.replace("_HCMC", "")
abs_Mekong["layer"] = abs_Mekong["aquifer"].replace(lookup_df.to_dict()["layer"])
abs_Mekong["aqf_nr"] = abs_Mekong["aquifer"].replace(lookup_df.to_dict()["aqf_nr"])
abs_Mekong["Q"] *= -1

tops_Mekong = xr.concat([imod.idf.open(top) for top in M_tops], dim="layer").sortby("layer")
bots_Mekong = xr.concat([imod.idf.open(bot) for bot in M_bots], dim="layer").sortby("layer")

bots_Mekong = bots_Mekong.assign_coords(layer = bots_Mekong.layer-1)

z = (tops_Mekong + bots_Mekong)/2

kwargs = dict(x     = xr.DataArray(abs_Mekong["X"],     dims="index", coords = {"index" : abs_Mekong.index }),
              y     = xr.DataArray(abs_Mekong["Y"],     dims="index"),
              layer = xr.DataArray(abs_Mekong["layer"], dims="index"), 
              method="nearest")

z_select   = z.sel(**kwargs)
bot_select = bots_Mekong.sel(**kwargs)
top_select = tops_Mekong.sel(**kwargs)

abs_Mekong["z"] = z_select.values *-1
abs_Mekong["top"] = top_select.values
abs_Mekong["bot"] = bot_select.values

#%%Process data for analysis
args_Nile = abs_Nile, "Screen_Mid", 321, "Irrigation", "Drinking", "Industrial", "Conjunctive_Use"
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

#%%Testing
a = abs_Mekong.loc[abs_Mekong["depth (m)"] == abs_Mekong["depth (m)"].cat.categories[6]]

#%%Plot
sns.set()

hbar_plot_df(Q_sum_depth_Nile, "rel GW abstracted (-)", os.path.join(figure_folder, "Nile_sum.png"))
hbar_plot_df(Q_count_depth_Nile, "rel number of wells (-)", os.path.join(figure_folder, "Nile_count.png"))
hbar_plot_df(Q_sum_depth_Mekong, "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum.png"))
hbar_plot_df(Q_count_depth_Mekong, "rel number of wells (-)", os.path.join(figure_folder, "Mekong_count.png"))
hbar_plot_df(Q_sum_aqf_Mekong,  "rel GW abstracted (-)", os.path.join(figure_folder, "Mekong_sum_aqf.png"))