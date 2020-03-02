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

def sum_Q_per_depth(df, *Q_cols):
    return(df.groupby("depth (m)").sum()[list(Q_cols)])

def count_Q_per_depth(df, *Q_cols):
    return(df.groupby("depth (m)").count()[list(Q_cols)])

def bin_to_mids(index):
    return([(a.left + a.right)/2 for a in index])

def get_Q_per_depth(df, depth_col, max_depth, *Q_cols, method=sum_Q_per_depth):
    assign_bin_depths(df, depth_col, max_depth)
    
    Q_sum_depth = sum_Q_per_depth(df, *Q_cols).sort_index(ascending=False)
    
    mid_depth = bin_to_mids(Q_sum_depth.index)

    #Scale
    Q_sum_depth /= Q_sum_depth.values.sum()
    
    #For some reason I had to reindex it twice befire the actual index was assigned
    Q_sum_depth = Q_sum_depth.reindex(mid_depth)
    Q_sum_depth = Q_sum_depth.reindex(mid_depth) 
    return(Q_sum_depth)

#%%Path management
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Nile_Delta\Data\Well_Data\Extractions_All.xlsx"
path_Mekong_folder = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mekong\Voor_joeri\VERSION_1"
path_out_Nile   = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions\Nile.png"
path_out_Mekong = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions\Mekong.png"

M_ipfs = glob(os.path.join(path_Mekong_folder, "WEL", "*.IPF"))
M_bots = glob(os.path.join(path_Mekong_folder, "BOT", "*.IDF"))
M_tops = glob(os.path.join(path_Mekong_folder, "TOP", "*.IDF"))

#%%Read
abs_Nile = pd.read_excel(path_Nile)

#%%Prepare data Nile
abs_Nile["Screen_Mid"] = (abs_Nile["Screen_Top"] + abs_Nile["Screen_Bottom"])/2

abs_Nile = abs_Nile[abs_Nile["Screen_Mid"] < 1000.]
abs_Nile = abs_Nile.dropna(subset=["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"])

args_Nile = abs_Nile, "Screen_Mid", 321, "Irrigation", "Drinking", "Industrial", "Conjunctive_Use"
Q_sum_depth_Nile = get_Q_per_depth(*args_Nile, method=sum_Q_per_depth,)

#%%Prepare data Mekong
import xarray as xr

lookup = {"QH" : 3,
          "QP3" : 5,
          "QP2-3" : 7,
          "QP23" : 7, 
          "QP1" : 9,
          "N22" : 11,
          "N21" : 13,
          "N13" : 15}

names = [os.path.splitext(os.path.basename(ipf))[0] for ipf in M_ipfs]
abs_Mekong = [imod.ipf.read(ipf) for ipf in M_ipfs]
abs_Mekong = pd.concat([df.assign(name=names[i]) for i, df in enumerate(abs_Mekong)])
abs_Mekong = abs_Mekong.reset_index(drop = True)
abs_Mekong = abs_Mekong.drop(columns = ["Aquifer", "Year"])
abs_Mekong["aquifer"] = abs_Mekong["name"].str.replace("2015", "")
abs_Mekong["aquifer"] = abs_Mekong["aquifer"].str.replace("_HCMC", "")
abs_Mekong["layer"] = abs_Mekong["aquifer"].replace(lookup)
abs_Mekong["Q"] *= -1

tops_Mekong = xr.concat([imod.idf.open(top) for top in M_tops], dim="layer").sortby("layer")
bots_Mekong = xr.concat([imod.idf.open(bot) for bot in M_bots], dim="layer").sortby("layer")

bots_Mekong = bots_Mekong.assign_coords(layer = bots_Mekong.layer-1)

z = (tops_Mekong + bots_Mekong)/2

z_select = z.sel(x=xr.DataArray(abs_Mekong["X"], coords = {"index" : abs_Mekong.index }, dims="index"),
      y=xr.DataArray(abs_Mekong["Y"], dims="index"),
      layer=xr.DataArray(abs_Mekong["layer"], dims="index"), method="nearest")

abs_Mekong["z"] = z_select.values *-1

args_Mekong = abs_Mekong, "z", 641, "Q"
Q_sum_depth_Mekong = get_Q_per_depth(*args_Mekong, method=sum_Q_per_depth)

#%%Plot
sns.set()
ax = Q_sum_depth_Nile.plot(kind="barh", stacked=True)
ax.set_xlabel("rel GW abstracted (-)")

plt.tight_layout()
plt.savefig(path_out_Nile, dpi=300)
plt.close()

ax = Q_sum_depth_Mekong.plot(kind="barh", stacked=True)
ax.set_xlabel("rel GW abstracted (-)")

plt.tight_layout()
plt.savefig(path_out_Mekong, dpi=300)
plt.close()
