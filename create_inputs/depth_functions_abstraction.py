# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:55:29 2020

@author: engelen
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

from glob import glob
import os
import imod


#%%Path management
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Nile_Delta\Data\Well_Data\Extractions_All.xlsx"
path_Mekong_folder = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\wells\Mekong\Voor_joeri\VERSION_1"
path_out = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions\Nile.png"
path_out_Mekong = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions\Mekong.png"

M_ipfs = glob(os.path.join(path_Mekong_folder, "WEL", "*.IPF"))
M_bots = glob(os.path.join(path_Mekong_folder, "BOT", "*.IDF"))
M_tops = glob(os.path.join(path_Mekong_folder, "TOP", "*.IDF"))

#%%Read
abs_Nile = pd.read_excel(path_Nile)

#%%Process data Nile
abs_Nile["Screen_Mid"] = (abs_Nile["Screen_Top"] + abs_Nile["Screen_Bottom"])/2

abs_Nile = abs_Nile[abs_Nile["Screen_Mid"] < 1000.]
abs_Nile = abs_Nile.dropna(subset=["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"])

abs_Nile["depth (m)"] = pd.cut(abs_Nile["Screen_Mid"], np.arange(0, 321, 20), duplicates="drop")

Q_sum_depth = abs_Nile.groupby("depth (m)").sum()[["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"]].reset_index()
Q_sum_depth = Q_sum_depth.set_index("depth (m)").sort_index(ascending=False)

mid_depth = [(a.left + a.right)/2 for a in Q_sum_depth.index]

#Scale
Q_sum_depth /= Q_sum_depth.values.sum()

Q_sum_depth = Q_sum_depth.reindex(mid_depth)
Q_sum_depth = Q_sum_depth.reindex(mid_depth)

#%%Process data Mekong
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

tops_Mekong = xr.concat([imod.idf.open(top) for top in M_tops], dim="layer").sortby("layer")
bots_Mekong = xr.concat([imod.idf.open(bot) for bot in M_bots], dim="layer").sortby("layer")

bots_Mekong = bots_Mekong.assign_coords(layer = bots_Mekong.layer-1)

z = (tops_Mekong + bots_Mekong)/2

z_select = z.sel(x=xr.DataArray(abs_Mekong["X"], coords = {"index" : abs_Mekong.index }, dims="index"),
      y=xr.DataArray(abs_Mekong["Y"], dims="index"),
      layer=xr.DataArray(abs_Mekong["layer"], dims="index"), method="nearest")

abs_Mekong["z"] = z_select.values
abs_Mekong["depth (m)"] = pd.cut(abs_Mekong["z"]*-1, np.arange(0, 641, 20), duplicates="drop")

Q_sum_depth_Mekong = abs_Mekong.groupby("depth (m)").sum()["Q"] * -1

mid_depth_Mekong = [(a.left + a.right)/2 for a in Q_sum_depth_Mekong.index]

Q_sum_depth_Mekong /= Q_sum_depth_Mekong.values.sum()

Q_sum_depth_Mekong = Q_sum_depth_Mekong.reindex(mid_depth_Mekong)
Q_sum_depth_Mekong = Q_sum_depth_Mekong.reindex(mid_depth_Mekong).sort_index(ascending=False)

#%%Exponential distribution
x = np.arange(0, 321, 20)
scale = 1/0.30
y = stats.expon.pdf(x/15, scale=scale)

#%%Plot
sns.set()
ax = Q_sum_depth.plot(kind="barh", stacked=True)
plt.plot(y, np.arange(len(x))[::-1]-1, color="k")
ax.set_xlabel("rel GW abstracted (-)")

plt.tight_layout()
plt.savefig(path_out, dpi=300)

#ax = Q_sum_depth_Mekong.plot(kind="barh", stacked=True)
#ax.set_xlabel("rel GW abstracted (-)")
#
#plt.tight_layout()
#plt.savefig(path_out_Mekong, dpi=300)
#%%Save
