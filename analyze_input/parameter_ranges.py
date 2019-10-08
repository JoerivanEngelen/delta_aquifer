# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:03:44 2019

@author: engelen
"""

import seaborn as sns
import pandas as pd
import numpy as np
from analyze_input.ci_range import pointplotrange
import os
import matplotlib.pyplot as plt
from functools import reduce
from collections import OrderedDict

def get_df_plot(df, var, i):
    df_min=df[["Delta", "%s_min"%var]].rename(columns={"%s_min"%var : "%s"%var}).assign(stat="min")
    df_max=df[["Delta", "%s_max"%var]].rename(columns={"%s_max"%var : "%s"%var}).assign(stat="max")
    
    df_plot = pd.concat([df_min, df_max], axis=0)
    df_plot["log%s"%var] = np.log10(df_plot["%s"%var])
    
    df_plot["n_pubs_%s"%i] = df_plot.groupby(by=["Delta"])["%s"%var].transform('count')/2
    df_plot["n_pubs_%s"%i] = df_plot["n_pubs_%s"%i].astype(np.int64)
    
    df_plot["n_pubs_%s"%i] = df_plot["n_pubs_%s"%i].replace(0, np.nan)
    df_plot = df_plot.sort_values(by="Delta")    
    return(df_plot)

#%%Path management
df_path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\Reftable.xlsx"

#%%Figsizes
agu_half  = (19/2.54, 11.5/2.54)
agu_whole = (19/2.54, 23/2.54)
agu_half_vert = (9.5/2.54, 23/2.54)

#%%Handle data
df = pd.read_excel(df_path, sheet_name="Hydrogeology_Raw", skiprows=[1])

to_keep = df["Aqt_explicit"] == 1
df['Kaqt_min'] = df['Kaqt_min'].where(to_keep, df['Kaqf_min']/df['Anisotropy_max'])
df['Kaqt_max'] = df['Kaqt_max'].where(to_keep, df['Kaqf_max']/df['Anisotropy_min'])

df['Anisotropy_min'] = df['Anisotropy_min'].where(to_keep, np.nan)
df['Anisotropy_max'] = df['Anisotropy_max'].where(to_keep, np.nan)

hdrglgy = ["Kaqf", "Kaqt", "Recharge", "Anisotropy"]
df_hdrglgy_ls = [get_df_plot(df, var, i) for i, var in enumerate(hdrglgy)]
df_hdrglgy = reduce(lambda x, y: pd.merge(x, y, on = 'Delta'), df_hdrglgy_ls)

sheets = ["BC_Raw", "Lithology_Raw", "Geometry_Raw"]

df_ls = list(pd.read_excel(df_path, sheet_name=sheets, skiprows=[1]).values())
df_all = reduce(lambda x, y: pd.merge(x, y, on = 'Delta'), df_ls)

#%%Plot hydrogeology ranges per delta
sns.set(style="darkgrid")
sns.set_context("paper")

fig, ax_grid = plt.subplots(nrows=1, ncols=4, figsize = agu_half, sharey=True)
ax_lst = ax_grid.flatten()

opts = dict(ci="range", join=False, errwidth=1.5,scale=0.4)

var2plot = OrderedDict(
            logKaqf = "$\log(K_{aqf} \; [m/d])$",
            logKaqt = "$\log(K_{aqt} \; [m/d])$", 
            Recharge = "$R \; [m/d]$", 
            logAnisotropy = "$\log(K_h/K_v \; [-])$"
            )

for i, var in enumerate(var2plot.keys()):
    pointplotrange(y="Delta", x=var, data=df_hdrglgy, ax=ax_lst[i], hue="n_pubs_%s"%i, **opts)
    labels = dict(xlabel=var2plot[var])
    if i != 0:
        labels["ylabel"] = ""
    ax_lst[i].set(**labels)

for i in [0, 1, 2]:
    ax_lst[i].get_legend().remove()

plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "logK.png"), dpi=300)
plt.close()


#%%Plot histogram all
fig, ax_grid = plt.subplots(nrows=3, ncols=4, figsize = agu_half, sharey=False)
ax_lst = ax_grid.flatten()

var2plot = ['l_a', 'beta', r'H_a/H_b', 'H_b', r'Mud/Total', 
            'l_conf', 'N_aqt', 'N_pal', 'l_tra', 'N_chan']

opts = dict(hist=True, rug=False, kde=False, bins='doane')

for i, var in enumerate(var2plot):
    sns.distplot(df_all[var].dropna(), ax=ax_lst[i], **opts)

ax_grid[-1, -1].axis('off')
ax_grid[-1, -2].axis('off')

plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "hist_all.png"), dpi=300)
plt.close()