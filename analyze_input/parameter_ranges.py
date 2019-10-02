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

def get_df_plot(df, var):
    df_min=df[["Delta", "%s_min"%var]].rename(columns={"%s_min"%var : "%s"%var}).assign(stat="min")
    df_max=df[["Delta", "%s_max"%var]].rename(columns={"%s_max"%var : "%s"%var}).assign(stat="max")
    
    df_plot = pd.concat([df_min, df_max], axis=0)
    df_plot["log%s"%var] = np.log10(df_plot["%s"%var])
    
    df_plot["n_pubs"] = df_plot.groupby(by=["Delta"])["%s"%var].transform('count')/2
    df_plot["n_pubs"] = df_plot["n_pubs"].astype(np.int64)
    
    df_plot["n_pubs"] = df_plot["n_pubs"].replace(0, np.nan)
    df_plot = df_plot.sort_values(by="Delta")    
    return(df_plot)

def _get_logK(df, stat):
    logK = np.log10(df_aqf.loc[df_aqf["stat"]==stat]["Kaqf"])
    logK = logK.rename(column={"Kaqf" : "logKaqf"})
    return(logK)

#%%Path management

df_path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\Reftable.xlsx"

#%%Handle data

df = pd.read_excel(df_path, sheet_name="Hydrogeology_Raw", skiprows=[1])

to_keep = df["Aqt_explicit"] == 1
df['Kaqt_min'] = df['Kaqt_min'].where(to_keep, df['Kaqf_min']/df['Anisotropy_max'])
df['Kaqt_max'] = df['Kaqt_max'].where(to_keep, df['Kaqf_max']/df['Anisotropy_min'])

df['Anisotropy_min'] = df['Anisotropy_min'].where(to_keep, np.nan)
df['Anisotropy_max'] = df['Anisotropy_max'].where(to_keep, np.nan)

df_aqf = get_df_plot(df, "Kaqf")
df_aqt = get_df_plot(df, "Kaqt")
df_rch = get_df_plot(df, "Recharge")
df_ani = get_df_plot(df, "Anisotropy")

#%%Plot
sns.set(style="darkgrid")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize = (8, 12), sharey=True)

pointplotrange(y="Delta", x="logKaqf", data=df_aqf, ci="range", hue="n_pubs", join=False, errwidth=1.5, ax=ax1)
pointplotrange(y="Delta", x="logKaqt", data=df_aqt, ci="range", hue="n_pubs", join=False, errwidth=1.5, ax=ax2)
pointplotrange(y="Delta", x="Recharge", data=df_rch, ci="range", hue="n_pubs", join=False, errwidth=1.5, ax=ax3)
pointplotrange(y="Delta", x="logAnisotropy", data=df_ani, ci="range", hue="n_pubs", join=False, errwidth=1.5, ax=ax4)

ax1.get_legend().remove()
ax3.get_legend().remove()
ax4.get_legend().remove()

ax1.set(xlabel="$\log(K_{aqf} \; [m/d])$")
ax2.set(xlabel="$\log(K_{aqt} \; [m/d])$", ylabel="")
ax3.set(xlabel="$R \; [m/d]$")
ax4.set(xlabel="$\log(K_h/K_v \; [-])$", ylabel="")

plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "logK.png"), dpi=300)
plt.close()
#%%cause logKrange?
#Not enough data/too much classes right now...
df_dig = df.loc[:, ["Source K", "Cause variability K"]]
df_dig["logKaqfrange"] = (np.log10(df["Kaqf_max"])-np.log10(df["Kaqf_min"])).replace(0., np.nan)

g1 = sns.PairGrid(df_dig, x_vars=["logKaqfrange"],
                 y_vars=["Source K", "Cause variability K"],
                 height=5, aspect=.5)

g1.map(sns.stripplot, alpha="0.6", size=9, jitter=True)
g1.fig.set_size_inches(4,6)
plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "logKrange.png"), dpi=300)
plt.close()
#%%
fig, ax1 = plt.subplots(nrows=1, figsize = (8, 6))

df_aqf = df_aqf.dropna(subset=["Kaqf"])

sns.distplot(_get_logK(df_aqf, "min"), hist=True, rug=False, kde=False, ax=ax1)
sns.distplot(_get_logK(df_aqf, "max"), hist=True, rug=False, kde=False, ax=ax1)

plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "histlogKaqf"), dpi=300)
plt.close()