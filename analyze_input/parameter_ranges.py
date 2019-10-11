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
import matplotlib.gridspec as gridspec
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

hdrglgy = ["Kaqf", "Kaqt", "Recharge", "Anisotropy", "Ss"]
df_hdrglgy_ls = [get_df_plot(df, var, i) for i, var in enumerate(hdrglgy)]
#TODO: 'Delta' has duplicate values, so marge creates huge Dataframe with lots of NaNs
#This does not influence the figures, but is confusing to the naive user.
df_hdrglgy = reduce(lambda x, y: pd.merge(x, y, on = 'Delta'), df_hdrglgy_ls)

sheets = ["BC_Raw", "Lithology_Raw", "Geometry_Raw"]

df_ls = list(pd.read_excel(df_path, sheet_name=sheets, skiprows=[1]).values())
df_all = reduce(lambda x, y: pd.merge(x, y, on = 'Delta'), df_ls)

#Filter out disproportionally large L_b values
df_all['l_a'].where((df_all['L_b'] < 500.) , other=np.nan, inplace=True)
df_all['l_a'].where((df_all['Delta'] != "Nakdong"), other=np.nan, inplace=True)

#%%Select variables to plot

#Ranges of hydrogeological parameters
hdrglgy2plot = OrderedDict((
            ("logKaqf",  "$\log(K_{h,aqf} \; [m/d])$"),
            ("logKaqt",  "$\log(K_{v,aqt} \; [m/d])$"), 
            ("Recharge", "$R \; [m/d]$"), 
            ("logAnisotropy", "$\log(K_h/K_v \; [-])$"),
            ("logSs", "$\log(S_s)$")
            ))

#Histograms for other parameters
var2plot = OrderedDict((("l_a",         "$l_a$"), 
                        ("beta",        r"$\beta$"),
                        (r"H_a/H_b",    "$f_H$"),
                        ('H_b',         "$H_b$"),
                        (r'Mud/Total',  "$f_{aqt}$"),
                        ('l_conf',      "$l_{conf}$"),
                        ('N_aqt',       "$N_{aqt}$"),
                        ('N_pal',       "$N_{pal}$"),
                        ('l_tra',       "$l_{tra}$"),
                        ('N_chan',      "$N_{chan}$")
                        ))

#%%Plot hydrogeology ranges per delta and histograms for the rest
#Set defaults
sns.set(style="darkgrid")
sns.set_context("paper")

n_range = len(hdrglgy2plot.keys())
n_hist = len(var2plot.keys())
n_hist_cols = np.ceil(n_hist/2).astype(int)

#Create figure and gridspec
fig = plt.figure(figsize = agu_whole)

gs0 = gridspec.GridSpec(2, 1, figure=fig)
gs00 = gridspec.GridSpecFromSubplotSpec(1, n_range+1, subplot_spec=gs0[0])
gs01 = gridspec.GridSpecFromSubplotSpec(2, n_hist_cols, subplot_spec=gs0[1], hspace=0.3, wspace=0.3)

hdrglgy_ax_lst = [fig.add_subplot(gs00[1])]
hdrglgy_ax_lst += [fig.add_subplot(gs00[i], sharey=hdrglgy_ax_lst[0]) for i in range(2,n_range+1)]

hist_ax_lst = [fig.add_subplot(gs01[i]) for i in range(n_hist)]

##Plot input ranges for hydrogeology
opts = dict(ci="range", join=False, errwidth=1.5,scale=0.4)

for i, var in enumerate(hdrglgy2plot.keys()):
    pointplotrange(y="Delta", x=var, data=df_hdrglgy, ax=hdrglgy_ax_lst[i], hue="n_pubs_%s"%i, **opts)
    labels = dict(xlabel=hdrglgy2plot[var])
    if i != 0:
        labels["ylabel"] = ""
    hdrglgy_ax_lst[i].set(**labels)

for i in range(n_range-1):
    #Get rid of automatically added legends
    hdrglgy_ax_lst[i].get_legend().remove()
    #We have to remove ticklabels manually, as gridspec does not do hide this
    plt.setp(hdrglgy_ax_lst[i+1].get_yticklabels(), visible=False)

##Plot histograms for the other parameters
opts = dict(hist=True, rug=False, kde=False, bins='doane')

for i, var in enumerate(var2plot.keys()):
    sns.distplot(df_all[var].dropna(), ax=hist_ax_lst[i], axlabel=var2plot[var], **opts)

plt.tight_layout()
plt.savefig(os.path.join(df_path, "..", "input_distributions.png"), dpi=300)
plt.savefig(os.path.join(df_path, "..", "input_distributions.pdf"))
plt.close()