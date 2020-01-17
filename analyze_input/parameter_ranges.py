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
from pkg_resources import resource_filename

def get_df_plot(df, var, i):
    df_min=df[["Delta", "%s_min"%var]].rename(columns={"%s_min"%var : "%s"%var}).assign(stat="min")
    df_max=df[["Delta", "%s_max"%var]].rename(columns={"%s_max"%var : "%s"%var}).assign(stat="max")
    
    df_min["n_pubs_%s"%i] = df_min.groupby(by=["Delta"])["%s"%var].transform('count').astype("Int64").replace(0, np.nan) #Int64 supports NaNs
    df_max["n_pubs_%s"%i] = df_min["n_pubs_%s"%i]

    df_min=df_min.groupby(by=["Delta"]).min()
    df_max=df_max.groupby(by=["Delta"]).max()
    
    df_plot = pd.concat([df_min, df_max], axis=0)
    df_plot["log%s"%var] = np.log10(df_plot["%s"%var])
    
    df_plot = df_plot.sort_values(by="Delta")    
    return(df_plot)

def lin(x, slope, intercept, clip):
    return(np.clip(x*slope+intercept, clip, None))

#%%Path management
datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))
df_path  = os.path.join(datafol, "Reftable.xlsx")
out_path = os.path.join(datafol, "..", "example", "scratch_figures")

#%%Figsizes
agu_small = (9.5/2.54, 11.5/2.54)
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
df_hdrglgy = reduce(lambda x, y: pd.merge(x, y, on = 'Delta'), df_hdrglgy_ls).reset_index(level=0)

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
#            ("logSs", "$\log(S_s)$")
            ))

#Histograms for other parameters
var2plot = OrderedDict((("l_a",         "$l_a$ [-]"), 
                        ("beta",        r"$\beta$ [rad]"),
                        (r"H_a/H_b",    "$f_H$ [-]"),
                        ('H_b',         "$H_b$ [m]"),
                        (r'Mud/Total',  "$f_{aqt}$ [-]"),
                        ('l_conf',      "$l_{conf}$ [-]"),
                        ('N_aqt',       "$N_{aqt}$ [-]"),
                        ('N_pal',       "$N_{pal}$ [-]"),
                        ('l_tra',       "$l_{tra}$ [-]"),
                        ('N_chan',      "$N_{chan}$ [-]")
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

leg = hdrglgy_ax_lst[n_range-1].get_legend()
leg.set_title("$N_{pubs}$")

##Plot histograms for the other parameters
opts = dict(hist=True, rug=False, kde=False, bins='doane')

for i, var in enumerate(var2plot.keys()):
    sns.distplot(df_all[var].dropna(), ax=hist_ax_lst[i], axlabel=var2plot[var], **opts)

plt.tight_layout()
plt.savefig(os.path.join(out_path, "input_distributions.png"), dpi=300)
plt.savefig(os.path.join(out_path, "input_distributions.svg"))
plt.savefig(os.path.join(out_path, "input_distributions.pdf"))
plt.close()

#%%Explore relation H_b and N_aqt
#fig = plt.figure(figsize = agu_small)
#gs0 = gridspec.GridSpec(1, 1, figure=fig)
#ax = fig.add_subplot(gs0[0])
#
#df_all["logH_b"] = np.log10(df_all["H_b"])
#sns.regplot(x="logH_b", y="N_aqt", data=df_all, robust=True, ax=ax)
#
##Seaborn does not return coefficients to 
#mean = ax.lines[0]
#verts = ax.collections[1]._paths[0]._vertices
##ASSUME: Second half of verts is upper line. Not sure if this is guarenteed.
#highest_ci = verts[int(verts.shape[0]/2):, :]
#
#plt.plot(highest_ci[:, 0], highest_ci[:, 1])
##plt.plot(highest_ci[:, 0], lin(highest_ci[:, 0], 8/0.95, -15.26, 2.))
#
#plt.tight_layout()
#plt.savefig(os.path.join(df_path, "..", "n_aqt_vs_logH_b.png"), dpi=300)
#plt.savefig(os.path.join(df_path, "..", "n_aqt_vs_logH_b.pdf"))
#plt.close()