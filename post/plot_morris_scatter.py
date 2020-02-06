# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:12:46 2019

@author: engelen
"""

#Requires adjustText https://github.com/Phlya/adjustText

import os
import re

import pandas as pd
from glob import glob
import numpy as np
from SALib.analyze.morris import analyze, compute_elementary_effects
from pkg_resources import resource_filename
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from collections import defaultdict

from adjustText import adjust_text

from itertools import combinations, product, chain

#%%Functions to parse
def tree(): return defaultdict(tree)

def get_model_id(file):     
    pattern=r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return(number)

def str_to_mathtxt(s):
    s_ls = s.split("_")
    var = s_ls[0]
    #Create subscript
    if len(s_ls) > 1:
        sub = s_ls[1:]
        sub = ",".join(sub)
        sub = "_{%s}" % sub
    else:
        sub=""
    s = "${}{}$".format(var, sub)
    return(s)

def convert_texts(texts):
    tex = {"alpha" : r"\alpha",
           "beta"  : r"\beta", 
           "phi_f"   : r"\phi_f",
           "Kh_Kv" : r"K_h/K_v",
           "Kh_aqf" : "K_h_aqf",
           "Kv_aqt" : "K_v_aqt"}
    
    texts = [tex[t] if t in tex.keys() else t for t in texts]
    texts = [str_to_mathtxt(s) for s in texts]
    return(texts)

def labeled_scatter(fig, ax, data, labels, legend_out=True, 
                    plot_y_label=True, plot_x_label=True, 
                    **exceptions):
    markers={"decrease": "v", "increase": "^"}
    if legend_out == True:
        legend_out = "brief"
    
    sns.scatterplot(x = "mu_star", y = "sigma", style="sign", ax = ax, 
                    data = data, hue = "group", palette = "Set2", alpha=1.0,
                    markers = markers, legend=legend_out, s=80)
    
    x_max=np.nanmax(data["mu_star"])
    y_max=np.nanmax(data["sigma"])
    x_max += 0.1 * x_max
    y_max += 0.1 * y_max
    ax.set_xlim(left=0.0, right=x_max)
    ax.set_ylim(bottom=0.0, top=y_max)
    
    x_line = np.linspace(0, x_max)
    ax.plot(x_line, x_line, ls=":", alpha=0.5, color="k")
    
    if plot_x_label:
        ax.set_xlabel("$\mu *$")
    else:
        ax.axes.get_xaxis().get_label().set_visible(False)
    
    if plot_y_label:
        ax.set_ylabel("$\sigma$")
    else:
        ax.axes.get_yaxis().get_label().set_visible(False)
    
    dx=np.diff(np.array(ax.get_xlim()))
    dy=np.diff(np.array(ax.get_ylim()))
    
    exceptions_off = [key for key, value in exceptions.items() if value==0]
    exceptions_on = [key for key, value in exceptions.items() if value==1]
    
    texts=[]
    
    for x, y, label in zip(data["mu_star"], data["sigma"], labels):
        if (
                (x > (0.3  * dx)) and (label not in exceptions_off)
                ) or (
                (y > (0.41 * dy)) and (label not in exceptions_off)
                ) or (
                (label in exceptions_on)
            ):
            texts.append(ax.text(x, y, label, va="center", ha="center"))

    adjust_text(texts, ax=ax, expand_points=(1.7, 2), force_text=(1.0, 1.0),  
                arrowprops=dict(arrowstyle="wedge", color="k", alpha=0.2, lw=0.1))
    
    if legend_out != False:
        leg = ax.get_legend()
        labels = [t._text for t in leg.texts]
        fig.legend(leg.legendHandles, labels, loc='upper left', 
                   bbox_to_anchor=(0.75, 0.5, 0.25, 0.5), mode="expand",
                   edgecolor="#ffffff")
        ax.get_legend().remove()
    
    texts = [t._text for t in texts]
    
    return(texts)

#%%Path management
#path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results\outputs_of_interest"
path = r"g:\synthdelta\results\outputs_of_interest"
traj_id_path = os.path.abspath(resource_filename("delta_aquifer", "../data/traj_id.csv"))
outf = os.path.join(path, "..", "morris")

#%%Process paths
#Unfinished trajectories
exclude=list(range(96, 120)) + list(range(216, 240))

paths = glob(os.path.join(path, "*.csv"))
paths.sort()
idxs = [get_model_id(path) for path in paths]
idxs = [idx for idx in idxs if idx not in exclude]
#Check duplicates:
dupl = set([x for x in idxs if idxs.count(x) > 1])
if len(dupl) != 0:
    print(dupl)
    raise ValueError("duplicate model ID's found, check which version is the correct one")

paths = [path for path in paths if get_model_id(path) not in exclude]

#%%Load and create dataframe
traj_id = pd.read_csv(traj_id_path, index_col=0)

dfs = [pd.read_csv(path, index_col=0).T.assign(mod_id = idx) for idx, path in zip(idxs, paths)]
df = pd.concat(dfs).set_index("mod_id")

#filter values without a delay. Not sure if this is a very good metric
df["delay"] = df["delay"].clip(0, None)

#%%Morris
lev=4
par_morris=traj_id.columns.to_list()
seed = 230

problem = {
        "num_vars" : len(par_morris),
        "names" : par_morris,
        "groups" : None,
        "bounds" : [[0, lev-1]] * len(par_morris)
        }

X = traj_id.iloc[idxs].values.astype(float)

#%%Get morris output
delta = lev / (2 * (lev - 1))

ee={}
for var in df.columns:
    Y = df[var].values
    ee[var] = compute_elementary_effects(  #For bookkeeping purposes
        X, Y, problem["num_vars"]+1, delta)
    
    ee[var] = pd.DataFrame(data = ee[var].T, columns = problem["names"])

#%%Correct errors
#In one instance Kh_aqf was lower than Kh_aqt (trajectory 9), 
#since we prescribed Kv_aqt and due to a high Kh_Kv, this got converted to a Kh_aqt higher than Kh_aqf.
#This flipper around the elementary effects for N_aqt and f_aqt, which we correct for here.
ee["old_water"]["N_aqt"][8] = ee["old_water"]["N_aqt"][8]*-1
ee["old_water"]["f_aqt"][8] = ee["old_water"]["f_aqt"][8]*-1

#
#ee["offshore_fw"]["alpha"][6] = 0

#%%Calculate Morris metrics
output = {}
for var in df.columns:
    output[var] = pd.DataFrame()
    output[var]["inputs"] = problem["names"]
    output[var] = output[var].set_index("inputs")
    output[var]["mu"] = np.average(ee[var], axis=0)
    output[var]["mu_star"] = np.average(np.abs(ee[var]), axis=0)
    output[var]["sigma"] = np.std(ee[var], ddof=1, axis=0)
    
#%%Check monotonicity
monotone=dict([(var, (np.abs(output[var]["mu"])/output[var]["mu_star"]).values) for var in df.columns])
monotone=pd.DataFrame(monotone, index=convert_texts(output[var].index))

#%%Figsizes
agu_small = (9.5/2.54, 11.5/2.54)
agu_half  = (19/2.54, 11.5/2.54)
agu_whole = (19/2.54, 23/2.54)
agu_half_vert = (9.5/2.54, 23/2.54)

#%%Select vars to plot
output.pop("max_fw_decrease", None)

order=["onshore_sw", "offshore_fw","old_water", "s_gradient"]
vars_paper=["$S_{on}$", "$FW_{off}$", "$S_{init}$", "$S_{on}'$"]

lookup = dict(zip(order, vars_paper))

monotone = monotone[order]

#%%Groups for inputs
groups = {"geometry"     : ["l_a", "H_b", "f_H", "alpha", "beta", "phi_f"],
          "lithology"    : ["f_aqt", "l_conf", "N_aqt", "N_pal", "s_pal"],
          "hydrogeology"  : ["Kh_aqf", "Kv_aqt", "f_Kh_pal", "Kh_Kv", "R"],
          "surface water": ["N_chan", "f_chan", "l_surf_end", "l_tra", "t_tra"],
          "solute transport" : ["n", "a_l"]}

#Flip around dictionary
groups = dict([(new_key, key) for key, values in groups.items() for new_key in values])
#%%Create overall figure outline
sns.set_style("white")
ncol=3
nrow=2
width_ratios=[2,2,1]

fig = plt.figure(figsize=agu_half)
gs = GridSpec(nrow, ncol, figure=fig, width_ratios=width_ratios)

ax_locs = list(product(*([[0,1]]*2)))
axes_scatter = [fig.add_subplot(gs[i, j]) for i, j in ax_locs]
ax_mono = fig.add_subplot(gs[1, -1])

sns.despine(right=True, top=True)

#%%Plot scatter
plot_labels = list(product(["plot_y_label", "plot_x_label"], [True, False]))
plot_labels = list(product(plot_labels[:2], plot_labels[2:]))
plot_labels = plot_labels[1::2] + plot_labels[::2] #shuffle in the right order
plot_labels = [dict(d) for d in plot_labels]

exceptions = {"offshore_fw": {"$n$" : 1},
              "s_gradient" : {"$f_{aqt}$" : 0},
              "old_water" : {},
              "onshore_sw": {}}

texts_all = []
for i, var in enumerate(order):
    ax=axes_scatter[i]
    output[var]["sign"] = np.where(output[var]["mu"] < 0, "decrease", "increase")
    output[var] = output[var].assign(group = output[var].index).replace(groups)
    
    labels = convert_texts(output[var].index)
    
    kwargs = plot_labels[i]
    kwargs.update(exceptions[var])
    
    texts_all.append(labeled_scatter(fig, ax, output[var], labels, legend_out= (var == 'offshore_fw'),
                                     **kwargs))
    ax.set_title(lookup[var])

#%%Select what to plot for monotonicity
mask = pd.DataFrame(1, index=monotone.index, columns=monotone.columns)
for txts, var in zip(texts_all, order):
    for txt in txts:
        mask.loc[txt, var] = 0


#%%Select what to plot for monotonicity
texts = list(chain.from_iterable(texts_all))
texts = sorted(texts,key=texts.count,reverse=True)
texts = list(dict.fromkeys(texts)) #Get unique values, and preseverve order (requires a least Python 3.6)
monotone = monotone.loc[texts]
mask = mask.loc[texts]


#%%Plot heatmap monotonicity
monotone.columns = vars_paper
mask.columns = vars_paper
cmap = sns.light_palette("midnightblue", as_cmap=True, reverse=True)
g = sns.heatmap(monotone, cmap= cmap, linewidths=.5, ax=ax_mono, mask=mask)
g.set_yticklabels(g.get_yticklabels(), rotation=0)
g.set_xticklabels(g.get_xticklabels(), rotation=90)
g.set_title("$\epsilon$")

#%%Save
dpi_paper = 300
dpi_poster = 12.5/5 * 1.1 * dpi_paper

plt.tight_layout()
plt.savefig(outf + ".pdf")
plt.savefig(outf + ".png", dpi=dpi_paper)
plt.savefig(outf + ".svg")
plt.close()