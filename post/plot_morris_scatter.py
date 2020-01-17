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
import itertools

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

from adjustText import adjust_text

#%%Functions to parse
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

def labeled_scatter(ax, xs, ys, labels):
    ax.scatter(xs, ys, alpha=0.7)
    
    x_max=np.nanmax(xs)
    y_max=np.nanmax(ys)
    x_max += 0.1 * x_max
    y_max += 0.1 * y_max
    ax.set_xlim(left=0.0, right=x_max)
    ax.set_ylim(bottom=0.0, top=y_max)
    
    x_line = np.linspace(0, x_max)
    ax.plot(x_line, x_line, ls=":", alpha=0.5)
    
    ax.set_xlabel("$\mu *$")
    ax.set_ylabel("$\sigma$")
    
    dx=np.diff(np.array(ax.get_xlim()))
    dy=np.diff(np.array(ax.get_ylim()))
    
    texts=[]
    
    for x, y, label in zip(xs, ys, labels):
        if x > (0.3 * dx):
            texts.append(ax.text(x, y, label, va="center", ha="center"))
        elif y > (0.44 * dy):
            texts.append(ax.text(x, y, label, va="center", ha="center"))
    
    adjust_text(texts, ax=ax, expand_points=(1.7, 2), force_text=(1.0, 1.0),  
                arrowprops=dict(arrowstyle="wedge", color="k", alpha=0.2, lw=0.1))
#                arrowprops=dict(arrowstyle="wedge", color="lightgray", lw=0.1))
    
    texts = [t._text for t in texts]
    
    return(texts)

#%%Path management
path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results\outputs_of_interest"
traj_id_path = os.path.abspath(resource_filename("delta_aquifer", "../data/traj_id.csv"))
outf = os.path.join(path, "..", "morris")

#%%Process paths
#Unfinished trajectories
exclude=list(range(96, 120)) + list(range(216, 240))

paths = glob(os.path.join(path, "*.csv"))
paths.sort()
idxs = [get_model_id(path) for path in paths]
idxs = [idx for idx in idxs if idx not in exclude]
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

output={}
ee={}
for var in df.columns:
    Y = df[var].values
    output[var]=analyze(problem, X, Y, num_resamples=1000, conf_level=0.95, 
                           print_to_console=False, num_levels=lev, seed=seed)

    ee[var] = compute_elementary_effects(  #For bookkeeping purposes
        X, Y, problem["num_vars"]+1, delta)

#%%Check monotonicity
monotone=dict([(var, np.abs(output[var]["mu"])/output[var]["mu_star"]) for var in df.columns])
monotone=pd.DataFrame(monotone, index=convert_texts(output[var]["names"]))

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

#%%Plot scatter
sns.set_style("white")
ncol=3
nrow=2

fig = plt.figure(figsize=agu_half)
gs = GridSpec(nrow, ncol, figure=fig)

ax_locs = list(itertools.product(*([[0,1]]*2)))
axes_scatter = [fig.add_subplot(gs[i, j]) for i, j in ax_locs]
ax_mono = fig.add_subplot(gs[:, -1])

sns.despine(right=True, top=True)

texts_all = []

for i, var in enumerate(order):
    ax=axes_scatter[i]
    xs = output[var]["mu_star"]
    ys = output[var]["sigma"]
    negativ = output[var]["mu"] < 0
    
    labels = convert_texts(output[var]["names"])
    
    texts_all.append(labeled_scatter(ax, xs, ys, labels))
    ax.set_title(lookup[var])

#%%Select what to plot for monotonicity

mask = pd.DataFrame(1, index=monotone.index, columns=monotone.columns)
for txts, var in zip(texts_all, order):
    for txt in txts:
        mask.loc[txt, var] = 0


#%%Select what to plot for monotonicity
texts = list(itertools.chain.from_iterable(texts_all))
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
g.set_xticklabels(g.get_xticklabels(), rotation=45)
g.set_title("$\epsilon$")

#%%Save
dpi_paper = 300
dpi_poster = 12.5/5 * 1.1 * dpi_paper

plt.tight_layout()
plt.savefig(outf + ".pdf")
plt.savefig(outf + ".png", dpi=dpi_paper)
plt.savefig(outf + ".svg")
plt.close()