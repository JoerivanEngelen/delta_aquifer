# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:12:46 2019

@author: engelen
"""

import pandas as pd
import seaborn as sns
from glob import glob
import os
import re
from SALib.analyze.morris import analyze
from pkg_resources import resource_filename

import matplotlib.pyplot as plt
import numpy as np

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
           "Kh_Kv" : r"Kh/Kv"}
    
    texts = [tex[t] if t in tex.keys() else t for t in texts]
    texts = [str_to_mathtxt(s) for s in texts]
    return(texts)

def labeled_scatter(ax, xs, ys, labels):
    ax.scatter(xs, ys)
    
    dx=np.diff(np.array(ax.get_xlim()))
    dy=np.diff(np.array(ax.get_ylim()))
    
    for x, y, label in zip(xs, ys, labels):
        if x > (0.1 * dx):
            ax.text(x+0.01*dx, y+0.01*dy, label)
    return(ax)

#%%Path management
path = r"g:\synthdelta\results\outputs_of_interest"
traj_id_path = os.path.abspath(resource_filename("delta_aquifer", "../data/traj_id.csv"))
outf = os.path.join(path, "..", "morris_covariance.png")

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
output={}
for var in df.columns:
    Y = df[var].values
    output[var]=analyze(problem, X, Y, num_resamples=1000, conf_level=0.95, 
                           print_to_console=False, num_levels=lev, seed=seed)

#%%Figsizes
agu_small = (9.5/2.54, 11.5/2.54)
agu_half  = (19/2.54, 11.5/2.54)
agu_whole = (19/2.54, 23/2.54)
agu_half_vert = (9.5/2.54, 23/2.54)

#%%Plot
sns.set_style("whitegrid")
ncol=3
nrow=3

fig, axes = plt.subplots(ncol, nrow, figsize=agu_whole)

sns.despine(right=True, top=True)

#Onshore_sw
#df.loc[df["onshore_sw"].isnull()]

for i, var in enumerate(output.keys()):
    ax=axes.flatten()[i]
    xs = output[var]["mu_star"]
    ys = output[var]["sigma"]
    labels = convert_texts(output[var]["names"])
    
    labeled_scatter(ax, xs, ys, labels)
    ax.set_title(var)

for j in range(i+1, (ncol*nrow)):
    axes.flatten()[j].axis("off")

plt.tight_layout()
#%%Save

plt.savefig(outf)