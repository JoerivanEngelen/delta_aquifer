# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 13:52:56 2019

@author: engelen
"""

import pandas as pd
from pkg_resources import resource_filename
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#%%Path management
datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))

class_path = os.path.join(datafol, "model_classification.csv")
traj_id_path  = os.path.join(datafol, "traj_id.csv")
traj_real_path = os.path.join(datafol, "traj_real.csv")
constants_path = os.path.join(datafol, "fixed_pars.csv")

outf = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results\class_swarm"

#%%Read and select
clas = pd.read_csv(class_path, skiprows=10, index_col=0)
traj_id = pd.read_csv(traj_id_path, index_col=0)
traj_real = pd.read_csv(traj_real_path, index_col=0)
constants = pd.read_csv(constants_path, index_col=0)

filter_nan = ~np.isnan(clas["type"])
traj_real = traj_real[filter_nan]
traj_id = traj_id[filter_nan]
clas = clas[filter_nan]
clas = clas.astype(int)

df = pd.DataFrame()
df["type"] = clas["type"]

#%%Reclassify to change order

#Swap numbers, so that:
#1: Dispersion dominant
#2: Fingers
#3: Mix of 1 and 2
#4: Sharp interface
#5: Overtop
df["type"] = df["type"].replace({1:1, 4:2, 5:3, 3:4, 2:5})


#%%Calculate
df["rs"] = traj_real["H_b"] * traj_real["f_aqt"] / traj_real["N_aqt"] / traj_real["Kv_aqt"]
#df["rs"] = traj_real["H_b"] * traj_real["f_aqt"] / traj_real["Kv_aqt"]
df["rs"].loc[np.isinf(df["rs"])] = 0
df["rs"] += traj_real["H_b"] * (1-traj_real["f_aqt"]) / (traj_real["N_aqt"]+1) / (traj_real["Kh_aqf"] / traj_real["Kh_Kv"])
#df["rs"] += traj_real["H_b"] * (1-traj_real["f_aqt"]) / (traj_real["Kh_aqf"] / traj_real["Kh_Kv"])
df["rs"] = np.log10(df["rs"])

q_f = traj_real["alpha"] * traj_real["Kh_aqf"]
Q_f = traj_real["H_b"] * q_f
df["Q"] = np.log10(Q_f)

peclet = (constants["D_m"].values * traj_real["n"] + traj_real["a_l"]*q_f) / Q_f

df["peclet"] = np.log10(peclet)
df["H_b"] = np.log10(traj_real["H_b"])

#df["Kh_aqf"] = np.round(np.log10(traj_real["Kh_aqf"]), decimals=2)
df["Kh_aqf"] = np.round(traj_real["Kh_aqf"], decimals=1)

#%%Figsizes
agu_small = (9.5/2.54, 11.5/2.54)
agu_half  = (19/2.54, 11.5/2.54)
agu_whole = (19/2.54, 23/2.54)
agu_half_vert = (9.5/2.54, 23/2.54)

#%%Plot

dark = ".1"
light = ".95"

style_dict = {"axes.facecolor": dark, 
               "axes.labelcolor" : dark,
               'text.color': light, 
               'figure.facecolor': dark, 
               'xtick.color': light,
               'ytick.color': light}

with sns.axes_style("dark", style_dict):
    cmap = sns.color_palette("Paired", 5)     
    
    fig, ax1 = plt.subplots(figsize=agu_half)
    
    sns.swarmplot(x="Kh_aqf", y="rs", hue="type", data=df, 
                  dodge=False, palette = cmap,
                  linewidth=0.1, edgecolor=dark)
    
    ax1.set_xlabel("$K_{h, aqf}$ [m/d]")
    ax1.set_ylabel("$r_s$ [d]")

    ax1.get_yaxis().set_major_formatter(ticker.FormatStrFormatter("$10^{%d}$"))

    ax1.xaxis.label._color = light
    ax1.yaxis.label._color = light
    
    K_values = np.unique(df["Kh_aqf"])
    
    for x in range(len(K_values)-1):
        plt.axvline(x+0.5, color=light, linestyle=":") #Datacoordinates catagorical = 0, 1, 2, etc.
    
#    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1.02))
    ax1.legend(loc='upper left', bbox_to_anchor=(-0.5, 1.02))
    plt.tight_layout()
    
    plt.savefig(outf + ".png", facecolor=fig.get_facecolor(), edgecolor='none', dpi=300)
    plt.savefig(outf + ".svg", facecolor=fig.get_facecolor(), edgecolor='none')
