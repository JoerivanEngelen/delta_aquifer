# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:13:46 2020

@author: engelen
"""

import pandas as pd
import numpy as np

import os
from pkg_resources import resource_filename

#%%Path management

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)
outfol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\paper_40_deltas\inputs"
outfile = os.path.join(outfol, "table_inputs.csv")
#%%Process

par = par.drop(columns=["index", "dx", "dy", "riv_res", "Kh_Kv", "gamma"])

gb = par.groupby("Delta")

table = gb.mean()
table["Kh_aqf"] = pd.IntervalIndex([pd.Interval(i, j) for i, j in zip(gb.min()["Kh_aqf"], gb.max()["Kh_aqf"])])
table["Kv_aqt"] = pd.IntervalIndex([pd.Interval(i, j) for i, j in zip(gb.min()["Kv_aqt"], gb.max()["Kv_aqt"])])

table["La"] = table["L"] * table["l_a"]
table["Lb"] = table["L"] - table["La"]
table["H_a"] = table["H_b"] * table["f_H"]

table = table.drop(columns=["L", "l_a", "f_H"])

#%%Prepare in format useful for paper
table["La"] = table["La"]/1000.
table["Lb"] = table["Lb"]/1000.
table["t_tra"] = table["t_tra"]/1000.

table["phi_f"] = table["phi_f"]/np.pi

table.loc[:, ["La", "Lb", "H_b", "H_a"]] = table.loc[:, ["La", "Lb", "H_b", "H_a"]].round(decimals = 0)
table.loc[:, ["phi_f","l_tra", "l_surf_end"]] = table.loc[:, ["phi_f", "l_tra", "l_surf_end"]].round(decimals = 2)

#%%Reorder columns
col = ["La", "Lb", "alpha", "beta", "H_a", "H_b", "phi_f",
       "N_aqt", "f_aqt", "l_conf", "N_pal", "s_pal", "Kh_aqf", "Kv_aqt", "n", 
       "l_tra", "t_tra", "N_chan", "l_surf_end", "R"]

table = table[col]

#%%Write
table.to_csv(outfile)