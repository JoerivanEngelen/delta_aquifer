# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 09:43:24 2020

@author: engelen
"""

import xarray as xr
import cftime
import pandas as pd
import os
from pkg_resources import resource_filename
from glob import glob
import numpy as np
import re

import seaborn as sns
import matplotlib.pyplot as plt

#%%Functions to parse 
def get_model_id(file):     
    pattern=r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return(number)

#%%

#%%Path management
datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
res_fol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results_30deltas"

mas_f = glob(os.path.join(res_fol, "mas",  "mas_i*.nc"))
mas_pump_f = glob(os.path.join(res_fol, "mas",  "mas_pump_i*.nc"))
vol_f = glob(os.path.join(res_fol, "vols", "vol_i*.nc"))
vol_pump_f = glob(os.path.join(res_fol, "vols", "vol_pump_i*.nc"))

par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Create start and endpoints
delta_len = 9

to_finish = np.array([9, 18, 27, 63, 198, 207])

starts = np.array([0, 9, 18, 27, 36, 45, 54, 63, 81, 90, 153, 171, 198, 207, 216])
starts = starts[~np.isin(starts, to_finish)]
stops = starts+delta_len

from itertools import product

k_values = list(product(["min ", "max ", "mean "], ["max", "min", "mean"]))
k_values = pd.DataFrame([[i, j] for (i, j) in list(k_values)], columns=["Kh_aqf", "Kv_aqt"])

#%%Process to long format for seaborn
ds_deltas = {}
for start, stop in zip(starts, stops):
        delta_name = par.loc[start, "Delta"]
#        k_values = par.loc[slice(start, stop-1), ["Kh_aqf", "Kv_aqt"]].reset_index(drop=True)
        
        files_mod = [i for i in vol_f if get_model_id(i) in range(start, stop)]
        assert(len(files_mod)==delta_len)
        
        ds_deltas[delta_name] = [xr.open_dataset(f).sortby("time", ascending=True) for f in files_mod]
        ds_deltas[delta_name] = [ds.assign_coords(
                time=[t.year/1000 for t in ds.time.values[::-1]]
                ) for ds in ds_deltas[delta_name]]
        
        ds_deltas[delta_name] = [ds.to_dataframe().reset_index().assign(
                **k_values.loc[i].to_dict()
                ) for i, ds in enumerate(ds_deltas[delta_name])]
        
        ds_deltas[delta_name] = pd.concat(ds_deltas[delta_name])

df_deltas = pd.concat([ds.assign(delta=delta) for delta, ds in ds_deltas.items()])

#%%
sns.set(style="whitegrid")
paper_rc = {'lines.linewidth': 3}                  
sns.set_context(rc = paper_rc, font_scale=4)                                    

g = sns.relplot(x="time", y="fw_onshore", hue="Kh_aqf", style = "Kv_aqt", 
             col="delta", col_wrap = 4, data=df_deltas, kind="line", 
             palette = "Blues", hue_order = ["min ", "mean ", "max "],
             size=1, aspect = 1.05, 
             facet_kws=dict(sharey=False))

g.savefig(os.path.join(res_fol, "fw_onshore.png"), dpi=300)