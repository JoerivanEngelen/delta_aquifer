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


from itertools import product

import seaborn as sns
import matplotlib.pyplot as plt

#%%Functions to parse 
def get_model_id(file):     
    pattern=r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return(number)

#%%Function to process data
def t_to_ka(ds):
    return([t.year/1000 for t in ds.time.values[::-1]])

#%%Function to adapt plot
def label_offset(ax, axis="y"):
    if axis == "y":
        fmt = ax.yaxis.get_major_formatter()
        ax.yaxis.offsetText.set_visible(False)
        set_label = ax.set_ylabel
#        label = ax.get_ylabel()
        label = "$V_{fw}$"
    
    elif axis == "x":
        fmt = ax.xaxis.get_major_formatter()
        ax.xaxis.offsetText.set_visible(False)
        set_label = ax.set_xlabel
        label = ax.get_xlabel()

    def update_label(event_axes):
        offset = fmt.get_offset()
        if offset == '':
            set_label("{} (m3)".format(label))
        else:
            minus_sign = u'\u2212'
            expo = np.float(offset.replace(minus_sign, '-').split('e')[-1])
            set_label("%s ($\mathregular{10^{%d}}$ m3)" % (label, expo))    

    ax.callbacks.connect("ylim_changed", update_label)
    ax.callbacks.connect("xlim_changed", update_label)
    ax.figure.canvas.draw()
    update_label(None)

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

to_finish = np.array([9, 18, 63, 198, 207])

starts = np.array([0, 9, 18, 27, 36, 45, 54, 63, 81, 90, 153, 171, 198, 207, 216])
starts = starts[~np.isin(starts, to_finish)]
stops = starts+delta_len

k_values = list(product(["min ", "max ", "mean "], ["max", "min", "mean"]))
k_values = pd.DataFrame([[i, j] for (i, j) in list(k_values)], columns=["Kh_aqf", "Kv_aqt"])

#%%Process to long format for seaborn
ds_deltas = {}
for start, stop in zip(starts, stops):
        delta = par.loc[start, "Delta"]
        
        files_mod = [i for i in vol_f if get_model_id(i) in range(start, stop)]
        assert(len(files_mod)==delta_len)
        
        ds_deltas[delta] = [xr.open_dataset(f).sortby("time", ascending=True) for f in files_mod]
        ds_deltas[delta] = [ds.assign_coords(time=t_to_ka(ds)) for ds in ds_deltas[delta]]
        
        ds_deltas[delta] = [ds.to_dataframe().reset_index().assign(
                **k_values.loc[i].to_dict()
                ) for i, ds in enumerate(ds_deltas[delta])]
        
        ds_deltas[delta] = pd.concat(ds_deltas[delta])

df_deltas = pd.concat([ds.assign(delta=delta) for delta, ds in ds_deltas.items()])
df_deltas = df_deltas.rename(columns={"time" : "time (ka)"})

#%%
agu_whole = (19/2.54, 23/2.54)
col_wrap = 4

height = agu_whole[1]/col_wrap
aspect = agu_whole[0]/agu_whole[1]

sns.set(style="whitegrid")                               

g = sns.relplot(x="time (ka)", y="fw_onshore", hue="Kh_aqf", style = "Kv_aqt", 
             col="delta", col_wrap = col_wrap, data=df_deltas, kind="line", 
             palette = "Blues", hue_order = ["min ", "mean ", "max "],
             style_order = ["max", "mean", "min"], 
             height=height, aspect = aspect, 
             facet_kws=dict(sharey=False, legend_out=False, 
                            margin_titles=False, xlim=(125, 0)))

for ax in g.axes:
    ax.set_ylim(ymin=0)
    label_offset(ax, axis="y")
    ax.invert_xaxis()

g.axes[0].legend(loc='lower right', bbox_to_anchor=(0.975, 0.025), 
      bbox_transform=g.fig.transFigure)

plt.subplots_adjust(hspace=0.25, wspace=0.8)
g.set_titles(template="{col_name}")
g.savefig(os.path.join(res_fol, "fw_onshore.png"), dpi=300)