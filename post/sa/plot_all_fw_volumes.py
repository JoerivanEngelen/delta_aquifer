# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 14:26:12 2019

@author: engelen
"""

import re
import sys, os
from glob import glob
import xarray as xr
import matplotlib.pyplot as plt

def get_model_id(file):     
    pattern=r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return(number)

#%%Ideas
#Find radius of curvature of line, this says something of the smootheness of the line, a small radius means sharp angles.
#See response of Jugad here: https://www.reddit.com/r/Python/comments/3lmc30/is_there_a_way_to_calculate_the_smoothness_of_a/
#http://mathworld.wolfram.com/RadiusofCurvature.html

#%%Set paths
path = r"g:\synthdelta\results\trends"

#%%Process path
#Unfinished trajectories
exclude=list(range(96, 120)) + list(range(216, 240))

paths = glob(os.path.join(path, "*.nc"))
paths.sort()
idxs = [get_model_id(path) for path in paths]
idxs = [idx for idx in idxs if idx not in exclude]
paths = [path for path in paths if get_model_id(path) not in exclude]

#%%
ds_all = [xr.open_dataset(path) for path in paths]
fws = [ds["fw"].rename(i) for ds, i in zip(ds_all, idxs)]
fws = xr.merge(fws).drop(["dx", "dy"]).to_dataframe()

#%%

agu_small = (9.5/2.54, 11.5/2.54)
agu_half  = (19/2.54, 11.5/2.54)

fig, ax = plt.subplots(figsize=agu_half)

for i in fws.columns:
    fws[i].plot(color="black", alpha=0.1)

ax.invert_xaxis()
ax.set_ylabel("FW vol [-]")
ax.set_xlabel("time [ka]")

plt.savefig(os.path.join(path, "..", "fw_volume_trends.png"), dpi=300)
