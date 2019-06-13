# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:27:13 2019

@author: engelen
"""

import os, sys
from glob import glob
import xarray as xr
import numpy as np

if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt

import seaborn as sns

#%%Path management
fol = sys.argv[1]
#fol = r"g:\synthdelta\test_output\SD_i202\synth_SD_i202_m24_6378797"

nc_paths = glob(os.path.join(fol, "results_*.nc"))
nc_paths.sort()

#%%pre-process data for analysis
ds = xr.open_mfdataset(nc_paths, chunks={"time" : 2, "z": -1, "y": -1, "x": -1})
ds["vol"] = np.abs(ds.dx * ds.dy * ds.dz)

years = np.array([t.year for t in ds.time.values])
years = (years - years[-1])/1000 * -1

ds = ds.assign_coords(time=years)

#%%Calculate fw_vols
fw_vol  = xr.where((ds["conc"] > -1.) & (ds["conc"] < 1.),  ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("fw")
br_vol  = xr.where((ds["conc"] > 1.)  & (ds["conc"] < 30.), ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("bw")
sal_vol = xr.where(ds["conc"] > 30., ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("sw")
tot_vol = xr.where(ds["conc"] > -1., ds["vol"], 0).sum(dim=["x", "y", "z"])

frac_vols = xr.merge([fw_vol, br_vol, sal_vol])/tot_vol

frac_vols.to_netcdf(os.path.join(fol, "frac_vols.nc"))

#%%Plot
sns.set_style("darkgrid", {"axes.facecolor": ".4", "xtick.bottom": True})

lines = []
labels = ["fw", "bw", "sw"]
colors = [[0.5703125, 0.7734375, 0.87109375],
          [0.99609375, 0.87109375, 0.6015625],
          [0.9609375, 0.5625, 0.32421875]]

for i, lab in enumerate(labels):
    lines.extend(frac_vols[lab].plot())
    lines[i].set_label(lab)
    lines[i].set_color(colors[i])

ax = plt.gca()

ax.invert_xaxis()
ax.xaxis.grid(False)

plt.title(None)
plt.ylabel("frac. volume (-)")
plt.ylim((0, 1))
plt.xlabel("time (ka)")
plt.legend()

plt.savefig(os.path.join(fol, "fw_volume.png"))