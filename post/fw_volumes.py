# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:27:13 2019

@author: engelen
"""

import os, sys
from glob import glob
import xarray as xr
import numpy as np
from delta_aquifer import initial_conditions as ic

if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt

import seaborn as sns

#%%Path management
if len(sys.argv) > 1:
    fol  = sys.argv[1]
else:
    #Local testing on my own windows laptop
    fol = r"g:\synthdelta\results\test_output\synth_SD_i029_m24_7065039"

nc_paths = glob(os.path.join(fol, "results_*.nc"))
nc_paths.sort()

#%%pre-process data for analysis
ds = xr.open_mfdataset(nc_paths, 
                       chunks={"time" : 2, "z": -1, "y": -1, "x": -1},
                       combine="nested",
                       concat_dim="time")
ds["vol"] = np.abs(ds.dx * ds.dy * ds.dz)

years = np.array([t.year for t in ds.time.values])
years = (years - years[-1])/1000 * -1

ds = ds.assign_coords(time=years)

#%%Get coastline
active = (ds["conc1"] > -1.)
onshore = active.sel(z=0, method="nearest")
x_loc=xr.where(onshore, onshore.x, np.nan).max(dim="x")
x_loc=x_loc.fillna(x_loc.min()) #fillna

#%%Calculate fw_vols
fw_vols  = xr.where((ds["conc1"] > -1.) & (ds["conc1"] < 1.),  ds["vol"], 0)
tot_vols = xr.where(active, ds["vol"], 0)
fw_vols_offshore = xr.where((fw_vols!=0) & (ds.x > x_loc), fw_vols, 0)

brsal_onshore=xr.where(
        (ds["conc1"] > 0.) & (ds.x < x_loc), 
        ds["conc1"]*ds["vol"], 0
        ).sum(dim=["x", "y", "z"]
        ).rename("sbw_onshore")

fw_vol  = fw_vols.sum(dim=["x", "y", "z"]).rename("fw")
fw_vol_offshore = fw_vols_offshore.sum(dim=["x", "y", "z"]).rename("fw_offshore")
br_vol  = xr.where((ds["conc1"] > 1.)  & (ds["conc1"] < 30.), ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("bw")
sal_vol = xr.where(ds["conc1"] > 30., ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("sw")
old_vol = xr.where(ds["conc2"] > 0.5, ds["vol"], 0).sum(dim=["x", "y", "z"]).rename("ow")

tot_vol = tot_vols.sum(dim=["x", "y", "z"])
tot_masses = (ic.c2dens(ds["conc1"]*active)*tot_vols)
tot_mas = tot_masses.sum(dim=["x", "y", "z"])
tot_mass_onshore = ((ds.x > x_loc)*tot_masses).sum(dim=["x", "y", "z"])

frac_vols = xr.merge([fw_vol, br_vol, sal_vol, old_vol, fw_vol_offshore])/tot_vol
frac_mas = xr.merge([brsal_onshore])/tot_mass_onshore

frac_vols.to_netcdf(os.path.join(fol, "frac_vols.nc"))

#%%Outputs of interest
end_fw      = frac_vols["fw"].isel(time=-1).values
max_change  = frac_vols["fw"].max().values - end_fw
#old water should be old water mass normalized over total salt mass instead of volume ratio
old_water   = (frac_vols["ow"]/frac_vols["sw"]).isel(time=-1).values
onshore_sw  = frac_mas["sbw_onshore"].isel(time=-1).values
offshore_fw = frac_vols["fw_offshore"].isel(time=-1).values

#%%Plot
sns.set_style("darkgrid", {"axes.facecolor": ".4", "xtick.bottom": True})

lines = []
labels = ["fw", "bw", "sw", "ow"]
colors = [[0.5703125, 0.7734375, 0.87109375],
          [0.99609375, 0.87109375, 0.6015625],
          [0.9609375, 0.5625, 0.32421875],
          [0.1, 0.1, 0.1]]

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

plt.savefig(os.path.join(fol, "fw_volume.png"), dpi=300)