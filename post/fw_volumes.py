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
import pandas as pd

if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt

import seaborn as sns

def spat_sum(da, name):
    return(da.sum(dim=["x", "y", "z"]).rename(name))

def cftime_to_ka(ds):
    years = np.array([t.year for t in ds.time.values])
    years = (years - years[-1])/1000 * -1
    
    return(ds.assign_coords(time=years))

def coord_of_max(da):
    return(da.where(da==da.max(), drop=True).squeeze())

#%%Ideas
#Upward salt flux? -> Check if can even occur in this model
#How much recharge goes in, how much discharged?
#   -> qz*A/rch

#%%Path management
if len(sys.argv) > 1:
    fol  = sys.argv[1]
else:
    #Local testing on my own windows laptop
    fol = r"g:\synthdelta\results\test_output\synth_SD_i029_m24_7065039"

nc_paths = glob(os.path.join(fol, "results_*.nc"))
nc_paths.sort()

bc = os.path.join(fol, "input", "data", "bcs.nc")

#%%pre-process data for analysis
sea_level = xr.open_dataset(bc)["sea_level"]
sea_level = cftime_to_ka(sea_level)

ds = xr.open_mfdataset(nc_paths, 
                       chunks={"time" : 2, "z": -1, "y": -1, "x": -1},
                       combine="nested",
                       concat_dim="time")
ds["vol"] = np.abs(ds.dx * ds.dy * ds.dz)
ds = cftime_to_ka(ds)

#%%Get coastline
active = (ds["conc1"] > -1.)
#To have a moving coastline through time, I have not enough info
#For this I probably have to store coastline_x in bcs.nc and read this
#Right now not a problem because I'm looking at the endstate offshore/onshore stuff
onshore = active.isel(time=-1).sel(z=0, method="nearest")
x_loc=xr.where(onshore, onshore.x, np.nan).max(dim="x")
x_loc=x_loc.fillna(x_loc.min()) #fillna

#%%Calculate vols and masses
fw_vols  = xr.where((ds["conc1"] > -1.) & (ds["conc1"] < 1.),  ds["vol"], 0)
fw_vols_offshore = xr.where((fw_vols!=0) & (ds.x > x_loc), fw_vols, 0)

tot_vols = xr.where(active, ds["vol"], 0)
tot_masses = ic.c2dens(ds["conc1"]*active)*tot_vols
sal_masses=xr.where(
        (ds["conc1"] > 0.), 
        ds["conc1"]*ds["vol"], 0
        )

mas={}
mas["sal"]         = spat_sum(sal_masses, "sal")
mas["ol_sal"]      = spat_sum(sal_masses*ds["conc2"], "ol_sal")
mas["sal_onshore"] = spat_sum(xr.where((ds.x < x_loc), sal_masses, 0), "sal_onshore")
mas["tot"]         = spat_sum(tot_masses, "tot")
mas["tot_onshore"] = spat_sum((ds.x < x_loc)*tot_masses, "tot_onshore")
mas=xr.Dataset(data_vars=mas)

vol={}
vol["tot"]         = spat_sum(tot_vols,"tot")
vol["fw"]          = spat_sum(fw_vols, "fw")
vol["fw_offshore"] = spat_sum(fw_vols_offshore, "fw_offshore")
vol["bw"]          = spat_sum(xr.where((ds["conc1"] > 1.)&(ds["conc1"] < 30.), 
                               ds["vol"], 0), "bw")
vol["sw"]          = spat_sum(xr.where(ds["conc1"] > 30., ds["vol"], 0), "sw")
vol["ow"]          = spat_sum(xr.where(ds["conc2"] > 0.5, ds["vol"], 0), "ow")
vol=xr.Dataset(data_vars=vol)

frac_vols = vol/vol["tot"]
frac_mas = mas[["sal", "ol_sal", "sal_onshore"]]
frac_mas["sal"] = mas["sal"]/mas["tot"]
frac_mas["ol_sal"] = mas["ol_sal"]/mas["sal"]
frac_mas["sal_onshore"] = mas["sal_onshore"]/mas["tot_onshore"]

frac_vols.to_netcdf(os.path.join(fol, "frac_vols.nc"))
frac_mas.to_netcdf( os.path.join(fol, "frac_mas.nc"))

#%%Read again to remove chunking
frac_vols = xr.open_dataset(os.path.join(fol, "frac_vols.nc"))
frac_mas  = xr.open_dataset(os.path.join(fol, "frac_mas.nc"))

#%%Differentiate
grad_fw = frac_vols["fw"].differentiate("time") * -1 #Multiply with -1 because the time axis is decreasing
grad_sl = sea_level.differentiate("time") * -1 

#%%Outputs of interest
oi = {}
oi["end_fw"]          = frac_vols["fw"].isel(time=-1).values
oi["offshore_fw"]     = frac_vols["fw_offshore"].isel(time=-1).values
oi["max_fw_decrease"] = frac_vols["fw"].max().values - oi["end_fw"] 
oi["old_water"]       = frac_mas["ol_sal"].isel(time=-1).values
oi["onshore_sw"]      = frac_mas["sal_onshore"].isel(time=-1).values/(35./1025.)
oi["fw_gradient"]     = (grad_fw).isel( 
                            time=slice(-3, None)
                            ).mean()
oi["delay"] = coord_of_max(grad_sl).time - coord_of_max(grad_fw*-1).time 

keys, values = list(zip(*oi.items()))
oi = pd.DataFrame(data={"var" : keys, "value" : values}).set_index("var")

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