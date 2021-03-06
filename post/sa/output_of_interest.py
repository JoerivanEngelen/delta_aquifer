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
import re

def spat_sum(da, name):
    return(da.sum(dim=["x", "y", "z"]).rename(name))

def cftime_to_ka(ds):
    years = np.array([t.year for t in ds.time.values])
    years = (years - years[-1])/1000 * -1
    
    return(ds.assign_coords(time=years))

def coord_of_max(da):
    return(da.where(da==da.max(), drop=True)[-1].time)

#%%Ideas
#Upward salt flux? -> Check if can even occur in this model
#How much recharge goes in, how much discharged?
#   -> qz*A/rch

#Average slope of the interface?

#%%Path management
if len(sys.argv) > 1:
    fol  = sys.argv[1]
else:
    #Local testing on my own windows laptop
    fol = r"g:\synthdelta\results\test_output\synth_SD_i023_m24_7071784"

nc_paths = glob(os.path.join(fol, "results_*.nc"))
nc_paths.sort()

bc_path = os.path.join(fol, "input", "data", "bcs.nc")

match = re.search(r"(i[0-9]{3}).+([0-9]{7})", os.path.basename(fol))
oi_name = match.group(0)

fracs_path = os.path.join(fol, "fracs_%s.nc" % oi_name)
oi_path = os.path.join(fol, "oi_%s.csv" % oi_name)

#%%pre-process data for analysis
sea_level = xr.open_dataset(bc_path)["sea_level"]
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
if np.all(onshore==0.):
    onshore = active.isel(time=-1).sel(z=0, method="backfill")
x_loc=xr.where(onshore, onshore.x, np.nan).max(dim="x")
x_loc=x_loc.fillna(x_loc.min()) #fillna

#%%Calculate vols and masses
fw_vols  = xr.where((ds["conc1"] > -1.) & (ds["conc1"] < 10.),  ds["vol"], 0)
fw_vols_offshore = xr.where((fw_vols!=0) & (ds.x > x_loc), fw_vols, 0)
bw_vols = xr.where((ds["conc1"] > 10.)&(ds["conc1"] < 30.), 
                               ds["vol"], 0)
bw_vols_onshore = xr.where((bw_vols!=0) & (ds.x < x_loc), bw_vols, 0)

tot_vols = xr.where(active, ds["vol"], 0)
tot_offshore_vols = (ds.x > x_loc) * tot_vols

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
vol["tot_offshore"]= spat_sum(tot_offshore_vols, "tot_offshore")
vol["tot_onshore"] = vol["tot"] - vol["tot_offshore"]
vol["fw"]          = spat_sum(fw_vols, "fw")
vol["fw_offshore"] = spat_sum(fw_vols_offshore, "fw_offshore")
vol["fw_onshore"]  = vol["fw"] - vol["fw_offshore"]
vol["bw"]          = spat_sum(bw_vols, "bw")
vol["bw_onshore"]  = spat_sum(bw_vols_onshore, "bw_onshore") 
vol["sw"]          = spat_sum(xr.where(ds["conc1"] > 30., ds["vol"], 0), "sw")
vol["ow"]          = spat_sum(xr.where(ds["conc2"] > 0.5, ds["vol"], 0), "ow")
vol=xr.Dataset(data_vars=vol)

frac_vols = vol[["fw", "bw", "sw", "ow"]]/vol["tot"]
frac_vols["fw_offshore"] = vol["fw_offshore"]/vol["tot_offshore"]
frac_vols["fw_onshore"] = vol["fw_onshore"]/vol["tot_onshore"]
frac_vols["bw_onshore"] = vol["bw_onshore"]/vol["tot_onshore"]

frac_mas = {}
frac_mas["m_sal"] = mas["sal"]/mas["tot"]
frac_mas["m_ol_sal"] = mas["ol_sal"]/mas["sal"]
frac_mas["m_sal_onshore"] = mas["sal_onshore"]/mas["tot_onshore"]
frac_mas=xr.Dataset(data_vars=frac_mas)

#%%
fracs = xr.merge([frac_mas, frac_vols])
fracs.to_netcdf(fracs_path)