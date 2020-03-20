# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 10:20:32 2020

@author: engelen
"""

import xarray as xr
import os, sys
from glob import glob
import cftime
import pandas as pd
import numpy as np

#%%Functions
def check_sub_fol(results_fol, sub_fol):
    """Check if only 1 instance of sub_fol exists, 
    and if so return it.
    """
    
    sub_fols = glob(os.path.join(results_fol, sub_fol))
    n_results = len(sub_fols)
    
    if n_results != 1:
        raise ValueError("Should be 1 folder, instead got {} folders. " +
                         "Namely: {}".format((n_results, sub_fols)))
    
    return(sub_fols[0])

def get_coastline(active):
    #To have a moving coastline through time, I have not enough info
    #For this I probably have to store coastline_x in bcs.nc and read this
    #Right now not a problem because I'm looking at the endstate offshore/onshore stuff
    onshore = active.sel(z=0, method="nearest")
    if np.all(onshore==0.):
        onshore = active.sel(z=0, method="backfill")
    x_loc=xr.where(onshore, onshore.x, np.nan).max(dim="x")
    x_loc=x_loc.fillna(x_loc.min()) #fillna
    return(x_loc)

def spat_sum(da, name):
    return(da.sum(dim=["x", "y", "z"]).rename(name))

def get_mas_and_vol(ds, onshore, active):
    fw_vols  = xr.where(active & (ds["conc1"] < 5.),
                        ds["vol"], 0)
    
    fw_vols_offshore = fw_vols.where(~onshore, 0)
    
    tot_vols = xr.where(active, ds["vol"], 0)
    tot_offshore_vols = (~onshore) * tot_vols
    
    sal_masses = xr.where(
            (ds["conc1"] > 0.), 
            ds["conc1"]*ds["vol"], 0
            )
    
    ol_masses = sal_masses*ds["conc2"]
    
    mas={}
    mas["sal"]         = spat_sum(sal_masses, "sal")
    mas["ol_sal"]      = spat_sum(ol_masses, "ol_sal")
    mas["sal_onshore"] = spat_sum(xr.where(onshore, 
                           sal_masses, 0), "sal_onshore")
    mas["ol_sal_onshore"] = spat_sum(xr.where(onshore, 
                               ol_masses, 0), "ol_sal_onshore")
    mas["sal_offshore"] = mas["sal"] - mas["sal_onshore"]
    mas["ol_sal_offshore"] = mas["ol_sal"] - mas["ol_sal_onshore"]
    mas = xr.Dataset(data_vars=mas)
    
    vol={}
    vol["tot"]         = spat_sum(tot_vols,"tot")
    vol["tot_offshore"]= spat_sum(tot_offshore_vols, "tot_offshore")
    vol["tot_onshore"] = vol["tot"] - vol["tot_offshore"]
    vol["fw"]          = spat_sum(fw_vols, "fw")
    vol["fw_offshore"] = spat_sum(fw_vols_offshore, "fw_offshore")
    vol["fw_onshore"]  = vol["fw"] - vol["fw_offshore"]
    vol = xr.Dataset(data_vars=vol)
    return(mas, vol)

#%%Input control
if len(sys.argv)>1:
    res_fol = sys.argv[1]
    sim_nr = int(sys.argv[2])
else:
    res_fol = r"c:\Users\engelen\test_imodpython\test_results_30_deltas"
    sim_nr = 157

#%%Path management
nat_fol=check_sub_fol(res_fol, "synth_RD_i{:03d}_m24_*".format(sim_nr))
ant_fol=check_sub_fol(res_fol, "synth_AD_i{:03d}_m24_*".format(sim_nr))

nat_results=glob(os.path.join(nat_fol, "results_[0-9][0-9][0-9].nc"))
ant_results=glob(os.path.join(ant_fol, "results_[0-9][0-9][0-9].nc")) 

#%%Read
ds = xr.open_mfdataset(nat_results, 
                       chunks={"time" : 2, "z": -1, "y": -1, "x": -1},
                       combine="nested",
                       concat_dim="time")

ds_ant = xr.open_dataset(ant_results[0])

#%%Fix times anthropocene
#Override times anthropocene to correct for rounding 
#of to decades in update_timesteps.py 
dyear_ant = 5 #Resampled over 5 years
dyears_ant = np.cumsum([dyear_ant]*ds_ant.time.size)

timedeltas = [pd.Timedelta(dyear_ant*365, unit="D") for dyear_ant in dyears_ant]
times = [ds.time.values[-1] + timedelta for timedelta in timedeltas]

ds_ant = ds_ant.assign_coords(time=times)

#%%Calculate vols and masses
ds_ant["vol"] = np.abs(ds_ant.dx * ds_ant.dy * ds_ant.dz)
ds["vol"] = np.abs(ds.dx * ds.dy * ds.dz)

active = (ds["conc1"] > -1.).isel(time=-1)
active_ant = (ds_ant["conc1"] > -1.).isel(time=-1)
onshore = (ds["conc1"].x <= get_coastline(active))
onshore_ant = (ds_ant["conc1"].x <= get_coastline(active_ant))

mas, vol = get_mas_and_vol(ds, onshore, active)
mas_ant, vol_ant = get_mas_and_vol(ds_ant, onshore_ant, active_ant)

mas = xr.concat([mas, mas_ant], dim="time")
vol = xr.concat([vol, vol_ant], dim="time")

#%%Metadata for times
units = "days since 0000-01-01"
calendar = '365_day'

mas["time"].encoding["units"] = units
mas["time"].encoding["calendar"] = calendar
vol["time"].encoding["units"] = units
vol["time"].encoding["calendar"] = calendar

#%%Save
mas.to_netcdf(os.path.join(ant_fol, "mas_i{:03d}.nc".format(sim_nr)))
vol.to_netcdf(os.path.join(ant_fol, "vol_i{:03d}.nc".format(sim_nr)))