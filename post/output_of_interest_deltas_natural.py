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

def get_mas_and_vol(ds, coastline, active):
    onshore = (ds["conc1"].x <= coastline)
    fw_vols  = xr.where(active & (ds["conc1"] < 5.),
                        ds["vol"], 0)
   
    #Calculate distance from contemporary coastline, positive is intrusion onshore
    distance_from_coastline = (ds.x - coastline.drop(["time", "z", "layer"])) * -1
    intrusion_length = xr.where(active & (ds["conc1"] > 17.5),
                        distance_from_coastline, np.nan)
    
    fw_vols_offshore = fw_vols.where(~onshore, 0)
    
    tot_vols = xr.where(active, ds["vol"], 0)
    tot_offshore_vols = (~onshore) * tot_vols
    
    sal_masses = xr.where(
            (ds["conc1"] > 0.), 
            ds["conc1"]*ds["vol"], 0
            )
    
    ol_masses = sal_masses*ds["conc2"]
    
    mas={}
    mas["x_sal"]       = intrusion_length.max(dim=["x", "y", "z"])
    mas["sal"]         = spat_sum(sal_masses, "sal")
    mas["ol_sal"]      = spat_sum(ol_masses, "ol_sal")
    mas["sal_onshore"] = spat_sum(xr.where(onshore, 
                           sal_masses, 0), "sal_onshore")
    mas["ol_sal_onshore"]  = spat_sum(xr.where(onshore, 
                               ol_masses, 0), "ol_sal_onshore")
    mas["sal_offshore"]    = mas["sal"] - mas["sal_onshore"]
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

    #Metadata for times
    units = "days since 0000-01-01"
    calendar = '365_day'
    
    mas["time"].encoding["units"] = units
    mas["time"].encoding["calendar"] = calendar
    vol["time"].encoding["units"] = units
    vol["time"].encoding["calendar"] = calendar

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

nat_results=glob(os.path.join(nat_fol, "results_[0-9][0-9][0-9].nc"))

#%%Read
ds = xr.open_mfdataset(nat_results, 
                       chunks={"time" : 10, "z": -1, "y": -1, "x": -1},
                       combine="nested",
                       concat_dim="time")

#%%Calculate vols and masses
ds["vol"] = np.abs(ds.dx * ds.dy * ds.dz)

active = (ds["conc1"] > -1.).isel(time=-1)
coastline = get_coastline(active)

mas, vol = get_mas_and_vol(ds, coastline, active)
mas = mas * 2
mas["x_sal"] = mas["x_sal"]/2
vol = vol * 2

mas_pump, vol_pump = get_mas_and_vol(ds.sel(z=slice(None, -300)), 
                                     coastline, active.sel(z=slice(None, -300)))
mas_pump = mas_pump * 2
mas_pump["x_sal"] = mas_pump["x_sal"]/2
vol_pump = vol_pump * 2

#%%Save
mas.to_netcdf(os.path.join(nat_fol, "mas_i{:03d}.nc".format(sim_nr)))
vol.to_netcdf(os.path.join(nat_fol, "vol_i{:03d}.nc".format(sim_nr)))
mas_pump.to_netcdf(os.path.join(nat_fol, "mas_pump_i{:03d}.nc".format(sim_nr)))
vol_pump.to_netcdf(os.path.join(nat_fol, "vol_pump_i{:03d}.nc".format(sim_nr)))