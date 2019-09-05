# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:52:04 2019

@author: engelen
"""

import xarray as xr
from glob import glob
import re, os, sys
import numpy as np
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
from imod import idf
import configparser

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def combine_all(ds_list):
    ds_tot = ds_list[0]
    for ds in ds_list[1:]:
        ds_tot = ds_tot.combine_first(ds)
    
    ds_tot = ds_tot.sel(y=slice(None, None, -1))    
    return(ds_tot)

def move_window_mean(da, dimmap):
    """Modified moving window that ignores nans of shifted window
    """
    return((da+da.shift(dimmap).fillna(0.))/2)

def calc_fresh_water_head(head, conc, dense_ref=1000., denselp=0.7143):
    """See Post, Kooi & Simmons (2007)
    Using hydraulic head measurements in variable-density ground water flow analyses
    """
    rho_i = dense_ref+conc*denselp
    return(rho_i/dense_ref * head - (rho_i - dense_ref)/dense_ref * head.z)

#%%Path management
if len(sys.argv) > 1:
    modelfol  = sys.argv[1]
    sim_nr = int(sys.argv[2])  
else:
    ##For Testing
    modelfol = r"g:\synthdelta\test_output\test_nspecies\SD_i123_nr00"
    mod_nr = 0

globpath=os.path.join(modelfol, "results", "results_{:03d}_*.nc".format(mod_nr))
run_path=os.path.join(modelfol, "*{}.run".format(mod_nr))

r = re.compile("([a-zA-Z]+_i[0-9]+_nr)([0-9]+)")
mname =  r.match(os.path.basename(modelfol)).group(1)

#%%Get z values from runfile
runf = glob(run_path)
if len(runf) != 1:
    raise ValueError("More than 1 runfile in this folder")

cfg = configparser.ConfigParser(delimiters=["="])
cfg.read(runf[0])
top=cfg.getfloat("dis", "top")
bots=[float(bot) for _, bot in cfg.items("dis") if "botm_" in _]
tb=np.array([top]+bots)
z=(tb[1:] + tb[:-1])/2
dz=tb[:-1] - tb[1:]

#%%Process
nc_paths = natural_sort(glob(globpath))

ds_list = [xr.broadcast(xr.open_dataset(
        nc, use_cftime=True, chunks = {"x" : -1, "y": -1, "layer": -1, "time": 1},
        drop_variables=["bdgsto", "bdgbnd"],
        ).assign(subdomain=i), 
                        exclude=["species", "time", "layer"])[0] for i, nc in enumerate(nc_paths)]
mids = [[np.mean(ds.x).values, np.mean(ds.y).values] for ds in ds_list]

ds_tot = combine_all(ds_list).transpose("species", "time", "layer", "y", "x")
ds_tot = ds_tot.compute()

subdomain = ds_tot["subdomain"]
ds_tot = ds_tot.drop(["subdomain"])

ds_tot = ds_tot.where(ds_tot["conc"].sel(species=1) < 1e20)

#%%Add z
ds_tot = ds_tot.assign_coords(z=xr.ones_like(ds_tot.layer) * z)
ds_tot = ds_tot.assign_coords(dz=xr.ones_like(ds_tot.layer) * dz)

#%%Calc velocities
ds_tot["vz"] = move_window_mean(ds_tot["bdgflf"], dimmap = {"layer" : 1})
ds_tot["vy"] = move_window_mean(ds_tot["bdgfff"], dimmap = {"y" : 1})
ds_tot["vx"] = move_window_mean(ds_tot["bdgfrf"], dimmap = {"x" : 1})
ds_tot = ds_tot.drop(["bdgflf", "bdgfff", "bdgfrf"])

ds_tot["vz"] = ds_tot["vz"] /(ds_tot["dx"] * ds_tot["dy"] * -1) #*-1 because dy is negative
ds_tot["vy"] = ds_tot["vy"] /(ds_tot["dx"] * ds_tot["dz"])
ds_tot["vx"] = ds_tot["vx"] /(ds_tot["dy"] * ds_tot["dz"]* -1) #*-1 because dy is negative

#%%Calc fresh water head
ds_tot["fhead"] = calc_fresh_water_head(ds_tot["head"], ds_tot["conc"].sel(species=1))
ds_tot = ds_tot.compute()

#%%Create initial heads and concentrations for next run.
#Load into memory as otherwise saving is really slow 
#since dask has to load all individual idfs into memory seperately.
ds_ini = ds_tot.isel(time=-1)[["conc", "head"]].load()

idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "bas", "head"), ds_ini["head"])
idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "btn", "conc"), ds_ini["conc"], 
         pattern = r"{name}_{time:%Y%m%d%H%M%S}_c{species}_l{layer}{extension}")

#%%Split species into seperate variables
species = ds_tot.species.values
for s in species:
    ds_tot["conc{}".format(s)] = ds_tot["conc"].sel(species=s)

ds_tot = ds_tot.drop("conc").drop("species")

#%%Save to netcdf
ds_tot = xr.where(np.isfinite(ds_tot["conc1"]), ds_tot, -9999.)
ds_tot["subdomain"] = subdomain

ds_tot = ds_tot.swap_dims({"layer" : "z"})

ds_tot = ds_tot.compute()
ds_tot.to_netcdf(os.path.join(globpath, "..", "results_{:03d}.nc".format(mod_nr)))

#%%Plot subdomains

qm = ds_tot["subdomain"].plot(cmap="prism", add_colorbar = False)

for i, mid in enumerate(mids):
    plt.text(mid[0], mid[1], str(i), horizontalalignment = 'center', verticalalignment = 'center')

plt.savefig(os.path.join(modelfol, "results", "subdomains.png"), dpi=200)