# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:52:04 2019

@author: engelen
"""

import xarray as xr
from glob import glob
import re, os, sys
import numpy as np
import matplotlib.pyplot as plt
from imod import idf

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def combine_all(ds_list):
    ds_tot = ds_list[0]
    for ds in ds_list[1:]:
        ds_tot = ds_tot.combine_first(ds)
    
    return(ds_tot)

#%%Path management
##For Testing
#modelfol = r"g:\synthdelta\test_idf_output"
#mod_nr = 2

modelfol = sys.argv[1]
mod_nr = sys.argv[2]

globpath=os.path.join(modelfol, "results", "results_{:03d}_*.nc".format(mod_nr))

r = re.compile("([a-zA-Z]+)([0-9]+)")
mname =  r.match(os.path.basename(modelfol)).group(1)

#%%Process
nc_paths = natural_sort(glob(globpath))

ds_list = [xr.broadcast(xr.open_dataset(
        nc, use_cftime=True, chunks = {"x" : -1, "y": -1, "layer": -1, "time": 1}
        ).assign(subdomain=i), 
                        exclude=["time", "layer"])[0] for i, nc in enumerate(nc_paths)]
mids = [[np.mean(ds.x).values, np.mean(ds.y).values] for ds in ds_list]

ds_tot = combine_all(ds_list).transpose("time", "layer", "y", "x")
ds_tot.to_netcdf(os.path.join(globpath, "..", "results_{:03d}.nc".format(mod_nr)))

#%%Create initial heads and concentrations for next run.
#Load into memory as otherwise saving is really slow 
#since dask has to load all individual idfs into memory seperately.
ds_ini = ds_tot.isel(time=-1)[["conc", "head"]].load()

idf.save(os.path.join(modelfol, "..", mname+str(mod_nr+1), "bas", "head"), ds_ini["head"], use_cftime=True)
idf.save(os.path.join(modelfol, "..", mname+str(mod_nr+1), "btn", "conc"), ds_ini["conc"], use_cftime=True)

#%%Plot subdomains
qm = ds_tot["subdomain"].plot(cmap="prism", add_colorbar = False)

for i, mid in enumerate(mids):
    plt.text(mid[0], mid[1], str(i), horizontalalignment = 'center', verticalalignment = 'center')
    
plt.savefig(os.path.join(globpath, "..", "subdomains.png"), dpi=200)