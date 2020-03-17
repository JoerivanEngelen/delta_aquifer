# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 13:12:40 2020

@author: engelen
"""

from imod import idf
import xarray as xr
import os, sys
import re

#%%Path management
if len(sys.argv) > 1:
    modelfol  = sys.argv[1]
    mod_nr = int(sys.argv[2])
else:
    #Local testing on my own windows laptop
    modelfol = r"g:\test_UCN\RD_i153_nr00"
    mod_nr = 0

#%%Path management
path = os.path.join(modelfol, "results", "results_{:03d}.nc".format(mod_nr))

r = re.compile("([a-zA-Z]+_i[0-9]+_nr)([0-9]+)")
mname =  r.match(os.path.basename(modelfol)).group(1)

#%%Read
ds = xr.open_dataset(path, use_cftime=True)

ds_ini = ds.isel(time=-1)

#%%Format concentrations
ds_ini = ds_ini.rename(name_dict={"conc1" : "conc"})
ds_ini["conc"] = ds_ini["conc"].expand_dims(species=[1,2])
ds_ini["conc"] = xr.where(ds_ini["conc"].species==2, 
      ds_ini["conc2"], ds_ini["conc"])
ds_ini = ds_ini.drop(labels="conc2")

#%%

idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "bas", "head"), ds_ini["head"])
idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "btn", "conc"), ds_ini["conc"], 
         pattern = r"{name}_{time:%Y%m%d%H%M%S}_c{species}_l{layer}{extension}")
