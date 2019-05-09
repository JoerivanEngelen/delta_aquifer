# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:24:25 2019

@author: engelen
"""

import imod
import xarray as xr
import sys, os

modelfol = sys.argv[1]
mod_nr = int(sys.argv[2])
s_nr = int(sys.argv[3])
##For testing
#modelfol = r"g:\synthdelta\test_idf_output2"
#mod_nr = 1
#s_nr = 0

#%%Path management
path = os.path.join(modelfol, "results" ,"*", "*_p{:03d}.[iI][dD][fF]".format(s_nr))

#%%Process
ds = imod.idf.open_dataset(path, use_cftime=True, pattern=r"{name}_{time}_l{layer}_p\d*")
nc_path = os.path.join(modelfol, "results", "results_{:03d}_p{:03d}.nc".format(mod_nr, s_nr))
ds = xr.merge(ds.values())
ds.time.encoding["units"] = r"days since 1999-01-01"
ds.time.encoding["calendar"] = r"proleptic_gregorian"
ds.to_netcdf(nc_path, unlimited_dims=["time"])