# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:24:25 2019

@author: engelen
"""

import imod
import xarray as xr
import sys, os

modelfol = sys.arg[1]
mod_nr = sys.argv[2]
s_nr = sys.argv[3]
#modelfol = r"g:\synthdelta\test_idf_output"
#mod_nr = 1
#s_nr = 17

#%%Path management
mname = os.path.basename(modelfol)
path = os.path.join(modelfol, mname+str(mod_nr), "results" ,"*", "*_p{:03d}.IDF".format(s_nr))

#%%Process
ds = imod.idf.open_dataset(path, use_cftime=True, pattern=r"{name}_{time}_l{layer}_p\d*")
nc_path = os.path.join(modelfol, "results", "results_{:03d}_p{:03d}.nc".format(mod_nr, s_nr))
xr.merge(ds.values()).to_netcdf(nc_path)