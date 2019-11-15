# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 09:42:57 2019

@author: engelen
"""

import xarray as xr
import numpy as np
from imod import idf
import os

def get_conc(bcs):
    sea = (bcs["sea"]==1)
    bcs["conc"]  = xr.where(sea, bcs["sea_conc"] , bcs["riv_conc"])
    return(bcs)

def compare(path_old, path_new, idfs, pattern=None):
    da_old = idf.open(os.path.join(path_old, idfs), pattern=pattern)
    da_new = idf.open(os.path.join(path_new, idfs), pattern=pattern)
    equal=da_old.equals(da_new)
    print("%s is equal: %r" % (idfs, equal))

#%%Path management
path_old = r"c:\Users\engelen\test_imodpython\synth_delta_test\SD_i239_master\SD_i239_nr00"
path_new = r"c:\Users\engelen\test_imodpython\synth_delta_test\SD_i239_refactor\SD_i239_nr00"

#%%Check boundary conditions delta_aquifer
bcs_old = xr.open_dataset(os.path.join(path_old, "..", "input", "data", "bcs.nc"))
bcs_new = xr.open_dataset(os.path.join(path_new, "..", "input", "data", "bcs.nc"))

for var in list(bcs_old.keys()):
    flag = bcs_old[var].equals(bcs_new[var])
    print("%s: %r" % (var, flag))

#%%check IBOUND
compare(path_old, path_new, os.path.join("bas", "ibound_*.idf"))
compare(path_old, path_new, os.path.join("btn", "icbund_*.idf"))

#%%Check initial conditions
compare(path_old, path_new, os.path.join("btn", "starting_concentration_*.idf"), pattern=r"{name}_c{species}_l{layer}")
compare(path_old, path_new, os.path.join("bas", "starting_head_*.idf"))

#%%Check boundary conditions
compare(path_old, path_new, os.path.join("ghb", "concentration_*.idf"), pattern=r"{name}_c{species}_l{layer}")
compare(path_old, path_new, os.path.join("ghb", "conductance_*.idf"))
compare(path_old, path_new, os.path.join("ghb", "head_*.idf"))
compare(path_old, path_new, os.path.join("ghb", "density_*.idf"))

compare(path_old, path_new, os.path.join("rch", "rate_*.idf"))

#%%Check k values
compare(path_old, path_new, os.path.join("lpf", "k_horizontal_*.idf"))
compare(path_old, path_new, os.path.join("lpf", "k_vertical_*.idf"))

#%%Check runfiles
with open(os.path.join(path_old, "SD_i239_nr00.run"), mode="r") as r_old:
    lines_old = r_old.readlines()

with open(os.path.join(path_new, "SD_i239_nr00.run"), mode="r") as r_new:
    lines_new = r_new.readlines()

diff = [ x for x in lines_old if x not in lines_new]
print("Runfiles equal?")
print(len(diff)==0)

