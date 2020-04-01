# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 13:06:51 2020

@author: engelen
"""

import xarray as xr
from glob import glob
import os, sys
import pandas as pd
import re

def calc_ddim(da, dim):
    shifts = {dim : 1}
    ddim = da[dim] - da[dim].shift(**shifts)
    ddim[0] = ddim[1]
    return(ddim)

#%%Input control
if len(sys.argv)>1:
    res_fol = sys.argv[1]
else:
    res_fol = r"c:\Users\engelen\test_imodpython\test_results_30_deltas"

#%%Path management
paths = glob(os.path.join(res_fol, "synth_RD_i*_*", "input", "data", "bcs.nc"))
outfile = os.path.join(res_fol, "recharges.csv")
#%%Process
pattern = r"synth_RD_i(\d{3})_m24"
sim_nrs = [int(re.search(pattern, path).group(1)) for path in paths]

da_ls = [xr.open_dataset(path)["rch"].isel(time=-1) for path in paths]
cellsizes = [calc_ddim(da, "x") * calc_ddim(da, "y") for da in da_ls]
lith_ls = [xr.open_dataset(path)["lith"].isel(time=-1) for path in paths]
unconf_ls = [(lith != 2).min(dim="z") for lith in lith_ls]
da_ls = [(da * cellsize * unconf) for da, cellsize, unconf in zip(da_ls, cellsizes, unconf_ls)]

rch = [da.sum().values * 365.25 * 2 for da in da_ls] #Times two, because symmetrical domain cut in half

df = pd.DataFrame({"rch" : rch}, index=sim_nrs)

#%%Save
df.to_csv(outfile)