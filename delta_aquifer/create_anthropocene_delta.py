# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:44:52 2020

@author: engelen
"""

import numpy as np
import pandas as pd
from delta_aquifer import Model
from delta_aquifer import initial_conditions as ic
import os, sys
from glob import glob

from pkg_resources import resource_filename


par_names=['gamma', 'L', 'S_s', 'bc_res', 't_start', 't_end', 'trpt', 'trpv', 
           'D_m', 'C_f', 'C_s', 'dx', 'dy', 'nz', 'l_a', 'H_b', 'f_H', 'alpha', 
           'beta', 'phi_f', 'f_aqt', 'l_conf', 'N_aqt', 'N_pal', 's_pal', 'Kh_aqf', 
           'Kv_aqt', 'f_Kh_pal', 'Kh_Kv', 'N_chan', 'f_chan', 'l_surf_end', 'l_tra', 
           't_tra', 'R', 'n', 'a_l', "riv_res"]

#Set constants
fixed_pars = {"t_start" : 12,
              "t_end" : 0,
              "trpt" : 0.1,
              "trpv" : 0.01,
              "D_m" : 8.64e-5,
              "C_f" : 0.0,
              "C_s" : 35.0,
              "nz" : 100,
              "f_chan" : 1.0,
              "bc_res" : 10., 
              "S_s" : 6.31e-5,
              "a_l" : 2.}
                            

#%%Path management
if len(sys.argv) > 1: #For local, interactive testing purposes.
    model_fol  = sys.argv[1]
    sim_nr = int(sys.argv[2])
    init_path = sys.argv[3]
    write_first_only=False
else:
    #Local testing on my own windows laptop
    model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"
    sim_nr = 153
    init_path = r"g:\test_UCN\results\synth_RD_i153_m24_7666720"
    write_first_only=True
    
mname = "RD_i{:03d}".format(sim_nr)

figfol = os.path.join(model_fol, mname, "input", "figures")
ncfol  = os.path.join(model_fol, mname, "input", "data")

os.makedirs(figfol, exist_ok=True)
os.makedirs(ncfol,  exist_ok=True)

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))

spratt = os.path.join(datafol, "..", "spratt2016.txt")
abstraction_f = os.path.join(datafol, "abstractions", "{}.nc")

init_f = glob(os.path.join(init_path, "results_[0-9][0-9][0-9].nc"))
init_f.sort()
init_f = init_f[-1]

#%%Read inputs
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Attach constants
fixed_pars = pd.DataFrame(fixed_pars, index=par.index)
par = pd.concat([par, fixed_pars], axis=1)

#%%Select simulation parameters
sim_par = par.loc[sim_nr]
#%%Set ts
ts = np.arange(55, -1, -5)
ts = ts / 1000.
sim_par.at["t_tra"] = sim_par["t_tra"] / 1000 

#%%Groundwater exstractions last 50 years
abstraction_path = abstraction_f.format(sim_par["Delta"])

#%%Read inits
import xarray as xr

def _rev_dims(da, *dims):
    """Reverse dims, 
    alternative to all the slow sortby actions that slowed this package down previously
    """
    
    kwargs = dict([(dim, da[dim][::-1]) for dim in dims])
    return(da.reindex(**kwargs))

def _mirror_dims(da, *dims):
    kwargs = dict([(dim , da[dim]*-1) for dim in dims])
    mirror = da.assign_coords(**kwargs)
    return(_rev_dims(mirror, *dims))

inits = xr.open_dataset(init_f).isel(time=-1)
inits = xr.concat([inits, _mirror_dims(inits, "y")], dim="y")

#%%Solver settings
hclose = 2e-4

#Rule of thumb for 3D MODFLOW models is dx*dy*hclose. Since SEAWAT expresses
#its fluxes in mass, RCLOSE has to be multiplied with the reference density. 
rclose = sim_par["dx"] * sim_par["dy"] * hclose * ic.c2dens(sim_par["C_f"])

#%%Create synthetic model
M = Model.Synthetic(sim_par.to_dict(), ts, hclose, rclose, figfol, ncfol, 
                    spratt, abstraction_path=abstraction_path, init_salt = init_f, 
                    init_half2full=True, half_model=False)

#%%Depth 
import xarray as xr

lith = M.geo["aqt"]

aqf_nrs = xr.full_like(lith, 0, dtype=np.int64) #This step not necessary if this issue gets fixed https://github.com/pydata/xarray/issues/2017
axis = dict(axis = aqf_nrs.get_axis_num("z"))
aqf_nrs.values = np.flip(np.maximum.accumulate(np.flip(lith.values, **axis), **axis), **axis) #Flipping to do a reverse accumulate
aqf_nrs = xr.where(M.geo["Kh"].isel(time=-1)==(sim_par["Kh_aqf"]), aqf_nrs, np.nan)

onshore = M.geo["topsys"]>1

#%%Continue processing model

#M.prepare(init_salt = inits, half_model=False)   

#M.write_model(model_fol, mname, write_first_only=write_first_only)

#%%Test
w = M.wel
w1 = w.loc[w["time"] == 0.0025]["Q"].sum()*365.25/(sim_par["dx"]*sim_par["dy"])

import xarray as xr
a = xr.open_dataset(abstraction_path)["__xarray_dataarray_variable__"]
onshore = (np.isnan(M.bcs.sea.max(dim="z")) & (M.geo["IBOUND"]==1.))
a1 = a.sel(time=2014).where(onshore).sum()