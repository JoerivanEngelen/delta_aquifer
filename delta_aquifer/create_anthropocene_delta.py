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
    write_first_only=False
else:
    #Local testing on my own windows laptop
    model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"
    sim_nr = 80
    write_first_only=True
    
mname = "RD_i{:03d}".format(sim_nr)

figfol = os.path.join(model_fol, mname, "input", "figures")
ncfol  = os.path.join(model_fol, mname, "input", "data")

os.makedirs(figfol, exist_ok=True)
os.makedirs(ncfol,  exist_ok=True)

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))

spratt = os.path.join(datafol, "..", "spratt2016.txt")
abstraction_f = os.path.join(datafol, "abstractions", "{}.nc")

#%%Read inputs
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Attach constants
fixed_pars = pd.DataFrame(fixed_pars, index=par.index)
par = pd.concat([par, fixed_pars], axis=1)

#%%Select simulation parameters
sim_par = par.loc[sim_nr]
#%%Set ts
ts = np.array([30000, 25000, 20000, 15000, 13000, 12000, 11000, 10000, 9000,
               8000, 7000, 6000, 5000,  4000, 3000, 2000, 1000])

ts = np.concatenate((np.arange(125000, 30000, -8000), ts))  #Add Pleistocene
anthro_years = np.arange(55, -1, -5)

ts = np.concatenate((ts, anthro_years))                     #Add Anthropocene
ts = ts / 1000.

#%%Groundwater exstractions last 50 years
abstraction_path = abstraction_f.format(sim_par["Delta"])

#%%Solver settings
hclose = 2e-4

#Rule of thumb for 3D MODFLOW models is dx*dy*hclose. Since SEAWAT expresses
#its fluxes in mass, RCLOSE has to be multiplied with the reference density. 
rclose = sim_par["dx"] * sim_par["dy"] * hclose * ic.c2dens(sim_par["C_f"])

#%%Create synthetic model
M = Model.Synthetic(sim_par.to_dict(), ts, hclose, rclose, figfol, ncfol, 
                    spratt, abstraction_path=abstraction_path)
#M.prepare(init_salt = True, half_model=True)   

#M.write_model(model_fol, mname, write_first_only=write_first_only)

#%%Test
w = M.wel
w1 = w.loc[w["time"] == 0.0025]["Q"].sum()*365.25/(sim_par["dx"]*sim_par["dy"])

import xarray as xr
a = xr.open_dataset(abstraction_path)["__xarray_dataarray_variable__"]
onshore = (np.isnan(M.bcs.sea.max(dim="z")) & (M.geo["IBOUND"]==1.))
a1 = a.sel(time=2014).where(onshore).sum()