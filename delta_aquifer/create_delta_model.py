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
           't_tra', 'R', 'n', 'a_l']


fixed_pars = pd.DataFrame({"t_start" : 12,
                           "t_end" : 0,
                           "trpt" : 0.1,
                           "trpv" : 0.01,
                           "D_m" : 8.64e-5,
                           "C_f" : 0.0,
                           "C_s" : 35.0,
                           "nz" : 100,
                           "f_chan" : 1.0})

#%%Path management
if len(sys.argv) > 1:
    model_fol  = sys.argv[1]
    sim_nr = int(sys.argv[2])
else:
    #Local testing on my own windows laptop
    model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"
    sim_nr = 242
    
mname = "SD_i{:03d}".format(sim_nr)

figfol = os.path.join(model_fol, mname, "input", "figures")
ncfol  = os.path.join(model_fol, mname, "input", "data")

os.makedirs(figfol, exist_ok=True)
os.makedirs(ncfol,  exist_ok=True)

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))

spratt = os.path.join(datafol, "spratt2016.txt")

#%%Read inputs
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Clip off all ofshore >200 km
L_b = par["L"]*(1-par["l_a"])
L_a = par["L"]*par["l_a"]
L_b = L_b.clip(upper=200.)

par["L"] = L_a + L_b

#%%Determine discretization
par["dx"] = 1000
par["dx"] = par["dx"].where(par["L"] > 130., 500)
par["dx"] = par["dx"].where(par["L"] < 400., 2000)

par["dy"] = par["dx"]

#%%Set ts
ts = np.array([30000, 25000, 20000, 15000, 13000, 12000, 11000, 10000, 9000,
               8000, 7000, 6000, 5000,  4000, 3000, 2000, 1000, 0,])

ts = np.concatenate((np.arange(92000, 30000, -8000), ts)) / 1000.
    
#%%TODO
#Start initially completely salt
#ts 100k year



##%%Solver settings
#hclose = 2e-4
#
##Rule of thumb for 3D MODFLOW models is dx*dy*hclose. Since SEAWAT expresses
##its fluxes in mass, RCLOSE has to be multiplied with the reference density. 
#rclose = pars["dx"] * pars["dy"] * hclose * ic.c2dens(pars["C_f"])

##%%Create synthetic model
#M = Model.Synthetic(pars, ts, hclose, rclose, figfol, ncfol, spratt)
#M.prepare()
#if len(sys.argv)==1:
#    write_first_only=True
#else:
#    write_first_only=False
#
#M.write_model(model_fol, mname, write_first_only=write_first_only)