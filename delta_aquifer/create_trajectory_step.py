# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
import pandas as pd
from delta_aquifer import Model
from delta_aquifer import initial_conditions as ic
import os, sys

from pkg_resources import resource_filename

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

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "sensitivity_analysis")))

spratt = os.path.join(datafol, "..", "spratt2016.txt")

#%%Parameters
fixpars = pd.read_csv(os.path.join(datafol, "fixed_pars.csv"), index_col=0).iloc[[0]]
varpars = pd.read_csv(os.path.join(datafol, "traj_real.csv" ), index_col=0).iloc[[sim_nr]]
pars = pd.concat([fixpars.reset_index(), varpars.reset_index()], axis=1)
pars.T.to_csv(os.path.join(ncfol, "parameters.csv"), sep="\t")
pars = dict([(key, pars[key].values[0]) for key in pars.columns])

# Time discretization
ts = (
    np.array(
        [
            46000,
            38000,
            30000,
            25000,
            20000,
            15000,
            13000,
            12000,
            11000,
            10000,
            9000,
            8000,
            7000,
            6000,
            5000,
            4000,
            3000,
            2000,
            1000,
            0,
        ]
    )
    / 1000
)

#%%Solver settings
hclose = 2e-4

#Set hclose higher for few troublesome ids
if sim_nr in [171, 173, 174, 175, 176, 177, 242, 281, 283, 284, 285, 286]:
    hclose=1e-3
#Rule of thumb for 3D MODFLOW models is dx*dy*hclose. Since SEAWAT expresses
#its fluxes in mass, RCLOSE has to be multiplied with the reference density. 
rclose = pars["dx"] * pars["dy"] * hclose * ic.c2dens(pars["C_f"])

#%%Create synthetic model
M = Model.Synthetic(pars, ts, hclose, rclose, figfol, ncfol, spratt)
M.prepare()
if len(sys.argv)==1:
    write_first_only=True
else:
    write_first_only=False

M.write_model(model_fol, mname, write_first_only=write_first_only)