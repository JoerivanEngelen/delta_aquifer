# -*- coding: utf-8 -*-
"""
Created on Tue May 21 09:16:33 2019

@author: engelen
"""

import numpy as np
from SALib.sample.morris import sample
from SALib.plotting.morris import sample_histograms
from collections import OrderedDict
import pandas as pd
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()

import os
from pkg_resources import resource_filename

#%%Path management
traj_real_path = os.path.abspath(resource_filename("delta_aquifer", "../data/traj_real.csv"))
traj_id_path   = os.path.abspath(resource_filename("delta_aquifer", "../data/traj_id.csv"))

fig_out_path = os.path.abspath(resource_filename("delta_aquifer", "../data/morris_inp_dist.png"))

fixed_pars_path = os.path.abspath(resource_filename("delta_aquifer", "../data/fixed_pars.csv"))

#%%Morris parameters
lev = 4
grid_jump = 2

seed = 230

np.random.seed(seed=seed)

#%%Initial distribution
pars = OrderedDict()

# Domain geometry
pars["a"]           = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas
pars["b"]           = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas
pars["D"]           = np.logspace(np.log10(70), np.log10(1000), num=lev)
pars["dD"]          = np.linspace(0.2, 0.6, num=lev)
#For nile roughly: np.arctan(15/150000)
pars["alpha"]       = np.linspace(0.75e-4, 1.5e-4, num=lev)
#For nile roughly: np.arctan(60/75000) = 8e-4 rad * 180/np.pi = 0.046 degrees
pars["beta"]        = np.linspace(6e-4, 12e-4, num=lev)
pars["gamma"]       = 5e-2 # FIXED will be corrected based on thickness.
pars["phi"]         = np.linspace(0.125, 0.5, num=lev) * np.pi
pars["L"]           = 200000
# FIXED Check all delta lengths, so if this value is representative. 
#Also check leakage factors whether this long enough

# Internal geometry
pars["SM"]          = np.linspace(0.1, 0.4, num=lev)
pars["clay_conf"]   = np.linspace(0.2, 1.0, num=lev)
pars["n_clay"]      = np.linspace(0, 3, num=lev, dtype=int)
pars["N_pal"]       = np.linspace(1, 4, num=lev, dtype=int)
pars["s_pal"]       = np.linspace(0, 1.0, num=lev)

# Hydrogeological parameters
pars["kh"]          = np.logspace(0, 2, num=lev)
pars["kh_conf"]     = np.logspace(-3, 0, num=lev)
pars["kh_mar"]      = np.logspace(-4, -2, num=lev)
pars["f_kh_pal"]    = np.linspace(0, 1, num=lev)
pars["ani"]         = 10.

#River system
pars["N_chan"]      = np.linspace(1, 4, num=lev, dtype=int)
pars["f_cond_chan"] = np.logspace(0, 2, num=lev)
pars["intrusion_L"] = np.linspace(0, 0.5, num=lev)
pars["bc_res"]      = 100.

#Recharge
pars["rch_rate"]    = np.linspace(0.0, 0.4, num=lev) #Deze nog checken met Perry/Joost

#Solute transport
pars["por"]         = np.linspace(0.1, 0.35, num=lev)
pars["al"]          = np.logspace(-0.6, 1, num=lev)
pars["trpt"]        = 0.1
pars["trpv"]        = 0.01
pars["diff"]        = 8.64e-5
pars["c_f"]         = 0.
pars["c_s"]         = 35.

# Transgression
pars["tra"]         = np.linspace(0.25, 1, num=lev)
pars["t_start"], pars["t_max"], pars["t_end"] = 14, 7, 0

# Model discretization
pars["dx"], pars["dy"], pars["nz"] = 1000, 1000, 100

par_morris = [
        key for key, var in pars.items() if isinstance(var, np.ndarray) if len(var) == lev
        ]

fixed_pars = OrderedDict([(key, var) for key, var in pars.items() if not isinstance(var, np.ndarray)])

#%%Create morris sample
problem = {
        "num_vars" : len(par_morris),
        "names" : par_morris,
        "groups" : None,
        "bounds" : [[0, lev-1]] * len(par_morris)
        }

param_values = sample(problem, N=12, grid_jump=grid_jump, num_levels=lev,
                      sample4uniformity = 1000).astype(np.int64)

#%%Plot
fig = plt.figure(figsize=(8, 6))
sample_histograms(fig, param_values, problem, {'color': 'y'})
plt.tight_layout()
plt.savefig(fig_out_path, dpi=100)

#%%Save as csv
traj_real = pd.DataFrame(OrderedDict([(par, pars[par][param_values[:, i]]) for i, par in enumerate(par_morris)]))
traj_id  = pd.DataFrame(OrderedDict([(par, param_values[:, i]) for i, par in enumerate(par_morris)]))
fixed_pars = pd.DataFrame(fixed_pars,  index=["fix"])

traj_real.to_csv(traj_real_path)
traj_id.to_csv(traj_id_path)
fixed_pars.to_csv(fixed_pars_path)
