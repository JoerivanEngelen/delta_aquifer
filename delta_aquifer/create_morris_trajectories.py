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
import matplotlib.pyplot as plt

#Morris parameters
lev = 4
grid_jump = 2

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

# Hydrogeological parameters
pars["kh"]          = np.logspace(0, 2, num=lev)
pars["kh_conf"]     = np.logspace(-3, 0, num=lev)
pars["kh_mar"]      = np.logspace(-4, -2, num=lev)
pars["ani"]         = np.logspace(0, 1.5, num=lev)
pars["bc-res"]      = 100
pars["riv_depth"]   = pars["D"][-1] #set this very low to linearize system

#Solute transport
pars["por"]         = np.linspace(0.1, 0.35, num=lev)
pars["al"]          = np.logspace(-0.6, 1, num=lev)

# Transgression
pars["tra"]         = np.linspace(0.25, 1, num=lev)
pars["t_start"], pars["t_max"], pars["t_end"] = 14, 7, 0

# Model discretization
pars["dx"], pars["dy"], pars["nz"] = 1000, 1000, 100

par_morris = [
        key for key, var in pars.items() if isinstance(var, np.ndarray) if len(var) == lev
        ]

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
fig2 = plt.figure()
sample_histograms(fig2, param_values, problem, {'color': 'y'})
plt.tight_layout()
plt.show()

#%%To 
traj_real = pd.DataFrame(OrderedDict([(par, pars[par][param_values[:, i]]) for i, par in enumerate(par_morris)]))
traj_ref  = pd.DataFrame(OrderedDict([(par, param_values[:, i]) for i, par in enumerate(par_morris)]))