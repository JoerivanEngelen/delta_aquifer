# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
from . import geometry

#%%Path management
figfol=r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Figures"
netcdf=r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Data\synth_model.nc"

#%%Parameters
#Morris parameters
lev = 4

#Geometry input parameters
L = 200000 #FIXED Check all delta lengths, so if this value is representative. Also check leakage factors whether this long enough
b = np.linspace(0.3, 0.6, num=lev) #Check this for deltas
a = np.linspace(0.3, 0.6, num=lev) #Check this for deltas

D = np.array([70, 250, 500, 1000])
dD = np.linspace(0.2, 0.6, num=lev)

alpha  = np.linspace(0.75e-4, 1.5e-4, num=lev) #Check this for deltas
beta = np.linspace(6e-4, 12e-4, num=lev) #Check this for deltas
gamma = 5e-2 #FIXED will be corrected based on thickness.

#alpha = 8e-4 #in rads! For nile roughly: np.arctan(60/75000) = 8e-4 rad * 180/np.pi = 0.046 degrees
#beta = 1e-4 #in rads! For nile roughly: np.arctan(15/150000)

phi = np.linspace(0.125, 0.5, num=lev) * np.pi

clay_conf = np.linspace(0.2, 1., num=lev)
n_clay = np.linspace(0, 3, num=lev, dtype=int)
SM = 0.3 #FIXED FOR NOW ELSE np.linspace(0.1, 0.6, num=4)

#Model discretization
dx, dy = 1000, 1000
nz = 100

#Hydrogeological parameters
kh = np.logspace(0, 2, num=lev)
c_conf = np.logspace(0, 3, num=lev)
c_mar = np.logspace(1, 4, num=lev)
ani = np.logspace(0, 1.5, num=lev)

#%%For testing
pars = {}

#Domain geometry
pars["a"]  = a[1]
pars["b"]  = b[1]
pars["D"]  = D[1]
pars["dD"] = dD[0]
pars["alpha"] = alpha[2]
pars["beta"]  = beta[2]
pars["gamma"] = gamma
pars["phi"] = phi[2]
pars["L"]   = L
pars["SM"]  = SM

#Internal geometry
pars["clay_conf"] = clay_conf[2]
pars["n_clay"] = n_clay[3]

#Hydrogeological parameters
pars["kh"] = kh[1]
pars["c_conf"] = c_conf[2]
pars["c_mar"]  = c_mar[2]
pars["ani"]    = ani[2]

#Discretization
pars["dx"], pars["dy"], pars["nz"] = dx, dy, nz

#%%
d3, _ = geometry.get_geometry(figfol=figfol, netcdf=netcdf, **pars)