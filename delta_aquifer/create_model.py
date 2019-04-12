# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
from delta_aquifer import geometry, defaults
from delta_aquifer import boundary_conditions as bc
from collections import OrderedDict

#%%Path management
figfol = (
    r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Figures"
)
ncfol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Data"

#model_fol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Model"
model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"

#%%Parameters
# Morris parameters
lev = 4

# Geometry input parameters
L = (
    200000
)  # FIXED Check all delta lengths, so if this value is representative. Also check leakage factors whether this long enough
b = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas
a = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas

D = np.array([70, 250, 500, 1000])
dD = np.linspace(0.2, 0.6, num=lev)

alpha = np.linspace(0.75e-4, 1.5e-4, num=lev)  # Check this for deltas
beta = np.linspace(6e-4, 12e-4, num=lev)  # Check this for deltas
gamma = 5e-2  # FIXED will be corrected based on thickness.

# alpha = 8e-4 #in rads! For nile roughly: np.arctan(60/75000) = 8e-4 rad * 180/np.pi = 0.046 degrees
# beta = 1e-4 #in rads! For nile roughly: np.arctan(15/150000)

phi = np.linspace(0.125, 0.5, num=lev) * np.pi

clay_conf = np.linspace(0.2, 1.0, num=lev)
n_clay = np.linspace(0, 3, num=lev, dtype=int)
#SM = 0.3  # FIXED FOR NOW ELSE np.linspace(0.1, 0.6, num=4)
SM = np.linspace(0.1, 0.4, num=lev)

# Hydrogeological parameters
kh = np.logspace(0, 2, num=lev)
c_conf = np.logspace(0, 3, num=lev)
c_mar = np.logspace(1, 4, num=lev)
ani = np.logspace(0, 1.5, num=lev)

#Solute transport parameters
por = np.linspace(0.1, 0.35, num=lev)
al =  np.logspace(-0.6, 1, num=lev) #Maybe lower bound too low to get convergence

# Transgression
tra = np.linspace(0.25, 1, num=lev)
t_end = 0
t_start = 14
t_max = 7

# Model discretization
dx, dy = 1000, 1000
nz = 100
#nz=10
ts = (
    np.array(
        [
            50000,
            49000,
#            40000,
#            30000,
#            25000,
#            20000,
#            15000,
#            13000,
#            12000,
#            11000,
#            10000,
#            9000,
#            8000,
#            7000,
#            6000,
#            5000,
#            4000,
#            3000,
#            2000,
#            1000,
            0,
        ]
    )
    / 1000
)

#%%For testing
pars = {}

# Domain geometry
pars["a"] = a[1]
pars["b"] = b[1]
pars["D"] = D[1]
pars["dD"] = dD[1]
pars["alpha"] = alpha[2]
pars["beta"] = beta[2]
pars["gamma"] = gamma
pars["phi"] = phi[2]
pars["L"] = L

# Internal geometry
pars["SM"] = SM[2]
pars["clay_conf"] = clay_conf[3]
pars["n_clay"] = n_clay[2]

# Hydrogeological parameters
pars["kh"] = kh[1]
pars["c_conf"] = c_conf[2]
pars["c_mar"] = c_mar[2]
pars["ani"] = ani[2]
pars["bc-res"] = 100

#Solute transport
pars["por"] = por[-1]
pars["al"] = al[-2]

# Transgression
pars["tra"] = tra[2]
pars["t_end"] = t_end 
pars["t_start"] = t_start
pars["t_max"] = t_max 

# Discretization
pars["dx"], pars["dy"], pars["nz"] = dx, dy, nz
#%%Get geometry
geo = geometry.get_geometry(figfol=figfol, ncfol=ncfol, **pars)

#%%Create boundary conditions
# Path management
spratt = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\delta_aquifer\data\spratt2016.txt"

bcs = bc.boundary_conditions(spratt, ts, geo, figfol=figfol, ncfol=ncfol, **pars)

start_year = 2000 #Must be minimum 1900 for iMOD-SEAWAT

#%%
import imod
import xarray as xr
import os
import cftime

tb=bc._mid_to_binedges(geo["z"].values)

geo = geo.swap_dims({"z" : "layer"}).drop("z").sortby("y", ascending=False)
bcs = bcs.swap_dims({"z" : "layer"}).drop("z").sortby("y", ascending=False)

bcs = bcs.assign_coords(time = [cftime.DatetimeProlepticGregorian(t, 1, 1) for t in (bcs.time.values[::-1] * 100 + start_year)])

bcs["sea"] = bcs["sea"].where(bcs["sea"]==1)

geo = geo.sortby("layer")
bcs = bcs.sortby("layer")
tb = tb[::-1]

#%%

model = OrderedDict()

model["bnd"] = (geo["IBOUND"] * xr.full_like(bcs.time, 1)).astype(np.int16)
model["icbund"] = geo["IBOUND"]
model["khv"] = geo["Kh"]
model["kva"] = geo["Kh"]/pars["ani"]
#model["khv"] = xr.full_like(geo["IBOUND"], 10.0) * geo["IBOUND"]
#model["kva"] = xr.full_like(geo["IBOUND"], 10.0) * geo["IBOUND"]

model["top"] = xr.full_like(geo["IBOUND"], 1) * xr.DataArray(tb[:-1], coords = {"layer":geo.layer}, dims=["layer"])
model["bot"] = xr.full_like(geo["IBOUND"], 1) * xr.DataArray(tb[1:],  coords = {"layer":geo.layer}, dims=["layer"])
model["thickness"] = model["top"] - model["bot"]

shd = bcs["river_stage"].max(dim="layer").isel(time=0, drop=True).fillna(bcs["sea_level"].isel(time=0, drop=True)) * geo["IBOUND"]
model["shd"] = xr.where(geo["IBOUND"] == 1, shd, -9999)
model["sconc"] = xr.where(geo["IBOUND"]==1., 0.0, -9999.0)

model["ghb-head"] = xr.where(bcs["sea"]==1, bcs["sea_level"], np.nan)
model["ghb-cond"] = bcs["sea"] * pars["dx"] * pars["dy"] / pars["bc-res"]
model["ghb-dens"] = bcs["sea"] * 1025
model["ghb-conc"] = bcs["sea"] * 35.

model["riv-stage"] = bcs["river_stage"]
model["riv-cond"]  = xr.where(np.isfinite(bcs["river_stage"]), pars["dx"] * pars["dy"] / pars["bc-res"], bcs["river_stage"])
model["riv-bot"]   = bcs["river_stage"] - 10.
model["riv-dens"]  = xr.where(np.isfinite(bcs["river_stage"]), 1000., bcs["river_stage"])
model["riv-conc"]  = xr.where(np.isfinite(bcs["river_stage"]), 0., bcs["river_stage"])

run_pars = defaults.get_seawat_default_runfile()
run_pars["por"] = pars["por"]
run_pars["al"] = pars["al"]
run_pars["dt0"] = 10.0

imod.seawat_write(os.path.join(model_fol, "test_small"), model, runfile_parameters=run_pars)