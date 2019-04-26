# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
from delta_aquifer import geometry, defaults
from delta_aquifer import boundary_conditions as bc
from delta_aquifer import non_convergence as ncg
from delta_aquifer import initial_conditions as ic
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
kh_conf = np.logspace(-3, 0, num=lev)
kh_mar = np.logspace(-4, -2, num=lev)
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
pars["n_clay"] = n_clay[1]

# Hydrogeological parameters
pars["kh"] = kh[1]
pars["kh_conf"] = kh_conf[1]
pars["kh_mar"] = kh_mar[1]
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
bcs["sea"] = bcs["sea"].where(bcs["sea"]==1)

geo = geo.sel(y=slice(0, geo.y.max()))
bcs = bcs.sel(y=slice(0, geo.y.max()))

#%%
import imod
import xarray as xr
import os
import cftime
    
start_year = 2000 #Must be minimum 1900 for iMOD-SEAWAT
approx_init = True

c_f = 0.0
c_s = 35.

rho_f, rho_s = ic.c2rho(c_f), ic.c2rho(c_s)

tb=bc._mid_to_binedges(geo["z"].values)[::-1]
z=geo.z

shd, sconc = ic.get_ic(bcs, geo, c_f, c_s, approx_init=approx_init)

geo = geo.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
bcs = bcs.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
sconc = sconc.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
shd = shd.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)

#bcs = bcs.assign_coords(time = [cftime.DatetimeProlepticGregorian(t, 1, 1) for t in (bcs.time.values[::-1] * 100 + start_year)])
bcs = bcs.assign_coords(time = [np.datetime64("2000-01-01"), np.datetime64("2001-01-01")])

#%%

m = imod.wq.SeawatModel("test_half_new")
m["bas"] = imod.wq.BasicFlow(ibound=geo["IBOUND"].assign_coords(dx=dx, dy=-dy), 
                             top=tb[0], 
                             bottom=xr.DataArray(tb[1:], {"layer": geo.layer}, ("layer")), 
                             starting_head=xr.where(geo["IBOUND"]==1,  shd, -9999.0))

m["lpf"] = imod.wq.LayerPropertyFlow(
    k_horizontal=geo["Kh"], k_vertical=geo["Kh"]/pars["ani"], specific_storage=0.0
)


m["btn"] = imod.wq.BasicTransport(
    icbund=geo["IBOUND"], 
    starting_concentration=xr.where(geo["IBOUND"]==1., sconc, -9999.0), 
    porosity=pars["por"]
)

m["adv"] = imod.wq.AdvectionTVD(courant=0.9)
m["dsp"] = imod.wq.Dispersion(
        longitudinal=pars["al"], 
        ratio_horizontal=0.1,
        ratio_vertical=0.01,
        diffusion_coefficient=8.64e-5
)

m["vdf"] = imod.wq.VariableDensityFlow(density_concentration_slope=0.7143)

m["ghb"] = imod.wq.GeneralHeadBoundary(head = xr.where(bcs["sea"]==1, bcs["sea_level"], np.nan),
                                       conductance=bcs["sea"] * pars["dx"] * pars["dy"] / pars["bc-res"],
                                       density=bcs["sea"] * rho_s,
                                       concentration=bcs["sea"] * c_s)

m["riv"] = imod.wq.River(stage = bcs["river_stage"],
                         conductance = xr.where(np.isfinite(bcs["river_stage"]), pars["dx"] * pars["dy"] / pars["bc-res"], np.nan),
                         bottom_elevation = bcs["river_stage"] - 10.,
                         density = xr.where(np.isfinite(bcs["river_stage"]), rho_f, np.nan),
                         concentration = xr.where(np.isfinite(bcs["river_stage"]), 0., np.nan))

m["pksf"] = imod.wq.ParallelKrylovFlowSolver(1000, 100, 0.0001, 100., 0.98,
                                             partition="uniform",
                                             solver="pcg",
                                             preconditioner="ilu",
                                             deflate=False,
                                             debug=False,)

m["pkst"] = imod.wq.ParallelKrylovTransportSolver(1000, 30, 
                                             cclose=1e-6,
                                             relax=0.98,
                                             partition="uniform",
                                             solver="bicgstab",
                                             preconditioner="ilu",
                                             debug=False,)

m["oc"] = imod.wq.OutputControl(save_head_idf=True, save_concentration_idf=True)

m.time_discretization(endtime="2005-01-01")
m.write(directory = r"c:\Users\engelen\test_imodpython\synth_delta_test\test_obj")

#%%
pointer_grid = geo["IBOUND"].sum(dim="layer")

with open(os.path.join(model_fol, "test_half", "pointer_grid.asc"), "w") as pg:
    for row in pointer_grid:
        pg.write(("{:8.2f}  " * (len(row)-1) + "{:8.2f}\n").format(*tuple(row.values)))
    

##This does not work because the pointer grid should be .ASC, not .IDF (with a header)
#imod.idf.save(os.path.join(model_fol, "test_half", "pointer_grid.idf"), geo["IBOUND"].sum(dim="layer"))

#%%non_conv_analyser
#cell1 = (11, 177, 102)
#ncg1, xyz1 = ncg.look_around(model, cell1, n=2, var=["ghb-head", "riv-stage", "khv", "icbund"])
#
#cell2 = (17, 193, 128)
#ncg2, xyz2 = ncg.look_around(model, cell2, n=2, var=["ghb-head", "riv-stage", "khv", "icbund"])
