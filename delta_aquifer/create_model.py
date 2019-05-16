# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
from delta_aquifer import geometry, defaults, time_util
from delta_aquifer import boundary_conditions as bc
from delta_aquifer import non_convergence as ncg
from delta_aquifer import initial_conditions as ic
import os

import imod
import xarray as xr
import cftime

#%%Path management
figfol = (
    r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Figures"
)
ncfol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Data"

#model_fol = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Model"
model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"

mname = "test_dynamic_conf"

#%%Parameters
# Morris parameters
lev = 4

# Geometry input parameters
L = (
    200000
)  # FIXED Check all delta lengths, so if this value is representative. Also check leakage factors whether this long enough
b = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas
a = np.linspace(0.3, 0.6, num=lev)  # Check this for deltas

#D = np.array([70, 250, 500, 1000])
D = np.logspace(np.log10(70), np.log10(1000), num=lev)
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
#pars["riv_depth"] = 10.
pars["riv_depth"] = pars["D"] #set this very low to see if rivers are causing non-convergence

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

topbot=bc._mid_to_binedges(geo["z"].values)[::-1]

#%%Create boundary conditions
# Path management
spratt = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\delta_aquifer\data\spratt2016.txt"

bcs = bc.boundary_conditions(spratt, ts, geo, figfol=figfol, ncfol=ncfol, **pars)
bcs["sea"] = bcs["sea"].where(bcs["sea"]==1)

#%%Dynamic geology
geo = geometry.dynamic_confining_layer(geo, bcs["sea"], pars["t_max"])
geo = geometry.create_Kh(geo, **pars)

#%%Data processing for model
#Cut model in half to speed up calculations
geo = geo.sel(y=slice(0, geo.y.max()))
bcs = bcs.sel(y=slice(0, geo.y.max()))

#Cut off unused x and y cells 
#otherwise writing the initial conditions for the next model run is 
#problematic due to the RCB algorithms completely leaving out usused rows and columns
geo = geo.dropna("x", how="all", subset=["tops"]).dropna("y", how="all", subset=["tops"])
bcs = bcs.dropna("x", how="all").dropna("y", how="all")
#%%Create initial conditions
approx_init = True

c_f = 0.0
c_s = 35.

rho_f, rho_s = ic.c2rho(c_f), ic.c2rho(c_s)
shd, sconc = ic.get_ic(bcs, geo, c_f, c_s, approx_init=approx_init)

#%%Some extra processing to make iMOD-python accept these DataArrays
geo = geo.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
bcs = bcs.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
sconc = sconc.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
shd = shd.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)

#%%Time management
start_year = 1999 #Must be minimum 1900 for iMOD-SEAWAT
t_kyear = -1 * (ts * 1000 - ts[0] * 1000)
max_perlen = 8000
sub_ts, sub_ends, sub_splits = time_util.subdivide_time(t_kyear, max_perlen)

#%%Non convergence
crashed_model = 4
cell1 = (25,31,152)

#%%
for mod_nr, (i_start, i_end) in enumerate(zip(sub_splits[:-1], sub_splits[1:])):
    print(".........processing model nr. {}..........".format(mod_nr))
    
    timesteps_mod = [cftime.DatetimeProlepticGregorian(t, 1, 1) for t in (sub_ts[mod_nr]+start_year)]
    endtime=cftime.DatetimeProlepticGregorian(sub_ends[mod_nr]+start_year, 1, 1)
    
    bcs_mod = bcs.isel(time=slice(i_start, i_end)).assign_coords(
            time = timesteps_mod)
    
    geo_mod = geo.isel(time=slice(i_start, i_end)).assign_coords(
            time = timesteps_mod)

    #Select for each timestep 
    time_step_min_conf = geo_mod["lith"].where(geo_mod["lith"] == 2).sum(dim=["x", "y", "layer"]).argmin()
    kh = geo_mod["Kh"].isel(time=time_step_min_conf).drop("time")
    
    #Remove empty layers to reduce amount of .idfs written
    river_stage = bcs_mod["river_stage"].dropna(dim="layer", how="all")
    sea = bcs_mod["sea"].dropna(dim="layer", how="all")

    mname_sub = mname + str(mod_nr)
    
    m = imod.wq.SeawatModel(mname_sub)
    
    if mod_nr == 0:
        starting_head = xr.where(geo["IBOUND"]==1,  shd, -9999.0)
        starting_conc = xr.where(geo["IBOUND"]==1., sconc, -9999.0)
    else:
        year_str = cftime.DatetimeProlepticGregorian(
                sub_ends[mod_nr-1]+start_year, 1, 1).strftime("%Y%m%d%H%M%S")
        starting_head = "bas/head_{}_l?.idf".format(year_str)
        starting_conc = "btn/conc_{}_l?.idf".format(year_str)
    
    #TODO: Refer to model0 for static idfs: IBOUND, ICBUND. Can save 1000 idfs.
    #Does not work for IBOUND?
    m["bas"] = imod.wq.BasicFlow(ibound=geo["IBOUND"].assign_coords(dx=dx, dy=-dy), 
                                 top=topbot[0], 
                                 bottom=xr.DataArray(topbot[1:], {"layer": geo.layer}, ("layer")), 
                                 starting_head=starting_head)
    
    m["lpf"] = imod.wq.LayerPropertyFlow(
        k_horizontal=kh, k_vertical=kh/pars["ani"], specific_storage=0.0
    )
    
    m["btn"] = imod.wq.BasicTransport(
        icbund=geo["IBOUND"], 
        starting_concentration=starting_conc, 
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
    
    m["ghb"] = imod.wq.GeneralHeadBoundary(head = xr.where(sea==1, bcs_mod["sea_level"], np.nan),
                                           conductance=sea * pars["dx"] * pars["dy"] / pars["bc-res"],
                                           density= rho_s, #sea * rho_s,
                                           concentration=sea * c_s)
    
    m["riv"] = imod.wq.River(stage = river_stage,
                             conductance = xr.where(np.isfinite(river_stage), pars["dx"] * pars["dy"] / pars["bc-res"], np.nan),
                             bottom_elevation = river_stage - pars["riv_depth"],
                             density = rho_f, #xr.where(np.isfinite(bcs["river_stage"]), rho_f, np.nan),  
                             concentration =  xr.where(np.isfinite(river_stage), 0., np.nan)) #0.) 
    
    m["pksf"] = imod.wq.ParallelKrylovFlowSolver(
                                                 max_iter=1000, 
                                                 inner_iter=100, 
                                                 hclose=0.0001, 
                                                 rclose=1000., 
                                                 relax=1.00,
                                                 partition="rcb",
                                                 solver="pcg",
                                                 preconditioner="ilu",
                                                 deflate=False,
                                                 debug=False,
                                                 )
    
    m["pkst"] = imod.wq.ParallelKrylovTransportSolver(
                                                 max_iter=1000, 
                                                 inner_iter=30, 
                                                 cclose=1e-6,
                                                 relax=0.98,
                                                 partition="rcb",
                                                 solver="bicgstab",
                                                 preconditioner="ilu",
                                                 debug=False,
                                                 )
    
    m["oc"] = imod.wq.OutputControl(save_head_idf=True, save_concentration_idf=True)
    
    n_timesteps_p1 = 8
    time_util.time_discretization(m, 1000., 
                                  endtime=endtime,
                                  n_timesteps_p1=n_timesteps_p1,
                                  timestep_multiplier=5.)
    
    m.write(directory = os.path.join(r"c:\Users\engelen\test_imodpython\synth_delta_test", mname))

    #%non_conv_analyser
    if mod_nr == crashed_model:
        ncg1, xyz1 = ncg.look_around(m, cell1, n=2)
