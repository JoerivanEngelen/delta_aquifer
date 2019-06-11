# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:15:21 2019

@author: engelen
"""
import numpy as np
import pandas as pd
from delta_aquifer import geometry, time_util
from delta_aquifer import boundary_conditions as bc
from delta_aquifer import non_convergence as ncg
from delta_aquifer import initial_conditions as ic
import os, sys

import imod
import xarray as xr
import cftime

from pkg_resources import resource_filename

#%%TODO Tests
#-Wat gebeurt er als zand boven de shelf legt?

#%%TODO Processen
#-Brijn onderin?
#
#-Recharge
#-Conductance waar grote rivier. 1 cel dik.
#->Recharge + conductance combinatie? 
#-Drain? Om te helpen (zoute) rivieren te laten infiltreren
#
#-Zoute rivieren (Savenije methode) -> Overal in waaier
#-Aantal riviertakken.
#
#-Conductiviteit preferente stroombaan door mariene klei (confining layer ook)
#-Breedte preferente stroombaan

#%%Path management
#model_fol = r"c:\Users\engelen\test_imodpython\synth_delta_test"
#sim_nr = 86

model_fol = sys.argv[1]
sim_nr = int(sys.argv[2])

mname = "SD_i{:03d}".format(sim_nr)

figfol = os.path.join(model_fol, mname, "input", "figures")
ncfol  = os.path.join(model_fol, mname, "input", "data")

os.makedirs(figfol, exist_ok=True)
os.makedirs(ncfol,  exist_ok=True)

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))

spratt = os.path.join(datafol, "spratt2016.txt")

#%%Parameters
fixpars = pd.read_csv(os.path.join(datafol, "fixed_pars.csv"), index_col=0).iloc[[0]]
varpars = pd.read_csv(os.path.join(datafol, "traj_real.csv" ), index_col=0).iloc[[sim_nr]]
pars = pd.concat([fixpars.reset_index(), varpars.reset_index()], axis=1)
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
hclose = 1e-4
rclose = pars["dx"] * pars["dy"] * hclose * 10.

#%%Get geometry
geo = geometry.get_geometry(figfol=figfol, ncfol=None, **pars)

topbot=bc._mid_to_binedges(geo["z"].values)[::-1]

#%%Create boundary conditions
# Path management
c_f = pars["c_f"]
c_s = pars["c_s"]

bcs = bc.boundary_conditions(spratt, ts, geo, c_s, c_f, 
                             conc_noise = 0.05, figfol=figfol, ncfol=None, **pars)
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
geo["active"] = geo["IBOUND"].where(geo["IBOUND"]==1.)
geo = geo.dropna("x", how="all", subset=["active"]).dropna("y", how="all", subset=["active"])
bcs = bcs.dropna("x", how="all", subset=["sea", "river_stage"]).dropna("y", how="all", subset=["sea", "river_stage"])
#%%Create initial conditions
approx_init = True

rho_f, rho_s = ic.c2rho(c_f), ic.c2rho(c_s)
shd, sconc = ic.get_ic(bcs, geo, c_f, c_s, approx_init=approx_init)

#%%Time management
start_year = 1999 #Must be minimum 1900 for iMOD-SEAWAT
t_kyear = -1 * (ts * 1000 - ts[0] * 1000)
max_perlen = 8000
sub_ts, sub_ends, sub_splits = time_util.subdivide_time(t_kyear[-(len(bcs.time)+1):], max_perlen)

#%%Save ncs
time_util.num2date_ds(t_kyear[-len(bcs.time):], bcs, geo)

bcs["lith"] = geo["lith"].astype(np.float64)
bcs.to_netcdf(os.path.join(ncfol, "bcs.nc"))
bcs = bcs.drop(["lith"])

#%%Some extra processing to make iMOD-python accept these DataArrays
geo = geo.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
bcs = bcs.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
sconc = sconc.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)
shd = shd.swap_dims({"z" : "layer"}).drop("z").sortby("layer").sortby("y", ascending=False)

#%%
bcs["heads"] = xr.where(bcs["sea"]==1, bcs["sea_level"], bcs["river_stage"])
bcs["conc"] = xr.where(np.isfinite(bcs["river_stage"]), 0., bcs["sea"] * bcs["sea_conc"])
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
#    time_step_min_conf = geo_mod["lith"].where(geo_mod["lith"] == 2).sum(dim=["x", "y", "layer"]).argmin()
#    kh = geo_mod["Kh"].isel(time=time_step_min_conf).drop("time")
    kh = geo["Kh"].isel(time=0).drop("time")

#    mname_sub = mname + str(mod_nr)
    mname_sub = "{}_nr{:02d}".format(mname, mod_nr)
    
    m = imod.wq.SeawatModel(mname_sub)
    
    if mod_nr == 0:
        starting_head = xr.where(geo["IBOUND"]==1,  shd,   -9999.0)
        starting_conc = xr.where(geo["IBOUND"]==1., sconc, -9999.0)
    else:
        year_str = cftime.DatetimeProlepticGregorian(
                sub_ends[mod_nr-1]+start_year, 1, 1).strftime("%Y%m%d%H%M%S")
        starting_head = "bas/head_{}_l?.idf".format(year_str)
        starting_conc = "btn/conc_{}_l?.idf".format(year_str)
    
    #TODO: Refer to model0 for static idfs: IBOUND, ICBUND. Can save 1000 idfs.
    #Does not work for IBOUND?
    m["bas"] = imod.wq.BasicFlow(ibound=geo["IBOUND"].assign_coords(dx=pars["dx"], dy=-pars["dy"]), 
                                 top=topbot[0], 
                                 bottom=xr.DataArray(topbot[1:], {"layer": geo.layer}, ("layer")), 
                                 starting_head=starting_head)
    
    m["lpf"] = imod.wq.LayerPropertyFlow(
        k_horizontal=kh, k_vertical=kh/pars["ani"], specific_storage=0.0,
        save_budget=True,
    )
    
    m["btn"] = imod.wq.BasicTransport(
        icbund=geo["IBOUND"], 
        starting_concentration=starting_conc, 
        porosity=pars["por"]
    )
    
    m["adv"] = imod.wq.AdvectionTVD(courant=0.9)
    m["dsp"] = imod.wq.Dispersion(
            longitudinal=pars["al"], 
            ratio_horizontal=pars["trpt"],
            ratio_vertical=pars["trpv"],
            diffusion_coefficient=pars["diff"]
    )
    
    m["vdf"] = imod.wq.VariableDensityFlow(density_concentration_slope=0.7143)
    
    m["ghb"] = imod.wq.GeneralHeadBoundary(head = bcs_mod["heads"],
                                           conductance=pars["dx"] * pars["dy"] / pars["bc-res"],
                                           density=ic.c2rho(bcs_mod["conc"]), 
                                           concentration=bcs_mod["conc"])
        
    m["pksf"] = imod.wq.ParallelKrylovFlowSolver(
                                                 max_iter=1000, 
                                                 inner_iter=100, 
                                                 hclose=hclose, 
                                                 rclose=rclose, 
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
    
    m["oc"] = imod.wq.OutputControl(save_head_idf=True, 
                                    save_concentration_idf=True, 
                                    save_budget_idf=True)
    
    n_timesteps_p1 = 10
    time_util.time_discretization(m, 1000., 
                                  endtime=endtime,
                                  n_timesteps_p1=n_timesteps_p1,
                                  timestep_multiplier=7.)
    
    m.write(directory = os.path.join(model_fol, mname))

#    #%non_conv_analyser
#    if mod_nr == crashed_model:
#        ncg1, xyz1 = ncg.look_around(m, cell1, n=2)
