# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:58:43 2019

@author: engelen
"""

import sys, os
import numpy as np
import pandas as pd
import xarray as xr
import imod

model_fol  = sys.argv[1]
mname = "dummy"

# Basic flow
layer = np.arange(1, 4)
z = np.arange(25.0, 0.0, -10.0)
y = np.arange(4.5, 0.0, -1.0)
x = np.arange(0.5, 5.0, 1.0)
ibound = xr.DataArray(
    np.full((3, 5, 5), 1.0),
    coords={"z": z, "layer": ("z", layer), "y": y, "x": x, "dx": 1.0, "dy": -1.0},
    dims=("z", "y", "x"),
)
starting_head = xr.full_like(ibound, 0.0)
top = 30.0
bot = xr.DataArray(
    np.arange(20.0, -10.0, -10.0), coords={"layer": layer}, dims=("layer",)
)
# BTN
icbund = ibound.copy()

# LPF
k_horizontal = ibound.copy()

# GHB
head = ibound.copy()

# RCH
datetimes = pd.date_range("2000-01-01", "2000-01-05")
rate = xr.DataArray(
    np.full((5, 5, 5), 1.0),
    coords={"time": datetimes, "y": y, "x": x, "dx": 1.0, "dy": -1.0},
    dims=("time", "y", "x"),
)

# DSP
longitudinal = ibound.copy()

# Fill model
m = imod.wq.SeawatModel("test_model")
m["bas6"] = imod.wq.BasicFlow(
    ibound=ibound, top=top, bottom=bot, starting_head=starting_head
)
m["lpf"] = imod.wq.LayerPropertyFlow(
    k_horizontal=k_horizontal,
    k_vertical=k_horizontal.copy(),
    horizontal_anisotropy=k_horizontal.copy(),
    interblock=k_horizontal.copy(),
    layer_type=k_horizontal.copy(),
    specific_storage=k_horizontal.copy(),
    specific_yield=k_horizontal.copy(),
    save_budget=False,
    layer_wet=k_horizontal.copy(),
    interval_wet=0.01,
    method_wet="wetfactor",
    head_dry=1.0e20,
)
m["ghb"] = imod.wq.GeneralHeadBoundary(
    head=head,
    conductance=head.copy(),
    concentration=head.copy(),
    density=head.copy(),
    save_budget=False,
)
m["riv"] = imod.wq.River(
    stage=head,
    conductance=head.copy(),
    bottom_elevation=head.copy(),
    concentration=head.copy(),
    density=head.copy(),
    save_budget=False,
)

m["rch"] = imod.wq.RechargeHighestActive(
    rate=rate, concentration=rate.copy(), save_budget=False
)

m["btn"] = imod.wq.BasicTransport(
    icbund=icbund,
    starting_concentration=icbund.copy(),
    porosity=icbund.copy(),
    n_species=1,
    inactive_concentration=1.0e30,
    minimum_active_thickness=0.01,
)
m["adv"] = imod.wq.AdvectionTVD(courant=1.0)
m["dsp"] = imod.wq.Dispersion(
    longitudinal=longitudinal,
    ratio_horizontal=longitudinal.copy(),
    ratio_vertical=longitudinal.copy(),
    diffusion_coefficient=longitudinal.copy(),
)
m["vdf"] = imod.wq.VariableDensityFlow(
    density_species=1,
    density_min=1000.0,
    density_max=1025.0,
    density_ref=1000.0,
    density_concentration_slope=0.71,
    density_criterion=0.01,
    read_density=False,
    internodal="central",
    coupling=1,
    correct_water_table=False,
)


m["pksf"] = imod.wq.ParallelKrylovFlowSolver(
                                             max_iter=1000, 
                                             inner_iter=100, 
                                             hclose=1e-5, 
                                             rclose=1., 
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

m.write(directory = os.path.join(model_fol, mname))