# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:06:34 2019

@author: engelen
"""

import xarray as xr
import numpy as np

def c2rho(c, denselp=0.7143, rho_ref=1000.):
    return(rho_ref + c * denselp) #0.7142859 more accurate for denselp

def correct_head(sconc, shd, bcs, c_f, c_s):
    "Correcting point water heads by assuming fresh water head is equal with depth."
    
    rho_f= c2rho(c_f)
    srho = c2rho(sconc)
    
    shd = rho_f/srho*shd + (srho-rho_f)/srho * bcs.z
    return(shd)

def ghyb_herz(riv_stage_2d, bcs, shd, c_f, c_s):
    sea_level = bcs["sea_level"].isel(time=0, drop=True)
    sea_loc_xy = bcs["sea"].isel(time=0, drop=True).max(dim="z")
    
    xsea_min = xr.where(sea_loc_xy==1, sea_loc_xy.x, np.nan).min()
    ysea_max = xr.where(sea_loc_xy.sel(x=xsea_min)==1., sea_loc_xy.y, np.nan).max()
    
    rho_f, rho_s = c2rho(c_f), c2rho(c_s)
    z_interface = sea_level + rho_f/(rho_s-rho_f) * (sea_level - riv_stage_2d)
    below_interface = (bcs.z <= z_interface)
    
    sconc = xr.where( below_interface | (~np.isnan(sea_loc_xy) | ((np.abs(sea_loc_xy.y) >= ysea_max) & (sea_loc_xy.x >= xsea_min))), c_s, c_f)
    return(shd, sconc)

def get_ic(bcs, geo, c_f, c_s, approx_init=False):
    """Get initial conditions. Defaults to point water and salinities of zero. 
    If approx_init = True, approximate initial salinities with a ghyben herzberg relationship, sharp interface.
    Consequently correct heads for these initial salinities.
    """
    riv_stage_2d = bcs["riv_stage"].max(dim="z").isel(time=0, drop=True)
    shd = riv_stage_2d.fillna(bcs["sea_level"].isel(time=0, drop=True)) * geo["IBOUND"]
    
    if approx_init == True:
        shd, sconc = ghyb_herz(riv_stage_2d, bcs, shd, c_f, c_s)
        shd = correct_head(sconc, shd, bcs, c_f, c_s)
    else:
        sconc = c_f
    
    return(shd, sconc)