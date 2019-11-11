# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:06:34 2019

@author: engelen
"""

import xarray as xr
import numpy as np

def c2dens(c, denselp=0.7143, dens_ref=1000.):
    return(dens_ref + c * denselp) #0.7142859 more accurate for denselp

def correct_head(sconc, shd, bcs, C_f, C_s):
    "Correcting point water heads by assuming fresh water head is equal with depth."
    
    dens_f= c2dens(C_f)
    sdens = c2dens(sconc)
    
    shd = dens_f/sdens*shd + (sdens-dens_f)/sdens * bcs.z
    return(shd)

def ghyb_herz(riv_stage_2d, bcs, shd, C_f, C_s):
    sea_level = bcs["sea_level"].isel(time=0, drop=True)
    sea_loc_xy = bcs["sea"].isel(time=0, drop=True).max(dim="z")
    
    xsea_min = xr.where(sea_loc_xy==1, sea_loc_xy.x, np.nan).min()
    ysea_max = xr.where(sea_loc_xy.sel(x=xsea_min)==1., sea_loc_xy.y, np.nan).max()
    
    dens_f, dens_s = c2dens(C_f), c2dens(C_s)
    z_interface = sea_level + dens_f/(dens_s-dens_f) * (sea_level - riv_stage_2d)
    below_interface = (bcs.z <= z_interface)
    
    sconc = xr.where( 
            below_interface | (~np.isnan(sea_loc_xy) | (
            (np.abs(sea_loc_xy.y) >= ysea_max) & (sea_loc_xy.x >= xsea_min)
            )
        ), C_s, C_f)
    return(shd, sconc)

def get_ic(bcs, geo, C_f=None, C_s=None, approx_init=False, deep_salt=None, 
           **kwargs):
    """Get initial conditions. Defaults to point water and salinities of zero. 
    
    Parameters
    ----------
    approx_init : bool 
        If True, approximate initial salinities with a ghyben herzberg 
        relationship, sharp interface. Consequently correct heads for 
        these initial salinities.
    
    deep_salt : float
        If not None, initially everything below this depth is 
    """
    riv_stage_2d = bcs["riv_stage"].max(dim="z").isel(time=0, drop=True)
    shd = riv_stage_2d.fillna(bcs["sea_level"].isel(time=0, drop=True)) * geo["IBOUND"]
    
    if approx_init == True:
        shd, sconc = ghyb_herz(riv_stage_2d, bcs, shd, C_f, C_s)
        if deep_salt is not None:
            assert(type(deep_salt) in [float, xr.core.dataarray.DataArray, np.ndarray])
            sconc = xr.where(sconc.z<deep_salt, C_s, sconc)
        shd = correct_head(sconc, shd, bcs, C_f, C_s)
    else:
        sconc = xr.DataArray(C_f)
    
    return(shd, sconc)