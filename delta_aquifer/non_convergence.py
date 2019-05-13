# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 09:46:32 2019

@author: engelen
"""

import xarray as xr

def get_out_bounds_slice(dim, cell_i, n):
    """return out of bounds slice
    """
    
    lower_bound = max(cell_i-n, 0) #max, because lower_bound may not be lower than 0
    upper_bound = min(cell_i+n+1, dim) 
    
    return(slice(lower_bound, upper_bound))

def isel_out_bounds(ds, shape, cell, n):
    """Look n cells around specific cell, but ignore if cell outside range.
    """
    
    ddims = {}
    for i, dim in enumerate(["layer", "y", "x"]):
        ddims[dim] = get_out_bounds_slice(shape[i], cell[i], n)
    
    return(ds.isel(**ddims))      

def look_around(model, cell_fortran, n=2, 
                var={
         "ghb" : "head", 
         "riv" : "stage", 
         "lpf" : "k_horizontal", 
         "btn" : "icbund"
         }, z = None):
    """Look around cell that gave non-convergence. 
    
    Insert the indexes of the troublesome cell in Fortran (1 based) based indexing!
    (As printed in the model prompt output)
    """

    cell = (cell_fortran[0]-1, cell_fortran[1]-1, cell_fortran[2]-1)
    
    ds = xr.merge([model[pack][v] for pack, v in var.items()])
    
    shape = [ds.dims[dim] for dim in ["layer", "y", "x"]]
    ds_sel = isel_out_bounds(ds, shape, cell, n)
    
    try:
        xyz = [ds_sel.x[n], ds_sel.y[n], ds_sel.layer[n]]
    except IndexError:
        xyz = [ds_sel.x[0], ds_sel.y[0], ds_sel.layer[0]]
    
    if z is not None:
        xyz[2] = z[xyz[2]-1]
    
    return(ds_sel, xyz)