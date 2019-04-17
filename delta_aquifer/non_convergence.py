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

def look_around(model, cell_fortran, n=2, var=["ghb-head", "riv-stage", "khv", "IBOUND"]):
    """Look around cell that gave non-convergence. 
    
    Insert the indexes of the troublesome cell in Fortran (1 based) based indexing!
    (As printed in the model prompt output)
    """

    cell = (cell_fortran[0]-1, cell_fortran[1]-1, cell_fortran[2]-1)
    
    ds = xr.merge([model[v] for v in var])
    
    shape = ds[var[0]].shape
    ds_sel = isel_out_bounds(ds, shape, cell, n)
    return(ds_sel)