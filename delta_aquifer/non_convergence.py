# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 09:46:32 2019

@author: engelen
"""

import xarray as xr

dim_order=["layer", "y", "x"]

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
    for i, dim in enumerate(dim_order):
        ddims[dim] = get_out_bounds_slice(shape[i], cell[i], n)
    
    return(ds.isel(**ddims))      

def look_around(model, cell_fortran, n=2, 
                var={
         "ghb" : "conductance",
         "lpf" : "k_horizontal", 
         "btn" : "icbund"
         }, z = None):
    """Look around cell that gave non-convergence. 
    
    Parameters
    ----------
    model : imod.wq.model.SeawatModel
        The model object with all its boundary conditions etc.
    
    cell_fortran : tuple
        Tuple with non-converging cell. Should be exactly as reported in the prompt
        output, so with the Fortran 1-based indexing.
        
    n : int
        Number of cells to look around the cell. 
        n=0 returns just the cell, n=1 a 3x3x3 cube, n=2 a 5x5x5 cube etc.
        
        If the cell is located at the model boundary a smaller cube will be returned.
        
    var : dict
        dictionary with variables to look at as values and the package they belong to as key. 
        
        Downside of this implementation is that you can only have one variable per package.
    """

    cell = tuple([i-1 for i in cell_fortran])
            
    ds = xr.merge([model[pack][v] for pack, v in var.items()])

    if z is not None:
        ds = ds.assign_coords(z = ("layer", z))
    
    shape = [ds.dims[dim] for dim in dim_order]
    ds_sel = isel_out_bounds(ds, shape, cell, n)
    
    xyz = ds.isel(**dict([i for i in zip(dim_order, cell)])).coords
    
    return(ds_sel, xyz)