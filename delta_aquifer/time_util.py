# -*- coding: utf-8 -*-
"""
Created on Wed May  1 13:58:51 2019

@author: engelen
"""

import numpy as np
import imod
import xarray as xr
from datetime import timedelta


def add_timesteps(max_perlen, times, nstp_extra):
    nu_times = []
    for i, nstp in enumerate(nstp_extra):
        start_time = times[i]
        nu_times += [start_time + timedelta(j*max_perlen*365.25) for j in range(0, nstp+1)]
    nu_times.append(times[-1])   
    return(nu_times)

def time_discretization(model, max_perlen, endtime, starttime=None):
    """
    Collect all unique times and subdivde. Adapted from the function in imod/wq/model.py
    """
    
    model.use_cftime = model._use_cftime()

    times = []
    for pkg in model.values():
        if "time" in pkg.coords:
            times.append(pkg["time"].values)

    # TODO: check that endtime is later than all other times.
    times.append(imod.wq.timeutil.to_datetime(endtime, model.use_cftime))
    if starttime is not None:
        times.append(imod.wq.timeutil.to_datetime(starttime, model.use_cftime))

    # np.unique also sorts
    times = np.unique(np.hstack(times))    
    duration = imod.wq.timeutil.timestep_duration(times, use_cftime=True)
    
    # Update times, ensuring that max_perlen is not exceeded.
    nstp_extra = [int(d/(max_perlen*365.25)) for d in duration]
    nu_times = add_timesteps(max_perlen, times, nstp_extra)    
    nu_duration = imod.wq.timeutil.timestep_duration(nu_times, use_cftime=True)
    # Generate time discretization, just rely on default arguments
    # Probably won't be used that much anyway?    
    timestep_duration = xr.DataArray(
        nu_duration, coords={"time": np.array(nu_times)[:-1]}, dims=("time",)
    )
    model["time_discretization"] = imod.wq.TimeDiscretization(
        timestep_duration=timestep_duration
    )

def split_list(ls, splits):
    return([ls[i : j] for i, j in zip([0] + 
              splits, splits + [None])])

def find_cumsum_splits(ls, max_sum):
    split_ls = []
    comb_sum = ls[0]
    for i, p in enumerate(ls[1:]):
        if (comb_sum + p) <= max_sum:
            comb_sum += p
        else:
            split_ls.append((i+1,comb_sum))
            comb_sum=p
    
    return list(zip(*split_ls))

def subdivide_time(t, max_perlen):
    #NOTE THIS FUNCTION DOES NOT WORK IF perlen > max_perlen
    perlens = t[1:] - t[:-1]
    splits, comb_perlens = find_cumsum_splits(perlens, max_perlen)
    
    sub_t = split_list(t, list(splits)+[len(perlens)])
    first_el = np.array([sub[0] for sub in sub_t])
    starts, ends = first_el[:-1], first_el[1:]
    
    sub_ends = ends-starts
    sub_t = [sub-start for sub, start in zip(sub_t, starts)]
    
    #TO DO: Needlessly complex, can better use splits with split.prepend(0) and split.append(len(perlen))
    sub_splits = np.cumsum(np.array([0]+[i.shape[0] for i in sub_t]))
    
    return(sub_t, sub_ends, sub_splits)