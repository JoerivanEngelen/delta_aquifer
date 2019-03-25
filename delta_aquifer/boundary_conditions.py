# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:04:19 2019

@author: engelen
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os

def get_sea_level(path_eustatic, ts, figfol=None):
    df = pd.read_csv(path_eustatic, sep = " ", header = 9).rename(columns={"Age(ka)" : "age"})
    sea_level = calc_weighted_mean(df, ts)

    if figfol is not None:
        start, end = int(ts[0]), int(ts[-1])
        plt.plot(sea_level["age"], sea_level["50%"], drawstyle='steps-mid')
        plt.scatter(sea_level["age"], sea_level["50%"])
        plt.plot(df["age"][end:start], df["50%"][end:start])
        ax = plt.gca()
        ax.invert_xaxis()
        plt.savefig(os.path.join(figfol, "sea_level_curve.png"))
        plt.close()

def calc_weighted_mean(df, ts):
    """Calculate weighted means of a dataframe for N timeslices from a array of timeslice bins (N+1). 
    
    Function accounts for timeslice bin edges that fall in between two datapoints of the dataframe 
    by interpolation and subsequently granting the interpolated datapoint a lower weight.
    """
    means = []
    ds = xr.Dataset.from_dataframe(df.set_index("age"))
    dt = ds.age - ds.age.shift(age=1) #Calc dts
    dt = dt.fillna(0)
    dt.name = "dt"
    
    for start, end in zip(ts[:-1], ts[1:]):
        weights_ends = (np.array([start, end]) % dt.sel(age=[start, end], method="backfill")).fillna(1.)
        
        sl = slice(*(weights_ends.age.values[::-1]))
        
        vals_ends = ds.interp(age=[start, end])
        weights_ends["age"]=vals_ends["age"]
     
        if weights_ends[0] > 0.:
            weights_ends[0] = 1-weights_ends[0]
    
        weights = dt.sel(age=sl)
        weights = xr.concat([weights_ends, weights], dim="age")
        
        vals = ds.sel(age=sl)
        vals = xr.concat([vals_ends, vals], dim="age")
        
        means.append((weights * vals).sum()/weights.sum())

    means = xr.concat(means, dim="age")
    means["age"] = (ts[:-1] + ts[1:])/2

    return(means)
