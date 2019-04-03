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
from delta_aquifer import geometry

#%%Override buggy xarray function with this hack. This hack probably only works within this context.

def _unique_value_groups(ar, sort=True):
    """Group an array by its unique values.

    Parameters
    ----------
    ar : array-like
        Input array. This will be flattened if it is not already 1-D.
    sort : boolean, optional
        Whether or not to sort unique values.

    Returns
    -------
    values : np.ndarray
        Sorted, unique values as returned by `np.unique`.
    indices : list of lists of int
        Each element provides the integer indices in `ar` with values given by
        the corresponding value in `unique_values`.
    """
    inverse, values = pd.factorize(ar, sort=sort)
#    groups = [[] for _ in range(len(values))]
    groups = [[] for _ in range(len(values)+1)] #Add extra list for NaNs
    for n, g in enumerate(inverse):
#        if g >= 0:
        if True:
            # pandas uses -1 to mark NaN, but doesn't include them in values
            # therefore NaNs are stored in the last (extra) list
            groups[g].append(n)
    
    if not groups[-1]:
        groups = groups[:-1]
    
    return values, groups

xr.core.groupby.unique_value_groups = _unique_value_groups

#%%Helper functions
def _isclose(a, b, atol=1e-05, rtol=1e-8):
    return abs(a - b) <= (atol + rtol * abs(b))

def _mid_to_binedges(mids):
    binedges=np.zeros(len(mids)+1)
    delta = mids[1:] - mids[:-1]
    binedges[1:-1] = mids[1:]-delta/2
    binedges[0] = mids[0]-delta[0]/2
    binedges[-1] = mids[-1]+delta[-1]/2
    return(binedges)

def _dumb(x):
    """This function is used as a hack to convert labeled groupby values back to a DataArray.
    """
    return(x)

#%%Main functions
def get_sea_level(path_eustatic, ts, figfol=None):
    df = pd.read_csv(path_eustatic, sep=" ", header=9).rename(
        columns={"Age(ka)": "time"}
    )
    sea_level = calc_weighted_mean(df, ts)

    if figfol is not None:
        plot_sea_level_curve(sea_level, df, ts, figfol)

    return sea_level


def calc_weighted_mean(df, ts):
    """Calculate weighted means of a dataframe for N timeslices from a array of timeslice bins (N+1). 
    
    Function accounts for timeslice bin edges that fall in between two datapoints of the dataframe 
    by interpolation and subsequently granting the interpolated datapoint a lower weight.
    """
    means = []
    ds = xr.Dataset.from_dataframe(df.set_index("time"))
    dt = ds.time - ds.time.shift(time=1)  # Calc dts
    dt = dt.fillna(0)
    dt.name = "dt"

    for start, end in zip(ts[:-1], ts[1:]):
        weights_ends = (
            np.array([start, end]) % dt.sel(time=[start, end], method="backfill")
        ).fillna(1.0)

        sl = slice(*(weights_ends.time.values[::-1]))

        vals_ends = ds.interp(time=[start, end])
        weights_ends["time"] = vals_ends["time"]

        if weights_ends[0] > 0.0:
            weights_ends[0] = 1 - weights_ends[0]

        weights = xr.concat([weights_ends, dt.sel(time=sl)], dim="time")

        vals = xr.concat([vals_ends, ds.sel(time=sl)], dim="time")

        means.append((weights * vals).sum() / weights.sum())

    means = xr.concat(means, dim="time")
    means["time"] = (ts[:-1] + ts[1:]) / 2

    return means


def coastlines(geo, sea_level, phi=None, L = None, a = None, 
                  figfol=None, t_start=None, t_max=None, t_end=None, 
                  tra=None, **kwargs):
    # top_coast
    top_coast = _isclose(geo["tops"], sea_level["50%"], atol=0.5, rtol=2e-1) & geo[
        "edges"
    ].sel(
        z=sea_level["50%"], method="pad"
    )  # Get sea cells in upper layer

    coastline_rho = xr.where(top_coast, top_coast.x, np.nan).min(dim="x").sel(y=0.0, method="nearest")

    weights_trans = np.clip((sea_level.time - t_start)/(t_max - t_start), 0, 1)
    weights_reg = np.clip((sea_level.time - t_max)/(t_end - t_max), 0, 1)
    weights = weights_trans - weights_reg
    
    coastline_rho = L * a * (1-tra) * weights + (1-weights) * coastline_rho
    phis = np.linspace(-phi/2, phi/2, num=top_coast["y"].shape[0])
    phis = xr.DataArray(phis, coords={"phi": phis}, dims=["phi"])
     
    coastline_loc = xr.Dataset(dict([i for i in zip(*[["xc", "yc"], geometry._pol2cart(coastline_rho, phis)])]))

    x_bins, y_bins = [_mid_to_binedges(geo[dim].values) for dim in ["x", "y"]]
    coastline_loc = coastline_loc.groupby_bins("xc", x_bins, labels=geo["x"].values).apply(_dumb)
    coastline_loc = coastline_loc.groupby_bins("yc", y_bins, labels=geo["y"].values).apply(_dumb)
    
    coastline = xr.where((top_coast.x == coastline_loc["xc_bins"]) & (top_coast.y == coastline_loc["yc_bins"]), 1, 0).max(dim="phi")

    #Get rid of all the polar coordinate nonsense and create a more useful coastline_loc
    coastline_loc = xr.Dataset(dict([(r"%s_loc" % dim, xr.where(coastline, coastline[dim], np.nan).max(dim=dim)) for dim in ["x", "y"]]))
    
    if figfol is not None:
        coastline.sum(dim="time").plot()
        plt.savefig(os.path.join(figfol, "coastlines.png"))
        plt.tight_layout()
        plt.close()

    return(coastline, coastline_loc, coastline_rho)

def river_grid(
    geo, sea_level, rho_onshore, phi=None, L=None, figfol=None, **kwargs
):
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    # Create river cells
    apex = geo["tops"].max()

    dhdx = (apex - sea_level) / rho_onshore

    n_inp = 200
    rhos = np.linspace(0, L, num=n_inp)
    rhos = xr.DataArray(rhos, coords=[rhos], dims=["rho"])
    phis = np.linspace(-phi / 2, phi / 2, num=n_inp - 20)
    h = apex - rhos * dhdx

    h = dict(
        [(float(time), np.meshgrid(h.sel(time=time), phis)[0]) for time in h.time.values]
    )
    h["x"], h["y"] = geometry._pol2cart(*np.meshgrid(rhos, phis))

    h_grid = {}
    h_grid["x"], h_grid["y"] = np.meshgrid(geo.x, geo.y)

    h_grid = geometry.ds_d2_grid(
        geometry.pol2griddata(h, np.zeros((n_inp - 20, n_inp)).astype(bool), h_grid)
    )
    h_grid = xr.concat(
        [h_grid[time].assign_coords(time=time) for time in dhdx.time.values], dim="time"
    )
    h_grid = h_grid.where(h_grid > sea_level)
    
    return h_grid

def sea_3d(geo, sea_level, coastline_loc):
    coastline_loc["x_loc"] = coastline_loc["x_loc"].fillna(coastline_loc["x_loc"].min(dim="y"))
    
    return(xr.where(
    (sea_level > geo["edges"].z) & (geo.x <= coastline_loc["x_loc"]), geo["edges"], 0
    ))

def river_3d(
    geo, sea_level, rho_onshore, phi=None, L=None, figfol=None, **kwargs        
):
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    
    h_grid = river_grid(geo, sea_level, rho_onshore, phi, L, figfol, **kwargs)
    
    return((h_grid * geo["edges"].where(geo["edges"] != 0)).dropna(dim="z", how="all"))

#%%Plotting functions
def plot_sea_level_curve(sea_level, df, ts, figfol):
    start, end = int(ts[0]), int(ts[-1])
    plt.plot(sea_level["time"], sea_level["50%"], drawstyle="steps-mid")
    plt.scatter(sea_level["time"], sea_level["50%"])
    plt.plot(df["time"][end:start], df["50%"][end:start])
    ax = plt.gca()
    ax.invert_xaxis()
    plt.savefig(os.path.join(figfol, "sea_level_curve.png"))
    plt.close()
