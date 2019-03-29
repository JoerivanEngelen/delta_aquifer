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

#%%Helper functions
def _isclose(a, b, atol=1e-05, rtol=1e-8):
    return abs(a - b) <= (atol + rtol * abs(b))


#%%Main functions
def get_sea_level(path_eustatic, ts, figfol=None):
    df = pd.read_csv(path_eustatic, sep=" ", header=9).rename(
        columns={"Age(ka)": "age"}
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
    ds = xr.Dataset.from_dataframe(df.set_index("age"))
    dt = ds.age - ds.age.shift(age=1)  # Calc dts
    dt = dt.fillna(0)
    dt.name = "dt"

    for start, end in zip(ts[:-1], ts[1:]):
        weights_ends = (
            np.array([start, end]) % dt.sel(age=[start, end], method="backfill")
        ).fillna(1.0)

        sl = slice(*(weights_ends.age.values[::-1]))

        vals_ends = ds.interp(age=[start, end])
        weights_ends["age"] = vals_ends["age"]

        if weights_ends[0] > 0.0:
            weights_ends[0] = 1 - weights_ends[0]

        weights = xr.concat([weights_ends, dt.sel(age=sl)], dim="age")

        vals = xr.concat([vals_ends, ds.sel(age=sl)], dim="age")

        means.append((weights * vals).sum() / weights.sum())

    means = xr.concat(means, dim="age")
    means["age"] = (ts[:-1] + ts[1:]) / 2

    return means


def get_coastline(geo, sea_level, phi=None, L = None, a = None, 
                  figfol=None, t_start=None, t_max=None, t_end=None, 
                  tra=None, **kwargs):
    # top_coast
    top_coast = _isclose(geo["tops"], sea_level["50%"], atol=0.5, rtol=2e-1) & geo[
        "edges"
    ].sel(
        z=sea_level["50%"], method="pad"
    )  # Get sea cells in upper layer

    coastline_x = xr.where(top_coast, top_coast.x, np.nan).min(dim="x")

    weights_trans = np.clip((sea_level.age - t_start)/(t_max - t_start), 0, 1)
    weights_reg = np.clip((sea_level.age - t_max)/(t_end - t_max), 0, 1)
    weights = weights_trans - weights_reg
    
    x_trans, y_trans  = geometry._pol2cart(L * a * (1-tra),
                              np.linspace(-phi/2, phi/2, num=coastline_x["y"].shape[0])) 
    x_trans = xr.DataArray(x_trans, coords={"y":y_trans}, dims=["y"]).reindex_like(coastline_x, method="nearest")
    
    coastline_x = weights * x_trans + (1-weights) * coastline_x    
    coastline_x = top_coast.sel(x = coastline_x, method = "nearest").x    
    coastline_x = coastline_x.where(coastline_x<coastline_x.max())#For some reason nans get the maximum x...

    # Find coastline
    coastline = xr.where(top_coast.x == coastline_x, 1, 0)
    # Clip off outer edges, outside cone calculated by polar coordinate system.
    rho_onshore = coastline_x.sel(y=0, method="nearest")
    x_edge, y_edge = geometry._pol2cart(rho_onshore, phi / 2)
    coastline = xr.where((coastline.y < y_edge) & (coastline.y > -y_edge), coastline, 0)

    if figfol is not None:
        coastline.sum(dim="age").plot()
        plt.savefig(os.path.join(figfol, "coastlines.png"))
        plt.tight_layout()
        plt.close()

    return (coastline, coastline_x, rho_onshore)


def get_river_grid(
    geo, sea_level, coastline, rho_onshore, phi=None, L=None, figfol=None, **kwargs
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
        [(float(age), np.meshgrid(h.sel(age=age), phis)[0]) for age in h.age.values]
    )
    h["x"], h["y"] = geometry._pol2cart(*np.meshgrid(rhos, phis))

    h_grid = {}
    h_grid["x"], h_grid["y"] = np.meshgrid(geo.x, geo.y)

    h_grid = geometry.ds_d2_grid(
        geometry.pol2griddata(h, np.zeros((n_inp - 20, n_inp)).astype(bool), h_grid)
    )
    h_grid = xr.concat(
        [h_grid[age].assign_coords(age=age) for age in dhdx.age.values], dim="age"
    )
    h_grid = h_grid.where(h_grid > sea_level)
    return h_grid


#%%Plotting functions
def plot_sea_level_curve(sea_level, df, ts, figfol):
    start, end = int(ts[0]), int(ts[-1])
    plt.plot(sea_level["age"], sea_level["50%"], drawstyle="steps-mid")
    plt.scatter(sea_level["age"], sea_level["50%"])
    plt.plot(df["age"][end:start], df["50%"][end:start])
    ax = plt.gca()
    ax.invert_xaxis()
    plt.savefig(os.path.join(figfol, "sea_level_curve.png"))
    plt.close()
