# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:04:19 2019

@author: engelen
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt
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
    
    #Have to delete empty list as subsequent combine will not work otherwise.
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
def get_sea_level(sl_curve, ts, qt="50%", figfol=None):
    df = pd.read_csv(sl_curve, sep=" ", header=9).rename(
        columns={"Age(ka)": "time"}
    ).set_index("time")
    
    df = df - df.iloc[0, :] #Correct sea levels so that the last timestep is 0.
    
    sea_level = calc_weighted_mean(df, ts, qt)

    if figfol is not None:
        plot_sea_level_curve(sea_level, df, ts, qt, figfol)

    return sea_level


def calc_weighted_mean(df, ts, qt):
    """Calculate weighted means of a dataframe for N timeslices from a array of timeslice bins (N+1). 
    
    Function accounts for timeslice bin edges that fall in between two datapoints of the dataframe 
    by interpolation and subsequently granting the interpolated datapoint a lower weight.
    """
    means = []
    ds = xr.Dataset.from_dataframe(df)[qt]
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
    dx = geo.x[-1]-geo.x[-2]
    # top_coast
    top_coast = (
            _isclose(geo["tops"], sea_level, atol=0.5, rtol=2e-1) | (geo.x >= (L-dx))
            ) & geo["edges"].sel(
        z=sea_level, method="pad"
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
    """Create grid with river stages. 
    """
    
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    
    # Create river cells
    apex = geo["tops"].max()

    dhdx = (apex - sea_level) / rho_onshore

    coords = {}
    coords["rho"], coords["phi"] = geometry._cart2pol(geo["x"], geo["y"])
    
    h_grid = apex - coords["rho"] * dhdx
    h_grid = h_grid.where((h_grid > sea_level) & (geo["IBOUND"].max(dim="z") == 1))
    
    return h_grid

def sea_3d(geo, sea_level, coastline_loc):
    coastline_loc["x_loc"] = coastline_loc["x_loc"].fillna(coastline_loc["x_loc"].min(dim="y"))
    
    sea_z = geo.sel(z= sea_level, method="pad").z

    sea_edge = xr.where(sea_z > geo["edges"].z, geo["edges"], 0)
    sea_trans = ((geo.x >= coastline_loc["x_loc"]) & (geo.z == sea_z)) * geo["IBOUND"]
    
    return((sea_edge|sea_trans).astype(np.int16), sea_z)

def perturb_sea_conc(sea, conc_sea, seed=100, noise_frac=0.01, conc_fresh = 0.):
    """Mathematical formula taken from:
        Simmons et al. [1999]
        On a test case for density-dependent 
        groundwater flow and solute transport models: 
            The salt lake problem
    """
    np.random.seed(seed)
    flat_sea = sea.max(dim=["z", "time"])
    noise = np.random.rand(*flat_sea.shape)
    return(flat_sea*conc_sea + noise_frac * (conc_sea-conc_fresh) * (noise - 0.5))
    
def river_3d(
    geo, sea_level, rho_onshore, phi=None, L=None, figfol=None, **kwargs        
):
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    
    h_grid = xr.Dataset({"h_grid" : river_grid(geo, sea_level, rho_onshore, phi, L, figfol, **kwargs)})
    z_bins = _mid_to_binedges(geo["z"].values)

    h_grid = h_grid.groupby_bins("h_grid", z_bins, labels=geo.layer).apply(_dumb).rename({"h_grid_bins" : "h_l"})

    top_active_layer = xr.where(geo["IBOUND"]==1, geo.layer, np.nan).min(dim="z")

    #Ensure river layer does not exceed IBOUND.
    h_grid["h_l"] = xr.where(h_grid["h_l"] < top_active_layer, top_active_layer, h_grid["h_l"]) 

    #Displace h_grid_bins 1 down, except the maximum which should not go lower than the sea
    #This to avoid h_l from going outside the model domain    
    ismax_hl = (h_grid["h_l"] == h_grid["h_l"].max(dim=["x", "y"]))
    h_grid["h_l"] = xr.where(ismax_hl, h_grid["h_l"], h_grid["h_l"]+1)
    
    riv = h_grid * geo["IBOUND"].where((geo["IBOUND"] == 1) & (geo.layer == h_grid["h_l"]))
    
    return(riv, z_bins)

#%%Master function
def boundary_conditions(sl_curve, ts, geo, conc_sea, conc_fresh, conc_noise = 0.01, 
                        qt = "50%", figfol=None, ncfol=None, **kwargs):
    # Get sea level
    sea_level = get_sea_level(sl_curve, ts, qt=qt, figfol=figfol)
    
    #Start from first timestep where sea level exceeds aquifer bottom again.
    ts_inactive = np.argwhere((sea_level<geo.z.isel(z=0)).values)
    if len(ts_inactive)>0:
        clip_nr=ts_inactive[-1][0]+1
        print("Clipped off everything before timestep nr {} as sea_level lower than aquifer bottom".format(clip_nr))
        sea_level = sea_level.isel(time=slice(clip_nr, None))
    
    # Find active sea cells where GHB's should be assigned.
    coastline, coastline_loc, rho_onshore = coastlines(
        geo, sea_level, figfol=figfol, **kwargs
    )
    sea_cells, sea_z = sea_3d(geo, sea_level, coastline_loc)
    
    #Create river stages
    rivers, z_bins   = river_3d(geo, sea_level, rho_onshore, figfol=figfol, **kwargs)
    
    #Combine to dataset
    bcs = xr.Dataset({"sea": sea_cells, "river_stage" : rivers["h_grid"], "sea_level" : sea_level})
    bcs["sea"] = xr.where(np.isnan(bcs["river_stage"]), sea_cells, 0) #Make sure there are no river cells overlapping sea cells
    bcs["sea_conc"] = perturb_sea_conc(bcs["sea"], conc_sea, noise_frac=conc_noise, conc_fresh=conc_fresh)
    
    bcs = bcs.transpose("time", "z", "y", "x")
    
    if ncfol is not None:
        #Paraview only accepts monotonically increasing times, we have a decreasing one, so has to be changed.
        time = bcs["time"] 
        bcs["time"] = np.arange(len(bcs["time"]))
        #It also requires some CFTIME unit. This one is bullshit.
        bcs["time"].attrs["units"] = "hours since 2000-01-01 00:00:00.0"
        bcs.to_netcdf(os.path.join(ncfol, "bcs.nc"), unlimited_dims = ["time"])
        #Switch back to time in ka again
        bcs["time"] = time
    
    return(bcs)

#%%Plotting functions
def plot_sea_level_curve(sea_level, df, ts, qt, figfol):
    start, end = int(ts[0]), int(ts[-1])
    plt.plot(sea_level["time"], sea_level, drawstyle="steps-mid")
    plt.scatter(sea_level["time"], sea_level)
    plt.plot(df.index[end:start].values, df[qt][end:start].values)   
    ax = plt.gca()
    ax.invert_xaxis()
    plt.savefig(os.path.join(figfol, "sea_level_curve.png"))
    plt.close()
