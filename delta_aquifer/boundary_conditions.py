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

#%%Override xarray function with this hack that does not work with NaNs 
#This hack probably only works within this context.

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

def _get_inv_d1(d1):
    
    #z=0.0 occurs twice, we have to remove this
    mean_rho_z0=d1["top"].where(d1["top"]==0.0, drop=True).rho.mean()
    z, index = np.unique(d1["top"], return_index=True)
    
    
    d1_inv = xr.Dataset({"rho" : (["z"], d1["rho"].values)}, coords={"z": d1["top"].values})

#%%Get sea level
def get_sea_level(sl_curve, ts, qt="50%", figfol=None, ext=".png"):
    df = pd.read_csv(sl_curve, sep=" ", header=9).rename(
        columns={"Age(ka)": "time"}
    ).set_index("time")
    
    df = df - df.iloc[0, :] #Correct sea levels so that the last timestep is 0.
    
    sea_level = calc_weighted_mean(df, ts, qt)

    if figfol is not None:
        plot_sea_level_curve(sea_level, df, ts, qt, figfol, ext=ext)

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

#%%Boundary condition location
def coastlines(geo, d1, sea_level, phi=None, L = None, L_a = None, 
                  figfol=None, t_start=None, t_max=None, t_end=None, 
                  tra=None, **kwargs):
    dx = geo.x[-1]-geo.x[-2]
    
    #Flip 
    d1_inv = xr.Dataset({"rho" : (["z"], d1["rho"].values)}, coords={"z": d1["top"].values})
    #two times z = 0.0 causes error.
    d1_inv.sel(z=sea_level, method="pad")
#    # top_coast
#    top_coast = (
#            #TODO: test if making this dependent on dz is more robust
#            _isclose(
#            geo["tops"], sea_level, atol=0.5, rtol=5e-1
#                    ) | (
#            geo.x >= (L-dx)
#                    )
#            # Get sea cells in upper layer
#            ) & geo["edges"].sel(z=sea_level, method="pad")
#    
#    coast_y = xr.where(top_coast, top_coast.x, np.nan).min(dim="x")
#            
#    coastline_rho = coast_y.sel(y=0.0, method="nearest")
#    #Sometimes at y=0.0 there are no cells located, so fill these with the minimum across y.
#    coastline_rho = xr.where(np.isnan(coastline_rho), 
#                             coast_y.min(dim="y"), coastline_rho)
#
#    if np.any(np.isnan(coastline_rho)):
#        raise ValueError("coastline_rho is nan at certain moments in time")

    weights_trans = np.clip((sea_level.time - t_start)/(t_max - t_start), 0, 1)
    weights_reg = np.clip((sea_level.time - t_max)/(t_end - t_max), 0, 1)
    weights = weights_trans - weights_reg
    
    coastline_rho = L_a * (1-tra) * weights + (1-weights) * coastline_rho
    phis = np.linspace(-phi/2, phi/2, num=top_coast["y"].shape[0])
    phis = xr.DataArray(phis, coords={"phi": phis}, dims=["phi"])
     
    coastline_loc = xr.Dataset(dict([i for i in zip(*[["xc", "yc"], geometry._pol2cart(coastline_rho, phis)])]))

    x_bins, y_bins = [_mid_to_binedges(geo[dim].values) for dim in ["x", "y"]]
    coastline_loc = coastline_loc.groupby_bins("xc", x_bins, labels=geo["x"].values).apply(_dumb)
    coastline_loc = coastline_loc.groupby_bins("yc", y_bins, labels=geo["y"].values).apply(_dumb)
    
    coastline = xr.where(
            (
            top_coast.x == coastline_loc["xc_bins"]
            ) & (
            top_coast.y == coastline_loc["yc_bins"]
            ), 1, 0).max(dim="phi")

    #Get rid of all the polar coordinate nonsense and create a more useful coastline_loc
    coastline_loc = xr.Dataset(dict([(r"%s_loc" % dim, xr.where(coastline, coastline[dim], np.nan).max(dim=dim)) for dim in ["x", "y"]]))
    
    if figfol is not None:
        coastline.sum(dim="time").plot()
        plt.savefig(os.path.join(figfol, "coastlines.png"))
        plt.tight_layout()
        plt.close()

    return(coastline, coastline_loc, coastline_rho)

def river_grid(
    geo, sea_level, coastline_rho, phi=None, L=None, figfol=None, **kwargs
):
    """Create grid with river stages. 
    """
    
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    
    # Create river cells
    apex = geo["tops"].max()

    dhdx = (apex - sea_level) / coastline_rho

    coords = {}
    coords["rho"], coords["phi"] = geometry._cart2pol(geo["x"], geo["y"])
    
    h_grid = apex - coords["rho"] * dhdx
    h_grid = h_grid.where(
            (h_grid > sea_level) & (geo["IBOUND"].max(dim="z") == 1)
            )
    
    return h_grid, dhdx

def sea_3d(geo, sea_level, coastline_loc):
    #This can probably be simplified just by using a coastline_rho instead of these coastlines...
    coastline_loc["x_loc"] = coastline_loc["x_loc"].fillna(coastline_loc["x_loc"].min(dim="y"))
    
    sea_z = geo.sel(z= sea_level, method="pad").z

    sea_edge = xr.where(sea_z > geo["edges"].z, geo["edges"], 0)
    sea_trans = ((geo.x >= coastline_loc["x_loc"]) & (geo.z == sea_z)) * geo["IBOUND"]
    
    return((sea_edge|sea_trans).astype(np.int16), sea_z)
    
def river_3d(
    geo, sea_level, coastline_rho, phi=None, L=None, figfol=None, **kwargs        
):
    assert type(sea_level) == xr.core.dataarray.DataArray
    assert type(geo) == xr.core.dataarray.Dataset
    
    h_grid, dhdx = river_grid(geo, sea_level, coastline_rho, phi, L, figfol, **kwargs)
    h_grid = xr.Dataset({"h_grid" : h_grid})
    z_bins = _mid_to_binedges(geo["z"].values)

    h_grid = h_grid.groupby_bins("h_grid", z_bins, labels=geo.layer).apply(_dumb).rename({"h_grid_bins" : "h_l"})

    top_active_layer = xr.where(geo["IBOUND"]==1, geo.layer, np.nan).min(dim="z")

    #Ensure river layer does not exceed IBOUND.
    h_grid["h_l"] = xr.where(h_grid["h_l"] < top_active_layer, top_active_layer, h_grid["h_l"]) 

    riv = h_grid * geo["IBOUND"].where((geo["IBOUND"] == 1) & (geo.layer == h_grid["h_l"]))
    
    return(riv, z_bins, dhdx)

def create_channel_mask(d2_ds, N_chan, phi=None, L = None, **kwargs):
    n_inp=200
    rhos = np.linspace(0, L, num=n_inp)
    channel = geometry._cake_cuts(N_chan, phi, rhos, f_offset = 0.)
    channel = [c.flatten() for c in channel]

    da = d2_ds.sel(x=xr.DataArray(channel[0], dims="foo"),
                   y=xr.DataArray(channel[1], dims="foo"), 
                   method="nearest")
    
    channel_mask = ((d2_ds.x == da.x) & (d2_ds.y == da.y)).max(dim="foo")
    return(channel_mask)

#%%Hydrogeology
def calc_riv_conductance(coastline_loc, f_cond_chan=None, N_chan=None, 
                     dx=None, dy=None, bc_res=None, **kwargs):
    #Create channels
    base_cond = dx*dy/bc_res
    chan_cond = base_cond * f_cond_chan
    channel_mask = create_channel_mask(coastline_loc, N_chan, **kwargs)
    conductance = xr.where(channel_mask, chan_cond, base_cond)
    return(conductance, base_cond)

#%%Salinity surface water   
def perturb_sea_conc(sea, c_s, seed=100, noise_frac=0.01, c_f = 0.):
    flat_conc=sea.max(dim=["z", "time"])*c_s
    return(perturb_conc(flat_conc, c_s, seed, noise_frac, c_f))

def perturb_riv_conc(riv_conc, c_s, seed=100, noise_frac=0.01, c_f=0.):
    riv_perturbed = perturb_conc(riv_conc,  
                        c_s, seed=100, noise_frac=0.01, c_f = 0.)
    
    return(xr.where(riv_conc > 0., riv_perturbed, riv_conc))

def perturb_conc(flat_conc, c_s, seed=100, noise_frac=0.01, c_f = 0.):
    """Mathematical formula taken from:
        Simmons et al. [1999]
        On a test case for density-dependent 
        groundwater flow and solute transport models: 
            The salt lake problem
    """
    np.random.seed(seed)
    noise = np.random.rand(*flat_conc.shape)
    return(flat_conc + noise_frac * (c_s-c_f) * (noise - 0.5))
    
def correct_salinity_intrusion(intrusion_length, dhdx):
    """Correct salinity intrusion in surface water for the changing gradient
    
    This equation is derived from substituting a factor in the Chezy formula.
    With this formula we can calculate Qf, on which the the intrusion length L
    depends. See formulas 5.47 and 5.48 in Huub Savenije's Salinity and Tides
    in Alluvial Estuaries. Freely available here:
    https://hubertsavenije.files.wordpress.com/2014/01/salinityandtides2_21.pdf
    
    """
    f = dhdx/dhdx.isel(time=-1)
    
    f_L = np.log(1/np.sqrt(f)+1)/np.log(2)
    #Normalize over np.log(2) as we want f_L at time=-1  to be 1.

    return(f_L * intrusion_length)

def estuary_profile_slope(intrusion_rho, coastline_rho, c_s, c_f):
    """Create linear salinity profile for the estuary,
    
    dumbing down the Savenije equation from 3 parameters to 1 parameter
    
    """
    return((c_s - c_f)/(coastline_rho - intrusion_rho))
    
def salinity_profile(rhos_2d, intrusion_rho, coastline_rho, c_s, c_f):
    """
    """

    est_slope = estuary_profile_slope(intrusion_rho, 
                                 coastline_rho,
                                 c_s, c_f)

    concs_rho =  est_slope * (rhos_2d-intrusion_rho)

    estuary_salinity = xr.where((rhos_2d > (intrusion_rho)) & (rhos_2d < coastline_rho), 
                         concs_rho, c_f)
    estuary_salinity = xr.where((rhos_2d > coastline_rho), c_s, estuary_salinity)
    return(estuary_salinity.drop(labels=["z", "layer"]))

#%%Recharge
def recharge(onshore_mask, rch_rate):
    return(onshore_mask * rch_rate)

#%%Master function
def boundary_conditions(sl_curve, ts, geo, c_s=None, c_f=None, 
                        bc_res=None, N_chan=None, f_cond_chan=None,
                        L_a=None, intrusion_L=None, rch_rate=None, 
                        conc_noise=0.01, qt="50%", 
                        figfol=None, ncfol=None, **kwargs):
    # Get sea level
    sea_level = get_sea_level(sl_curve, ts, qt=qt, figfol=figfol)
    min_sea_level = np.min(sea_level) #Can be used for initial conditions
    
    #Start from first timestep where sea level exceeds aquifer bottom again.
    ts_inactive = np.argwhere((sea_level<geo.z.isel(z=0)).values)
    if len(ts_inactive)>0:
        clip_nr=ts_inactive[-1][0]+1
        print("Clipped off everything before timestep nr {} as sea_level lower than aquifer bottom".format(clip_nr))
        sea_level = sea_level.isel(time=slice(clip_nr, None))
    
    # Find active sea cells where GHB's should be assigned.
    coastline, coastline_loc, coastline_rho = coastlines(
        geo, sea_level, figfol=figfol, L_a=L_a, **kwargs
    )
    
    sea_cells, sea_z = sea_3d(geo, sea_level, coastline_loc)
    
    #Create river stages
    rivers, z_bins, dhdx = river_3d(geo, sea_level, coastline_rho, figfol=figfol, **kwargs)
    
    #Salinity intrusion in surface water
    intrusion_rho = coastline_rho - correct_salinity_intrusion(intrusion_L * L_a, dhdx)
    coords = {} #Should add these somewhere in geometry.py as dependent coordinates
    coords["rho"], coords["phi"] = geometry._cart2pol(geo["x"], geo["y"])
    estuary_salinity = salinity_profile(coords["rho"], intrusion_rho, coastline_rho, c_s, c_f) 
    
    #Calculate river conductance
    riv_conductance, base_cond = calc_riv_conductance(coastline_loc, 
                    f_cond_chan=f_cond_chan, N_chan=N_chan, 
                    bc_res=bc_res, **kwargs)
    
    #Combine to dataset
    bcs = xr.Dataset({"sea" : sea_cells, 
                      "riv_stage" : rivers["h_grid"], 
                      "sea_level" : sea_level})
    
    #Ensure that there are no river cells overlapping or overlying sea cells
    riv_mask = np.isfinite(bcs["riv_stage"])
    bcs["river"] = xr.where(riv_mask, 1, np.nan).astype(np.int16)
    bcs["sea"] = xr.where((~riv_mask.sum(dim="z"))&(bcs["sea"]==1), sea_cells, np.nan) 

    #Peturb concentrations (above 0.)
    sea_salinity = perturb_sea_conc(bcs["sea"], c_s, 
        noise_frac=conc_noise, c_f=c_f)
    bcs["sea_conc"]  = xr.where(bcs["sea"]==1, sea_salinity, np.nan)
    bcs["sea_cond"]  = xr.where(bcs["sea"]==1, base_cond, np.nan)
    
    estuary_salinity = perturb_riv_conc(estuary_salinity, c_s,
        noise_frac=conc_noise, c_f=c_f).clip(min=0.0)
    bcs["riv_conc"] = xr.where(riv_mask, estuary_salinity, np.nan)
    bcs["riv_cond"] = xr.where(riv_mask, riv_conductance, np.nan)
    
    #Recharge
    onshore_mask = (
            geo["IBOUND"].max(dim="z")-sea_cells.max(dim="z").fillna(0.)
            )
    bcs["rch"] = recharge(onshore_mask, rch_rate)
    
    #Put dimensions in right order
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
    
    return(bcs, min_sea_level)

#%%Plotting functions
def plot_sea_level_curve(sea_level, df, ts, qt, figfol, ext=".png"):
    start, end = int(ts[0]), int(ts[-1])
    plt.plot(sea_level["time"], sea_level, drawstyle="steps-mid")
    plt.scatter(sea_level["time"], sea_level)
    plt.plot(df.index[end:start].values, df[qt][end:start].values)   
    ax = plt.gca()
    ax.invert_xaxis()
    plt.savefig(os.path.join(figfol, "sea_level_curve%s"%ext))
    plt.close()
