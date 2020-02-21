# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:56:56 2020

@author: engelen
"""

from delta_aquifer import geometry as gm
import xarray as xr
import rasterio as rio
import geopandas as gpd
from create_inputs import create_geo_inputs as cgi
from shapely.geometry import Polygon
import pandas as pd
from imod.prepare import spatial as spat
from imod.prepare import Regridder
import imod
import numpy as np

import os
from pkg_resources import resource_filename

def assign_i_coordinate(da, dim, i_dim):
    """Assign dependent index coordinate to keep track of row or column nrs in
    selections
    
    da : DataArray
    
    dim : string
        dimension
    
    i_dim : string
        name of new index coordinate
    
    """
    
    i_coord = np.arange(len(da[dim]))
    i_coord = xr.DataArray(i_coord, coords=da[dim].coords)
    coords = {i_dim : i_coord}
    return(da.assign_coords(**coords))

def calc_ddim(da, dim):
    shifts = {dim : 1}
    ddim = da[dim] - da[dim].shift(**shifts)
    ddim[0] = ddim[1]
    return(ddim)

def get_gcps_coords(L_a, phi_f):
    r = np.array([0., L_a, L_a, L_a])
    phi = np.array([0., -phi_f/2, 0, phi_f/2])
    return(gm._pol2cart(r, phi))

def get_targgrid(dcell, L_a, phi_f, da):
    """Modified version of the function in gm.get_targgrid
    In this function we compensate for when the regridded DataArrays are larger 
    than targgrid calculated by gm.get_targgrid
    
    """
    targgrid = {}
    targgrid["x"], targgrid["y"] = gm.get_targgrid(dcell, dcell, 
            L_a, phi_f/2)
    
    x_max = np.max([da["xt"].max().values, np.max(targgrid["x"])])
    y_max = np.max([da["yt"].max().values, np.max(targgrid["y"])])
    y_min = np.min([da["yt"].min().values, np.min(targgrid["y"])])   
    
    x_max = gm._round2cell(x_max, dcell) + 0.5*dcell
    x_out = np.arange(0.5 * dcell, x_max+1, dcell)
    y_max = gm._round2cell(y_max, dcell) + 0.5*dcell
    y_min = gm._round2cell(y_min, dcell) - 0.5*dcell
    y_out = np.arange(y_min, y_max+1, dcell)
    return(np.meshgrid(x_out, y_out))
    
    

#%%Path management
#Path to annual groundwater abstractions data here. File too large to incorporate in git repo now
path_ann_abs = r"g:\Global_Data\PCR-GLOB_output\2019\totalGroundwaterAbstraction_annuaTot_output.nc"
path_2D_nc = os.path.join(path_ann_abs, "..", "totalGroundwaterAbstraction_{}.nc")

datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
path_shp = os.path.join(datafol, "geometry", "delta_points.shp")
path_gm  = os.path.join(datafol, "geometry.csv")

path_out = os.path.join(datafol, "..","..", "..", "Data", "30_deltas_abs")

#%%Options
save_2D_nc = False #Save seperate rasters

#%%Constants
deg2Rad = np.pi/180.

#%%Read
print("...reading...")
ds = xr.open_dataset(path_ann_abs)
points = gpd.read_file(path_shp).drop(columns = ["id"])
points = points[points["Type"] != "shelf"]
geom = pd.read_csv(path_gm, index_col=0)

#%%For easy checking NetCdfs in qgis:
if save_2D_nc:
    years = ds.time.values.astype('datetime64[Y]').astype(int) + 1970
    
    for i, year in enumerate(years):
        ds.isel(time=i).to_netcdf(path_2D_nc.format(year))

#%%Convert to polygons
deltas = points["Delta"].unique()

df = pd.DataFrame(index = deltas, columns = ["L_a", "phi_f", "azi_mid", "dist_mid", "geometry"])

for delta in deltas:
    dist, azi = cgi.get_azi_and_distance(points, delta)
    df.loc[delta, "L_a"]     = np.average(dist)
    df.loc[delta, "phi_f"]   = cgi.calculate_phi_f(azi)
    df.loc[delta, "azi_mid"] = np.sort(cgi.azi_to_angle(azi))[1]
    df.loc[delta, "dist_mid"] = np.sort(dist)[1]
    df.loc[delta, "geometry"] = Polygon(cgi.sel_delta(points, delta)["geometry"].values)

gdf = gpd.GeoDataFrame(df, geometry="geometry")
gdf["geometry"] = gdf["geometry"].convex_hull #Fix polygons if constructed with points in the wrong order

#%%Prepare for iMOD-Python prepare
#rasterize wants the x dimension to be named "x".
rnm = {"longitude" : "x", "latitude" : "y"}
ds = ds.rename(rnm)

#Assign cellwidths so that iMOD-python can interpret the non-equidistant grid
coords = {"dy" :  calc_ddim(ds, "y"),
          "dx" :  calc_ddim(ds, "x")}
ds = ds.assign_coords(coords=coords)

#%%Select all deltas with bound boxes to reduce data
bnds = gdf["geometry"].bounds #long live this attribute
#Enlarge the bounding box with one cell on each side.
bnds["minx"] -= np.max(ds.dx).values
bnds["maxx"] += np.max(ds.dx).values
bnds["miny"] += np.min(ds.dy).values #dlatitude is negative
bnds["maxy"] -= np.min(ds.dy).values

var = "total_groundwater_abstraction"
ds[var] = ds[var].fillna(0.0)

das = dict([(delta, ds[var].sel(x = slice(b["minx"], b["maxx"]),
                                y = slice(b["maxy"], b["miny"]))
              ) for delta, b in bnds.iterrows()])


#%%Regrid
print("...regridding...")
re_das = {}
    
for delta, da in das.items():
    x = np.arange(bnds.loc[delta]["minx"], bnds.loc[delta]["maxx"], step=0.0083333336)
    y = np.arange(bnds.loc[delta]["maxy"], bnds.loc[delta]["miny"], step=-0.0083333336)
    like = xr.DataArray(None, coords={"x" : x, "y" : y}, 
                        dims=["x", "y"])
    re_das[delta] = Regridder(method="multilinear").regrid(da, like) #Flux is in m/m2/year, so can use multilinear interpolation

#%%Clip shapes (Probably better not to do in a later stage, but useful for testing to see what's going on.)
#ras_shapes = dict([(
#        delta, spat.rasterize(gdf.loc[[delta]], da.isel(time=-1))
#                    ) for delta, da in re_das.items()])
#
#print("...clipping...")
#re_das = dict([(delta, da*ras_shapes[delta]) for delta, da in re_das.items()])
##Cutoff extra borders on the side.
#re_das = dict([(
#        delta, da.dropna(dim="x", how="all").dropna(dim="y", how="all")
#            ) for delta, da in re_das.items()])

#%%Prepare points for transformation, later to be used as gcps
points_trans = points.loc[:, ["Delta", "Type", "mid_coast"]]
points_trans["x"] = points.geometry.x
points_trans["y"] = points.geometry.y   
apexes = points.loc[points["Type"]=="apex"].set_index("Delta")

for delta in deltas:
    apex = apexes.loc[delta].geometry
    points_trans.loc[points_trans["Delta"] == delta, "x"] -= apex.x
    points_trans.loc[points_trans["Delta"] == delta, "y"] -= apex.y

points_trans["r"], points_trans["phi"] = gm._cart2pol(points_trans["y"], points_trans["x"])

#%%Transform
print("...transforming...")
mid_coast = points.loc[points["mid_coast"]==1].set_index("Delta")

for delta, da in re_das.items():
    apex = apexes.loc[delta].geometry
    mid = mid_coast.loc[delta].geometry
    da = da.assign_coords(x = da.x-apex.x, y=da.y-apex.y) #Check if apex is really at y=0...
    r, phi   = gm._cart2pol(*np.meshgrid(da.y, da.x))
    r_mid, _ = gm._cart2pol(mid.y-apex.y, mid.x-apex.x)
    
    #Scale r
    r *= df.loc[delta]["dist_mid"]/r_mid
    points_trans.loc[points_trans["Delta"] == delta, "r"] *= df.loc[delta]["dist_mid"]/r_mid
    
    #Rotate delta, so it points northwards
    phi = (phi/deg2Rad-df.loc[delta]["azi_mid"])*deg2Rad
    phi_p = points_trans.loc[points_trans["Delta"] == delta, "phi"]
    points_trans.loc[points_trans["Delta"] == delta, "phi"] = (phi_p/deg2Rad-df.loc[delta]["azi_mid"])*deg2Rad
    
    r, phi = [xr.DataArray(
            data = i, coords = {"x" : da.x, "y": da.y}, dims=["x", "y"]
            ) for i in [r, phi]]
    
    xt, yt = gm._pol2cart(r, phi)
    
    da = da.assign_coords(r = r, phi = phi, xt=xt, yt=yt)
    re_das[delta] = da

points_trans["xt"], points_trans["yt"] = gm._pol2cart(points_trans["r"], points_trans["phi"])

#%%Interpolate to grid created for delta
print("...interpolating...")
model_data = {}

for delta, da in re_das.items():
    dcell = geom.loc[delta, "dx"]
    assert(dcell == geom.loc[delta, "dy"])
    
    targgrid = {}
    targgrid["x"], targgrid["y"] = get_targgrid(dcell,
            df.loc[delta, "L_a"], df.loc[delta, "phi_f"],
            da)
    
    poldata = {}
    poldata["x"], poldata["y"] = da.xt.values.T, da.yt.values.T
    poldata["abstraction"] = da.isel(time=-1).values
    nan_idx = np.isnan(poldata["abstraction"])
    
    griddata = gm.pol2griddata(poldata,nan_idx, targgrid)
    coords = {"x" : griddata["x"][0], "y" : griddata["y"][:, 0]}
    
    #TODO: Interpolate all timesteps instead of just the last one 
    model_data[delta] = xr.DataArray(data=griddata["abstraction"].T, 
              coords=coords, dims=["x", "y"])

    #Assign coordinates which we use for the control points
    model_data[delta] = assign_i_coordinate(model_data[delta], "x", "ix")
    model_data[delta] = assign_i_coordinate(model_data[delta], "y", "iy")

#%%Create ground control points:
print("...creating gcps...")
gcps = {}

#Fix order of coastal points
points_trans = points_trans.sort_values(["Delta", "Type", "yt"])

for delta in deltas:
    x_sel = points_trans.loc[points_trans["Delta"] == delta, "xt"]
    y_sel = points_trans.loc[points_trans["Delta"] == delta, "yt"]
    x_sel, y_sel = xr.DataArray(x_sel.values, dims="select"), xr.DataArray(y_sel.values, dims="select")
    select = model_data[delta].sel(x=x_sel, y=y_sel, method="nearest") #ix and iy are rows and cols for gcp
    xg, yg = get_gcps_coords(df.loc[delta, "L_a"], df.loc[delta, "phi_f"])
    gcps[delta] = [rio.control.GroundControlPoint(row=r, col=c, x=x, y=-y) for r,c,x,y in zip(select.iy.values, select.ix.values, xg, yg)]

#%%Warp
print("...warping...")
warp_data = {}

for delta, da in model_data.items():
    data = model_data[delta].drop(labels=["ix", "iy"]).transpose("y", "x") #imod-python does not want other dimensions
    data = data.assign_coords(y=data.y*-1) #Flip y-coorindates
    dst = xr.full_like(data,np.nan)
    warp_data[delta] = xr.full_like(data,np.nan)
    src_crs = dst_crs = rio.crs.CRS.from_epsg(32630)
    dst_transform = imod.util.transform(dst)
    warp_data[delta].values, _ = rio.warp.reproject(data.values, dst.values, src_crs=src_crs, 
             dst_crs=dst_crs, dst_transform=dst_transform, gcps=gcps[delta])

#%%Clip out delta
print("...clipping...")
for delta, da in warp_data.items():
    dcell = geom.loc[delta, "dx"]
    
    targgrid["x"], targgrid["y"] = gm.get_targgrid(dcell, dcell,
                df.loc[delta, "L_a"], df.loc[delta, "phi_f"]/2)
    coords = {"x" : targgrid["x"][0], "y" : targgrid["y"][:, 0]}
    like = xr.DataArray(np.ones(targgrid["x"].shape), coords=coords, dims=["y", "x"])
    da = da.sortby(da.y) #Ensure y is monotonically increasing
    warp_data[delta] = da * like
    warp_data[delta] = warp_data[delta].fillna(0.0)

    

#%%Save
print("...saving...")
for delta, da in warp_data.items():
    path_nc = os.path.join(path_out, "{}.nc".format(delta))
    da.to_netcdf(path_nc)
