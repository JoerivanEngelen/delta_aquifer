# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:47:59 2020

@author: engelen
"""

import xarray as xr
import pandas as pd
import numpy as np
import os

from pkg_resources import resource_filename

from imod.prepare import spatial

import geopandas as gpd
from shapely.geometry import Polygon

def get_porosity(par, delta):
    return(par.loc[par["Delta"]==delta]["n"].mean()) #9 times the same value, so mean() just used to reduce to scalar.

def read_poly(delta_points, delta, crs):
    
    points = delta_points.loc[
        (delta_points["Type"] != "shelf") & (delta_points["Delta"] == delta)
            ]
    points = points.to_crs(crs)
    
    poly = gpd.GeoDataFrame(
            pd.DataFrame(dict(delta=delta, 
                              geometry=Polygon(points["geometry"].values)
                              ), index=[0]), geometry="geometry"
            )
            
    #Fix polygons if constructed with points in the wrong order
    poly["geometry"] = poly["geometry"].convex_hull
    return(poly)

def _rev_dims(da, *dims):
    """Reverse dims, 
    alternative to all the slow sortby actions that slowed this package down previously
    """
    
    kwargs = dict([(dim, da[dim][::-1]) for dim in dims])
    return(da.reindex(**kwargs))

#%%Path management
paths_Nile = [
        r"g:\3D_Nile_Delta_keep\ND_Paper_Anthropocene_Closed_Brine_Bot_Paleo_v3_n2m48_5335333",
        r"g:\3D_Nile_Delta_keep\ND_Paper_Anthropocene_Closed_noCL_Brine_Top_Paleo_v3_n2m48_5335365",
        r"g:\3D_Nile_Delta_keep\ND_Paper_Anthropocene_Half_Open_Brine_Top_Paleo_v3_n2m48_5334560",
        r"g:\3D_Nile_Delta_keep\ND_Paper_Anthropocene_Half_Open_condCL_Brine_Top_Paleo_v3_n2m48_5335601",
        r"g:\3D_Nile_Delta_keep\ND_Paper_Anthropocene_Half_Open_noCL_Brine_Top_Paleo_v3_n2m48_5335745"
        ]

paths_Nile = [os.path.join(path, "results_001.nc") for path in paths_Nile]
names_Nile = ["C-M-B-P", "C-N-B-P", "H-M-T-P", "H-F-T-P", "H-N-T-P"]

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
outfol = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "validation")))
geom_shp = os.path.join(datafol, "geometry", "delta_points.shp")

path_cl = r"g:\NL_Data\3dchloride_for_crosssections.nc"

par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Open and process chloride concentrations
cl = xr.open_dataset(path_cl)
cl = cl["3d-chloride"]
cl = cl.swap_dims({"layer" : "z"})
cl = _rev_dims(cl, "z")

cl = cl.sel(z=slice(-300, 20)) #select 

#Convert chloride to TDS
tds_Rhine = cl * (35/16500) #convert mg CL/L to g TDS/L

#%%Read Nile volumes 
#.to_dataset(name=name)
tds_Nile = [xr.open_dataset(path, decode_times=False)["conc1"].isel(time=-1).assign_coords(model=name) for path, name in zip(paths_Nile, names_Nile)]
tds_Nile = xr.concat(tds_Nile, dim="model")
tds_Nile = tds_Nile.sel(z=slice(10, -300))
tds_Nile = _rev_dims(tds_Nile, "y")

#%%Get Rhine shape
delta_points = gpd.read_file(geom_shp)
Rhine_poly = read_poly(delta_points, "Rhine-Meuse", "EPSG:28992")
Nile_poly = read_poly(delta_points, "Nile", "EPSG:22992")

#%%Rasterize Rhine-Meuse delta
Rhine_raster = spatial.rasterize(Rhine_poly, cl.isel(z=0, percentile=1))
Nile_raster = spatial.rasterize(Nile_poly, tds_Nile.isel(z=0, model=0))

#%%Mask everything outside synthetic model
tds_Rhine = Rhine_raster * tds_Rhine
tds_Nile  = Nile_raster  * tds_Nile
#%%Calculate volumes
FW_vols_tot = {}

por_Rhine = get_porosity(par, "Rhine-Meuse")
vols_Rhine = cl.dx * cl.dy * -1 * cl.dz * por_Rhine
FW_vols_Rhine = xr.where(tds_Rhine<5, vols_Rhine, np.nan) 
FW_vols_tot["Rhine-Meuse"] = FW_vols_Rhine.sum(dim=["x", "y", "z"])

por_Nile = get_porosity(par, "Nile")
vols_Nile = 1000 * 1000 * 20 * por_Nile
FW_vols_Nile = xr.where(tds_Nile<5, vols_Nile, np.nan)
FW_vols_tot["Nile"] = FW_vols_Nile.sum(dim=["x", "y", "z"])

#%%Save
for name, ds in FW_vols_tot.items():
    ds.to_dataframe(name="FW_vols").to_csv(os.path.join(outfol, "FW_vol_{}.csv".format(name)))