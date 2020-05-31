# -*- coding: utf-8 -*-
"""
Created on Sun May 31 14:16:36 2020

@author: engelen

This script uses the Anthropogenic biomes dataset, which can be found here:
    http://sedac.ciesin.columbia.edu/es/anthropogenicbiomes.html
And created by:
    Ellis, E. C. and N. Ramankutty. 2008. Putting people in the map: anthropogenic biomes of the world. 
    Data distributed by the Socioeconomic Data and Applications Center (SEDAC)
"""

import xarray as xr
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, box, mapping
from pyproj import CRS

import numpy as np
from create_inputs import create_geo_inputs as cgi

import os
from pkg_resources import resource_filename

#%%Functions
def clip_da(delta, gdf, da):
    try:
        da = da.rio.clip(
                gdf.loc[[delta]].geometry.apply(mapping), gdf.crs, drop=True, invert=False
                )
    except rxr.exceptions.OneDimensionalRaster:
        pass
    return(da)

#%%Path management
datafol = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
path_shp            = os.path.join(datafol, "geometry", "delta_points.shp")
path_classification = os.path.join(datafol, "classification", "delta_classification.shp")

#TODO: Check if data present and if not, download from server. In this way do not have to distribute data.
biomes_path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\40_deltas\GlobalAnthro\anthromes_v1-2001-2006-geotif\gl_anthrome.tif"

crs = CRS.from_epsg(4326)

#%%Read
points = gpd.read_file(path_shp).drop(columns = ["id"])
points = points[points["Type"] != "shelf"]

biomes = rxr.open_rasterio(biomes_path).sel(band=1).reset_coords(drop=True) #Only one band
biomes = biomes.rio.set_crs(crs)
#%%Convert to polygons
#I just copied this from the groundwater_abstractions.py 
deltas = points["Delta"].unique()

df = pd.DataFrame(index = deltas, columns = ["geometry"])

for delta in deltas:
    df.loc[delta, "geometry"] = Polygon(cgi.sel_delta(points, delta)["geometry"].values)

gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
gdf["geometry"] = gdf["geometry"].convex_hull #Fix polygons if constructed with points in the wrong order

#%%Get cellsizes
dx, dy = (0.08333333333333329, 0.08333333333333329)

#%%Group Anthrome into:
#1: Dense settlements
#2: Villages (Rice villages & Irrigate villages)
#3: Croplands
#4: Rangelands
#5: Seminatural
#6: Wildlands
#255: NoData

biomes = biomes.where(biomes==255, (biomes/10).astype(np.int8))

#%%Select all deltas with bound boxes to reduce data
bnds = gdf["geometry"].bounds #long live this attribute
#Enlarge the bounding box with one cell on each side.
bnds["minx"] -= dx
bnds["maxx"] += dx
bnds["miny"] -= dy
bnds["maxy"] += dy

das = dict([(delta, biomes.sel(x = slice(b["minx"], b["maxx"]),
                               y = slice(b["maxy"], b["miny"])
                               ).rio.set_crs(crs)
              ) for delta, b in bnds.iterrows()])

#%%Clip shapes 
das = dict([(delta, clip_da(delta, gdf, das[delta])
    ) for delta in das.keys()])

#%%Count values and concatenate into one DataFrame
cell_counts = pd.concat([das[delta].to_series().value_counts().to_frame(delta) for delta in das.keys()], axis=1)

cell_counts = cell_counts.drop(255)

#%%Further group into:
#Urban: (Dense settlements)
#Agriculture: (Villages + Croplands)
#Nature: (Rangelands + Seminatural + Wildlands)

group_dict = {1: "Urban", 
              2: "Agriculture", 
              3: "Agriculture", 
              4: "Nature", 5:"Nature", 6:"Nature"}

cell_counts = (cell_counts.reset_index()
.replace({'index': group_dict})
.groupby('index', sort=False).sum()
)

#%%Calculate fractional 
cell_frac = cell_counts / cell_counts.sum(axis=0)
cell_frac = cell_frac.T

#%%Classify deltas for map
# Hub: More than 20% landsurface covered with "Urban"
# Breadbasket: More than 50% landsurface covered with Agriculture
# Nature: More than 33% landsurface covered with Nature
deltas_classification = pd.DataFrame()
deltas_classification["Hub"] = cell_frac["Urban"] > 0.20
deltas_classification["Breadbasket"] = cell_frac["Agriculture"] > 0.5
deltas_classification["Habitat"] = cell_frac["Nature"] > 0.33

deltas_classification = deltas_classification.astype(np.int8)

#Due to coarse data, we have to correct some small deltas:
deltas_classification.loc["Donana"] = [0,0,1]

#%%Create shapefile
deltas_classification["geometry"] = gdf.geometry.centroid

shp_class = gpd.GeoDataFrame(deltas_classification, geometry="geometry", crs=crs)
shp_class = shp_class.reset_index().rename(columns={"index": "Delta"})

shp_class.crs = {"init" : "epsg:4326"}

shp_class.to_file(path_classification)