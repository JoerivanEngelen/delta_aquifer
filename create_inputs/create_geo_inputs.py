# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:40:08 2020

@author: engelen
"""

from delta_aquifer import geo_util as gu
import geopandas as gpd
import numpy as np
import pandas as pd

def sel_delta(gdf, delta):
    return(gdf.loc[gdf["Delta"]==delta])

def split_type(gdf, type1, type2):
    return(gdf.loc[gdf["Type"]==type1], gdf.loc[gdf["Type"]==type2])

def get_azi_and_distance(gdf, delta, type1="apex", type2="coast"):
    """Calculate arc distances and azimuth for a delta
    """
    gdf_delta = sel_delta(gdf, delta)
    a, c = split_type(gdf_delta, type1, type2)
    lats_longs = [a.geometry.y.values, a.geometry.x.values,
                  c.geometry.y.values, c.geometry.x.values]
    
    dist = gu.getArcDistance(*lats_longs)
    azi =  gu.getAzimuth(*lats_longs)
    return(dist, azi)

def calculate_phi_f(azi, threshold=180.):
    """From an array of azimuths, calculate angles between lines from apex to coast
    Subsequently correct for coastal points that are both west and east of the North line (e.g. Nile Delta)
    """
    deg2Rad = np.pi/180.
    
    phis = np.subtract.outer(azi, azi-360) % 360
    phis = phis[phis < threshold]
    return(np.max(phis)*deg2Rad)

def azi_to_angle(azi, threshold=180.):
    """Convert an azimuth of e.g. 320 degrees to an angle of -40 degrees 
    while preserving 40 degrees azimuth to 40 degrees.
    """
    return((azi+threshold) % 360 - threshold)

#%%Path management
path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\geometry\delta_points.shp"

#%%Read shapefile
delta_p = gpd.read_file(path)

#%%
deltas = delta_p["Delta"].unique()
df = pd.DataFrame(index = deltas, columns = ["L_a", "L_b", "phi_f"])

for delta in deltas:
    dist, azi = get_azi_and_distance(delta_p, delta)
    dist_shelf, azi_shelf = get_azi_and_distance(delta_p, delta, "coast", "shelf")
    df.loc[delta, "L_a"]   = np.average(dist)
    df.loc[delta, "L_b"]   = np.average(dist_shelf)
    df.loc[delta, "phi_f"] = calculate_phi_f(azi)

#%%Calculate gridsize
df["L_b"] = df["L_b"].clip(upper=200_000.)

df["L"] = df["L_a"]+df["L_b"]

df["dx"] = 1000
df["dx"] = df["dx"].where(df["L"] > 130_000., 500)
df["dx"] = df["dx"].where(df["L"] < 400_000., 2000)

df["dy"] = df["dx"]
