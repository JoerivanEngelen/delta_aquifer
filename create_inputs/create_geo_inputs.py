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

def split_type(gdf):
    return(gdf.loc[gdf["Type"]=="apex"], gdf.loc[gdf["Type"]=="coast"])

def get_azi_and_distance(gdf, delta):
    """Calculate arc distances and azimuth for a delta
    """
    gdf_delta = sel_delta(gdf, delta)
    a, c = split_type(gdf_delta)
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
df = pd.DataFrame(index = deltas, columns = ["L", "phi_f"])

for delta in deltas:
    dist, azi = get_azi_and_distance(delta_p, delta)
    df.loc[delta, "L_a"]   = np.average(dist)
    df.loc[delta, "phi_f"] = calculate_phi_f(azi)

#%%TODO
#Determine L_b
#Determine dx, dy