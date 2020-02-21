# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:31:23 2020

@author: engelen
"""
import pandas as pd
import numpy as np

import os
from pkg_resources import resource_filename

def select_re(df, key):
    return(df.loc[:, df.columns.str.contains(key)])
    
def add_logspace(df, columns, new_column, n=1):
    a = np.log10(df.loc[:, columns])
    for i in range(n):
        a.insert(i+1, new_column+"_%s"%i, np.nan)
    a = 10**a.interpolate(axis=1)
    a = select_re(a, new_column)
    return(a)

def average_column(df, key, log=False):
    """Calculate average across column
    and remove columns
    """
    
    min_key = "%s_min" % key
    max_key = "%s_max" % key
    cols = [min_key, max_key]
    
    if log == True:
        df[key] = 10**np.log10(df.loc[:, cols]).mean(axis=1)
    else:
        df[key] = df.loc[:, cols].mean(axis=1)
    df = df.drop(columns = cols)
    return(df)

def discretize_range(df, inp, identifier):
    """Discretize min-max input range and convert this to long format."""
    df = pd.concat([df, add_logspace(df, 
                                     [f"{inp}_min", f"{inp}_max"], 
                                     f"{inp}_mid")], axis=1)
    df = pd.wide_to_long(df.reset_index(), stubnames=f"{inp}", 
                              i=identifier, j=f"{inp}_type", 
                              sep='_', suffix='\w+').reset_index()
    return(df)
    
def calculate_rs(df, Kv, Kh):
    r_s=df["H_b"]*(
            (df["f_aqt"]/(df["N_aqt"]*Kv)
              )+(
                      (1-df["f_aqt"])*df['Anisotropy']/((df["N_aqt"]+1)*df['Kaqf_max']))
            )
    return(r_s)
#%%Path management
datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))
fnames = ["hydrogeology", "lithology", "geometry", "boundary_conditions"]
df_paths  = [os.path.join(datafol, "literature", f+r".csv") for f in fnames]    

model_inputs = os.path.join(datafol, "model_inputs.csv")

#%%Sheets and parameters to use for models
inputs = ["Delta", "L", "l_a", "alpha", "beta", "gamma", "phi_f", "H_b", "f_H",
          "N_aqt", "f_aqt", "l_conf", "N_pal", "s_pal", "l_tra", "t_tra",
          "N_chan", "l_surf_end"]

#%%Read data
sheets = {}

for i, f in enumerate(fnames):
    shtname = os.path.splitext(f)[0]
    sheets[shtname] = pd.read_csv(df_paths[i], skiprows=[1], encoding="latin")

hdrglgy = sheets["hydrogeology"]

#%%Filter cases where reported anisotropies are effective 
#  anisotropies for the complete groundwater system.

ani_eff = hdrglgy["Aqt_explicit"]==0

hdrglgy["Kaqt_min"] = hdrglgy["Kaqt_min"].where(~ani_eff, hdrglgy["Kaqf_min"]/hdrglgy["Anisotropy_max"])
hdrglgy["Kaqt_max"] = hdrglgy["Kaqt_max"].where(~ani_eff, hdrglgy["Kaqf_max"]/hdrglgy["Anisotropy_min"])
hdrglgy["Anisotropy_max"] = hdrglgy["Anisotropy_max"].where(~ani_eff)
hdrglgy["Anisotropy_min"] = hdrglgy["Anisotropy_min"].where(~ani_eff)

#%%Select only columns for calculation (no remarks or references etc.)
hdrglgy = select_re(hdrglgy, "min|max|Delta")

d = {"min|Delta" : min, 
     "max|Delta" : max}

hdrglgy = pd.concat(
        [select_re(hdrglgy, key).groupby(
                by="Delta"
                ).aggregate(func) for key, func in d.items()],
        axis=1)

#Drop Ss, not of influence on salt transport
hdrglgy = hdrglgy.drop(columns=["Ss_min", "Ss_max"])

#%%The model is insensitive to these variables, so we average these
col_dic = {"Anisotropy": True,
           "n": False,
#           "Ss": True,
           "Recharge": True,
           "riv_res": True}

columns = col_dic.keys()

for col in columns:
    hdrglgy = average_column(hdrglgy, col, col_dic[col])

#Guadalfeo aquifer adjacent to Donana
#so we take the inputs from that area
hdrglgy.rename(index={'Guadalfeo':'Doñana'},inplace=True)

#Count NaNs as this is a measure of uncertainty
hy_nans = hdrglgy.isna().sum(axis=1)

#%%NaNs with the median value to limit amount of runs required
hdrglgy.loc[:, columns] = hdrglgy.loc[:, columns].fillna(hdrglgy.loc[:, columns].median())
hdrglgy["Kaqt_min"] = hdrglgy["Kaqt_min"].fillna(hdrglgy["Kaqt_min"].quantile(0.25))
hdrglgy["Kaqt_max"] = hdrglgy["Kaqt_max"].fillna(hdrglgy["Kaqt_max"].quantile(0.75))


#%%Discretize min-max input range and convert this to long format.
#Comment this out to calculate r_s
hdrglgy = discretize_range(hdrglgy, "Kaqf", "Delta")
hdrglgy = discretize_range(hdrglgy, "Kaqt", ["Delta", "Kaqf_type"])

hdrglgy = hdrglgy.rename(columns = {"Kaqf" : "Kh_aqf", "Kaqt" : "Kv_aqt", 
                          "Recharge" : "R", "Anisotropy" : "Kh_Kv"})

#%%Merge and filter
df = sheets[fnames[1]]
for shtname in fnames[2:]:
    df = pd.merge(df, sheets[shtname], on = ["Delta", "Country"])

df = df.loc[:, inputs]
df = df.set_index("Delta").drop("Nakdong") #Too small, too little data

#Count NaNs as this is a measure of uncertainty
df_nans = df.isna().sum(axis=1)

#%%Fill missing data records
#Fill NaN in t_max with delta on same coast with
#roughly similar climate and sediment source
df.loc["Mahanadi", "t_tra"] = df.loc["Krishna", "t_tra"] 

#Create many clay layers and many paleochannels in the GBM delta
#To mimic chaotic lithology
df.loc["Ganges-Brahmaputra", "N_aqt"] = 50
df.loc["Ganges-Brahmaputra", "N_pal"] = 30

#Fill rest of missing data values with median
nanputs = ["f_H", "N_aqt", "f_aqt", "l_conf", "N_pal", "s_pal"]
for inp in nanputs:  
    df[inp] = df[inp].fillna(df[inp].median())

#%%
intputs = ["N_aqt", "N_pal", "N_chan"]
for inp in intputs:
    df[inp] = df[inp].astype(np.int64)

#%%Combine dfs
nans_all = pd.concat([df_nans, hy_nans], axis=1).sum(axis=1)
df_out = pd.merge(df, hdrglgy, on="Delta")

#%%r_s
#rs_min = calculate_rs(df_out, df_out["Kaqt_max"], df_out["Kaqf_min"])
#rs_max = calculate_rs(df_out, df_out["Kaqt_min"], df_out["Kaqf_max"])

#%%Save
df_out.to_csv(model_inputs)