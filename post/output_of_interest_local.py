# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 10:08:53 2019

@author: engelen

This script is to rerun calculation of the outputs of interest, starting from fracs.nc.
Saves a lot of computational time if fracs.nc are calculated correctly.
"""

import xarray as xr
import pandas as pd
import re
import os
from glob import glob

def coord_of_max(da):
    return(da.where(da==da.max(), drop=True)[-1].time)

#%%Path management

fracs_fol = r"g:\synthdelta\results\trends"
oi_fol = r"g:\synthdelta\results\outputs_of_interest"
t_start = 12

paths = glob(os.path.join(fracs_fol, "fracs_*.nc"))

#%%Loop
for fracs_path in paths:
    #Path management
    match = re.search(r"(i[0-9]{3}).+([0-9]{7})", os.path.basename(fracs_path))
    oi_name = match.group(0)
    oi_path = os.path.join(oi_fol, "oi_%s.csv" % oi_name)
    
    #Read
    fracs = xr.open_dataset(fracs_path)
    
    #%%Differentiate
    grad_fw = fracs["fw"].differentiate("time") * -1 #Multiply with -1 because the time axis is decreasing
    
    #%%Outputs of interest
    oi = {}
    oi["end_fw"]          = fracs["fw"].isel(time=-1).values
    oi["offshore_fw"]     = fracs["fw_offshore"].isel(time=-1).values
    oi["max_fw_decrease"] = fracs["fw"].max().values - oi["end_fw"] 
    oi["old_water"]       = fracs["m_ol_sal"].isel(time=-1).values
    oi["onshore_sw"]      = fracs["m_sal_onshore"].isel(time=-1).values/(35./1025.)
    oi["fw_gradient"]     = (grad_fw).isel( 
                                time=slice(-3, None)
                                ).mean().values
    oi["delay"] = (t_start - coord_of_max(fracs["fw"]*-1).time ).values
    
    keys, values = list(zip(*oi.items()))
    oi = pd.DataFrame(data={"var" : keys, "value" : values}).set_index("var")
    
    oi.to_csv(oi_path)
    
    #%%Close files that are open
    fracs.close()