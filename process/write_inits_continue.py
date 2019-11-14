# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:36:12 2019

@author: engelen

Script to generate input for a continuation run based on previous output

"""

import xarray as xr
from glob import glob
from imod import idf
import sys, os
import numpy as np
import cftime

#%%Path management

if len(sys.argv) > 1:
    res_folder  = sys.argv[1]
    model_folder = sys.argv[2] 
else:
    #Local testing on my own windows laptop
    res_folder = r"g:\synthdelta\results\test_output\synth_SD_i029_m24_7065039"
    model_folder = r"c:\Users\engelen\test_imodpython\synth_delta_test\SD_i029"
    
files = glob(os.path.join(res_folder, "results_[0-9][0-9][0-9].nc"))
files.sort()

nr_mod = glob(os.path.join(res_folder, "nr_mod_[0-9][0-9].nr"))[0]
with open(nr_mod) as f:
    nr_mod=int(f.readline())

cont_nr = len(files)

mname = os.path.basename(model_folder)
sub_model_path = os.path.join(model_folder, mname+"_nr{:02d}".format(cont_nr))

tempfile=os.path.join(model_folder, "cont_nr.txt")

#%%Get appropriate time

#Created in create_model.py. Update this if changed in create_model.py
sub_ends=np.array([8000., 8000., 5000., 5000., 8000., 8000., 4000.])[-(nr_mod+1):]
n_ts = (sub_ends/1000.).astype(np.int64)
start_year = 1999

time = cftime.DatetimeProlepticGregorian(
        sub_ends[cont_nr-1]+start_year, 1, 1)

#%%Process data
#Load
ds_last = xr.open_dataset(files[cont_nr-1])
#Check if this file has all the timesteps, if not go one back and rename wrong file.
if len(ds_last.time) != n_ts[cont_nr-1]:
    ds_last.close()
    os.rename(files[cont_nr-1], files[cont_nr-1]+"_fail")
    cont_nr -= 1
    ds_last = xr.open_dataset(files[cont_nr-1])

ds_ini = ds_last.isel(time=-1)[["conc1","conc2", "head"]].load()

#Expand conc with species
ds_ini = ds_ini.rename(name_dict={"conc1" : "conc"})
ds_ini["conc"] = ds_ini["conc"].expand_dims(species=[1,2])
ds_ini["conc"] = xr.where(ds_ini["conc"].species==2, 
      ds_ini["conc2"], ds_ini["conc"])
ds_ini = ds_ini.drop(labels="conc2")

#Ensure appropriate time is assigned (the call to update_timesteps.py messed this up)
ds_ini = ds_ini.assign_coords(time=time)

#Write .IDFS
idf.save(os.path.join(sub_model_path, "bas", "head"), ds_ini["head"])
idf.save(os.path.join(sub_model_path, "btn", "conc"), ds_ini["conc"], 
         pattern = r"{name}_{time:%Y%m%d%H%M%S}_c{species}_l{layer}{extension}")

#%%Save cont_nr in temporary file
with open(tempfile, mode="w") as f:
    f.write(str(cont_nr))