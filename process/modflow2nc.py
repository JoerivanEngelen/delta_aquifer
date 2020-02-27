# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:40:54 2020

@author: engelen

Process MODFLOW files to .nc and initital .IDFs for iMOD-SEAWAT
Only useful for a private branch of IMOD-WQ where .IDFs are read and MODFLOW
files are written, to make up for the >100k of .IDFs written as modeloutput.

"""

import flopy
from itertools import islice
import re
import pandas as pd
import io
from glob import glob
import numpy as np
import xarray as xr
import cftime
from collections import defaultdict
from imod import idf
import os, sys

def idrange2coord(start, stop, dcell):
    """Convert range of ids to coordinates
    """
    return(np.arange(start,stop)*dcell+dcell/2.)

def find_line(pattern, header):
    return([line for line in header if re.match(pattern, line) is not None])

def pattern2df(pattern, path_listfile):
    with open(path_listfile) as listfile:
        header = islice(listfile, 0, 20000)
        lines = find_line(pattern, header)
    return(fwf2df(lines))

def fwf2df(str_list):
    fwf = io.StringIO("".join(str_list))
    return(pd.read_fwf(fwf, header = None))

def combine_all(ds_list, flip_y=True):
    ds_tot = ds_list[0]
    for ds in ds_list[1:]:
        ds_tot = ds_tot.combine_first(ds)
    
    if flip_y:
        ds_tot.assign_coords(
                y= ds_tot.y.sel(y=slice(None, None, -1))
                )
    
    return(ds_tot)

def parse_species(paths):
    return([int(re.search("MT3D00(\d).UCN", path).group(1)) for path in paths])

def add_metadata(da, z, dz, dcell):
    """Assign some extra metadata for iMOD-python
    """
    da = da.assign_coords(z=xr.ones_like(da.layer) * z)
    da = da.assign_coords(dz=xr.ones_like(da.layer) * dz)
    da["dx"] = dcell
    da["dy"] = -dcell
    return(da)

def calc_fresh_water_head(head, conc, dense_ref=1000., denselp=0.7143):
    """See Post, Kooi & Simmons (2007)
    Using hydraulic head measurements in variable-density ground water flow analyses
    """
    rho_i = dense_ref+conc*denselp
    return(rho_i/dense_ref * head - (rho_i - dense_ref)/dense_ref * head.z)

#%%TODO
#-Create workflow (path management)

#%%Path management
if len(sys.argv) > 1:
    modelfol  = sys.argv[1]
    mod_nr = int(sys.argv[2])  
else:
    ##For Testing
    modelfol = r"g:\test_UCN\RD_i153_nr00"
    mod_nr = 0

run_path=os.path.join(modelfol, "*{}.run".format(mod_nr))

r = re.compile("([a-zA-Z]+_i[0-9]+_nr)([0-9]+)")
mname =  r.match(os.path.basename(modelfol)).group(1)

globpath_ucn = os.path.join(modelfol, r"MT3D00[0-9].UCN.p*")
globpath_hds = os.path.join(modelfol, r"*.hds.p*")
path_l = os.path.join(modelfol, "seawat_tmp", r"{}{:02d}.list.p000".format(mname.lower(), mod_nr))
path_dummy_idf = os.path.join(modelfol, "bas", "ibound_l1.idf")

paths_hds = glob(globpath_hds)
paths_hds.sort()

paths_ucn = glob(globpath_ucn)
paths_ucn.sort()

species = parse_species(paths_ucn)

d1 = defaultdict(list)
for k, v in zip(species, paths_ucn):
    d1[k].append(v)

#%%Cellsize
dummy = idf.open(path_dummy_idf)
dcell = dummy["dx"].values

#%%Read UCN data
ucn_data = {}
for s, paths in d1.items():
    ucn_data[s] = dict([(
            i, flopy.utils.binaryfile.UcnFile(path).get_alldata()
            ) for i, path in enumerate(paths)])

#%%Read hds data
hds_data = [flopy.utils.binaryfile.HeadFile(path).get_alldata() for path in paths_hds]

#%%Parse listfile to extract subdomains
df_pks = pattern2df(r" p\d{03} :", path_l)
df_mod = df_pks.loc[:, :8] #Overlap 1 cell for implicit solvers
df_mt  = df_pks.loc[:, 10:] #Overlap 2 cells for advection package

#Convert FORTRAN 1-based indexing to Python 0-based indexing
df_mod.loc[:, [5,6]] -= 1
df_mod.loc[:, [7,8]] -= 1

df_mt.loc[:, [13, 14]] -= 1
df_mt.loc[:, [15, 16]] -= 1

#%%Parse listfile for other coordinates...
df_top = pattern2df(r" TOP ELEVATION OF LAYER 1 =", path_l)
df_btm = pattern2df(r"   MODEL LAYER BOTTOM EL. =", path_l)

topbots = pd.concat([df_top[6], df_btm[5]]).reset_index(drop=True)

z = (topbots[1:].values + topbots[:1].values)/2
dz=topbots[:-1].values - topbots[1:].values
layer = np.arange(len(z))+1
#%%Calculate times
start_year = 1999
times = flopy.utils.binaryfile.UcnFile(paths_ucn[0]).times
years = start_year + np.round(np.array(times)/365.25, 0)
years = [cftime.DatetimeProlepticGregorian(year, 1, 1) for year in years]

#%%Calculate proxy coordinates
#These are only used with concatenation
x = [idrange2coord(i, j+1, dcell) for i, j in df_mt[[13, 14]].values]
y = [idrange2coord(i, j+1, dcell) for i, j in df_mt[[15, 16]].values]

#%%Create DataArrays and combine these into a big one.
dims = ("time", "layer", "y", "x")

das_s = [(s, [xr.DataArray(ucn, dims = dims, 
                         coords={"x" : x[i], "y": y[i], "layer": layer, "time": years}
                         ) for i, ucn in ucn_data[s].items()]) for s in ucn_data.keys()]

das_h = [xr.DataArray(hds, dims = dims,
                      coords={"x" : x[i], "y": y[i], "layer": layer, "time": years}
                      ) for i, hds in enumerate(hds_data)]

da_all = [combine_all(das) for s, das in das_s]
da_all = xr.concat(da_all, dim="species").assign_coords(species = list(ucn_data.keys()))

da_hds = combine_all(das_h)
da_hds = da_hds.where(da_hds<1e20)

#Override coordinates with coordinates dummy idf, because we cannot infer from the UCN and .hds file
#which rows/columns had been removed during dropna when creating the model
da_all = da_all.assign_coords(x=dummy.x, y=dummy.y)
da_hds = da_hds.assign_coords(x=dummy.x, y=dummy.y)

#%%Assign some extra metadata for iMOD-python
da_all = add_metadata(da_all, z, dz, dcell)
da_hds = add_metadata(da_hds, z, dz, dcell)
da_all = da_all.rename("conc")

da_all = xr.where(np.isfinite(da_all), da_all, -9999)
da_hds = da_hds.where(da_hds<1e20, 1e30)

#%%Save inits
shead = da_hds.isel(time=-1)
sconc = da_all.isel(time=-1)

idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "bas", "head"), shead)
idf.save(os.path.join(modelfol, "..", mname+"{:02d}".format(mod_nr+1), "btn", "conc"), sconc, 
         pattern = r"{name}_{time:%Y%m%d%H%M%S}_c{species}_l{layer}{extension}")

#%%Split species into seperate variables
#Can only provide name
ds = da_all.to_dataset(dim="species")
ds = ds.rename(name_dict = dict([(s, "conc{}".format(s)) for s in np.unique(species)]))

#%%Deal with NaNs and assign head to dataset
ds["head"] = da_hds.where(da_hds<1e20)
ds["fhead"] = calc_fresh_water_head(ds["head"], ds["conc1"])

ds["head"]  = xr.where(np.isfinite(ds["head"]),  ds["head"],  1e30)
ds["fhead"] = xr.where(np.isfinite(ds["fhead"]), ds["fhead"], 1e30)

#%%Save to netcdf
ds.time.encoding["units"]    = r"days since 1999-01-01"
ds.time.encoding["calendar"] = r"proleptic_gregorian"

ds = ds.swap_dims({"layer" : "z"})

ds.to_netcdf(os.path.join(modelfol, "results_{:03d}.nc".format(mod_nr)),
                 unlimited_dims=["time"])