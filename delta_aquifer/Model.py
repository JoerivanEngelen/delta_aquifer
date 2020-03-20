# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:01:24 2020

@author: engelen
"""

import numpy as np
import pandas as pd
from delta_aquifer import geometry, time_util, hydro_util
from delta_aquifer import boundary_conditions as bc
from delta_aquifer import initial_conditions as ic
import os

import imod
import xarray as xr
import cftime

class Synthetic(object):
    """Synthetic model object
    
    ----------
    Parameters
    ----------
    pars : dict
        dictionary with all the parameters necessary to run a model
    
    ts : list of floats
        list with timesteps in ka (note that there is some inconsistency with times in this package...)
    
    hclose : float
        hclose convergence criterium for SEAWAT
    
    rclose : float
        rclose convergence criterium for SEAWAT
    
    figfol : str
        path to folder where figures can be dumped
    
    ncfol : str
        path to folder where model input can be saved as NetCDF
        useful for post-processing.
    
    sl_curve : str
        path to txtfile with sea level curve
    
    -------
    Options
    -------
        init_salt : bool or str
            if str, should be the path to a NetCDF with inital conditions; 
            
            if bool, sets whether to start entirely initially salt or not. Note
            that False does not mean initially fresh, but instead the locatino of the
            interface is calculated based on the ghyben-herzberg approximation
        
        init_half2full : bool
            Sets whether to mirror read initial conditions along the y dimension
            Creates symmetrical initial conditions, useful to convert data 
            calculated with a halved model to initial conditions for a full model.
        
        half_model : bool
            Whether to half the model over the y-dimension, useful for 
            symmetrical models to save computational times
    """
    
    #%%Initialize
    def __init__(self, pars, ts, hclose, rclose, figfol, ncfol, sl_curve, 
                 half_model=False, init_salt=False, init_half2full=False, 
                 abstraction_path=None):
        self.hclose = hclose
        self.rclose = rclose
        self.ts = ts
        self.ncfol = ncfol
        self.pars = pars
        self.prepared = False
        self.wel = None
        
        d1, L_a = self.__initiate_geo(figfol)
        self.d1 = d1
        self.L_a = L_a
        self.__initiate_bcs(sl_curve, d1, L_a, figfol)
        self.__dynamize_geology()
        self.geo = geometry.create_Kh(self.geo, **self.pars)

        if half_model:
            self._half_model()
        self._clip_empty_cells()
        if type(init_salt) == str:
            self._read_initial_conditions(init_salt, half2full=init_half2full)
            self.assign_init_species=False
        else:
            self._create_initial_conditions(init_salt=init_salt)
            self.assign_init_species=True
    
        if abstraction_path is not None:
            if self.assign_init_species == False:
                sal = self.sconc.sel(species=1)
            else:
                sal = 0
            self.wel = bc.create_wells(abstraction_path, self.geo, self.bcs, 
                                       sal, figfol=figfol, **self.pars)
            
            self._assign_layers_to_wel()

            
    
    def __initiate_geo(self, figfol):
        """Initiate geometry & geology"""
        print("...initiating geometry & geology...")
        self.geo, d1, L_a = geometry.get_geometry(figfol=figfol, ncfol=None, **self.pars)

        self.topbot=bc._mid_to_binedges(self.geo["z"].values)[::-1]
        self.z=(self.topbot[:-1]+self.topbot[1:])/2
        
        return(d1, L_a)

    def __initiate_bcs(self, sl_curve, d1, L_a, figfol):
        """Initiate boundary conditions"""
        print("...initiating boundary conditions...")
        self.bcs, self.min_sea_level = bc.boundary_conditions(
                sl_curve, self.ts, self.geo, d1, conc_noise = 0.05,
                L_a=L_a,
                figfol=figfol, ncfol=None, **self.pars)

    def __dynamize_geology(self):
        """Create dynamic geology"""
        print("...dynamizing the geology...")
        #Holocene sedimentation model
        self.geo = geometry.dynamic_confining_layer(self.geo, self.bcs["sea"], self.pars["t_tra"])
        
        #Pleistocene erosion model
        is_bc = ((self.bcs["sea"]==1) | (self.bcs["river"]==1))
        self.geo = geometry.erosion_aquitards(self.geo, is_bc, self.bcs)

    def _assign_layers_to_wel(self):
        self.wel["layer"] = pd.cut(self.wel["well_z"], self.topbot[::-1],
                labels=self.geo.layer.values).astype(np.int64)

    def _half_model(self):
        """Cut model in half to speed up calculations"""
        print("...halving model...")
        self.geo = self.geo.sel(y=slice(0, self.geo.y.max()))
        self.bcs = self.bcs.sel(y=slice(0, self.geo.y.max()))
        
        #Since we cut the model in half, the middle channel (if existing) should half the conductance ~ y=0
        if (self.pars["N_chan"] % 2) == 1:
            y_loc = self.bcs.y.isel(y=0)
            self.bcs["riv_cond"].loc[dict(y = y_loc)] = self.bcs["riv_cond"].loc[dict(y = y_loc)]/2

    def _half_da_to_full(self, da):
        """Convert halved da back to a symmetrical full da along the y dimension
        """
        return(xr.concat([da, self._mirror_dims(da, "y")], dim="y"))
        
    def _clip_empty_cells(self):
        """Cut off unused x and y cells 
        otherwise writing the initial conditions for the next model run is 
        problematic due to the RCB algorithms completely leaving out unsused rows and columns
        """
        print("...clipping off empty cells...")
        active_plan = self.geo["IBOUND"].max(dim="z")
        active_plan = active_plan.where(active_plan==1.).dropna("x", how="all").dropna("y", how="all")
        
        self.geo = self.geo.sel(x=active_plan.x, y=active_plan.y)
        self.bcs = self.bcs.sel(x=active_plan.x, y=active_plan.y)

    def _read_initial_conditions(self, init_salt, half2full=False):
        """Read initial conditions, init_salt should be a path
        """
        assert(type(init_salt)==str)
        
        inits = xr.open_dataset(init_salt).isel(time=-1).drop("time")
        if half2full:
            inits = self._half_da_to_full(inits)
        
        self.shd = inits["head"]
        self.sconc = xr.concat(
                [inits["conc1"], inits["conc2"]], 
                "species").assign_coords(species=[1,2])

    def _create_initial_conditions(self, init_salt=False):
        print("...creating initial conditions...")
        assert(type(init_salt)==bool)
        if init_salt:
            salt_depth = 1000. #Set to depth at which salt begins way above 
                              #surface level so that everything is definitely salt
        else:
            salt_depth = self.min_sea_level #This was used in the global sensitivity analysis
        
        if salt_depth <= self.geo.z.min():
            self.approx_init = False
        else:
            self.approx_init = True
        
        self.shd, self.sconc = ic.get_ic(self.bcs, self.geo, approx_init=self.approx_init, 
                               deep_salt=salt_depth, **self.pars)

#%%Prepare functions
   
    def _calculate_dimensionless_numbers(self, pars, ncfol):
        print("...calculating dimensionless numbers...")
        dens_f, dens_s = ic.c2dens(pars["C_f"]), ic.c2dens(pars["C_s"])
        dimless = pd.DataFrame(np.array([hydro_util.rayleigh(
                        dens_f, dens_s, pars["H_b"], pars["D_m"], kv
                        ) for kv in [pars["Kh_aqf"]/pars["Kh_Kv"], pars["Kv_aqt"]]]),
                        index=["Ra", "Ra_mar"], columns=["value"])
            
        dimless.to_csv(os.path.join(ncfol, "dimless.csv"))

    def _prepare_time(self, max_sublen=8000):
        #Time management
        print("...preparing times...")
        self.start_year = 1999 #Must be minimum 1900 for iMOD-SEAWAT
        t_kyear = -1 * (self.ts * 1000 - self.ts[0] * 1000)
        
        #Adapt the times of bcs and geo, this is only used for save_inputs
        #In the writing process time will be overrided again by sub_ts
        time_util.num2date_ds(t_kyear[-len(self.bcs.time):], self.bcs, self.geo)
        
        
        sub_ts, sub_ends, sub_splits = time_util.subdivide_time(
                t_kyear[-(len(self.bcs.time)+1):], max_sublen
                )
        self.time_sub = {"ts_sub" : sub_ts,
                         "ends" : sub_ends,
                         "splits": sub_splits}

    def _save_inputs(self):
        print("...saving inputs...")
        self.bcs["lith"] = self.geo["lith"].astype(np.float64)
        self.bcs.to_netcdf(os.path.join(self.ncfol, "bcs.nc"))
        self.bcs = self.bcs.drop(["lith"])

    def _rev_dims(self, da, *dims):
        """Reverse dims, 
        alternative to all the slow sortby actions that slowed this package down previously
        """
        
        kwargs = dict([(dim, da[dim][::-1]) for dim in dims])
        return(da.reindex(**kwargs))

    def _mirror_dims(self, da, *dims):
        kwargs = dict([(dim , da[dim]*-1) for dim in dims])
        mirror = da.assign_coords(**kwargs)
        return(self._rev_dims(mirror, *dims))

    def revert_dim_for_iMOD(self, da):
        dx = da.x[1]-da.x[0]
        dy = da.y[1]-da.y[0]
        dlayer = da.layer[1] - da.layer[0]
        
        revert = ["x", dx < 0,], ["y", dy > 0,],["layer", dlayer < 0]
        args = [dim for dim, rev in revert if rev]
        return(self._rev_dims(da, *args))

    def _swap_and_reverse_dims(self):
        """extra processing to make iMOD-python accept these DataArrays
        and follow the .IDF format.
        """
        print("...swapping dimensions and reversing...")
        
        swap = {"z" : "layer"}
        
        self.geo = self.revert_dim_for_iMOD(self.geo.swap_dims(swap).drop("z"))
        self.bcs = self.revert_dim_for_iMOD(self.bcs.swap_dims(swap).drop("z"))
        self.shd = self.revert_dim_for_iMOD(self.shd.swap_dims(swap).drop("z"))
        
        if self.sconc.size > 1: #Can be a scalar
            self.sconc = self.revert_dim_for_iMOD(self.sconc.swap_dims(swap).drop("z"))

    def _combine_bcs(self):
        """Combine river and sea, as in the end we put both in the GHB anyway.
        """
        print("...combining boundary conditions...")
        sea = (self.bcs["sea"]==1)
        self.bcs["heads"] = xr.where(sea, self.bcs["sea_level"], self.bcs["riv_stage"])
        self.bcs = self.bcs.drop("riv_stage")
        self.bcs["conc"]  = xr.where(sea, self.bcs["sea_conc"] , self.bcs["riv_conc"])
        self.bcs = self.bcs.drop(["riv_conc", "sea_conc"])
        self.bcs["cond"]  = xr.where(sea, self.bcs["sea_cond"] , self.bcs["riv_cond"])
        self.bcs = self.bcs.drop(["riv_cond", "sea_cond"])
        

    def _add_species(self, assign_init_species=True):
        """Add extra species to keep track of the initial salt"""
        print("...adding species...")
        self.species = [1,2]
        self.bcs["conc"] = self.bcs["conc"].expand_dims(species=self.species)
        self.bcs["conc"] = xr.where(
                ((self.bcs["conc"].species == 2) & (np.isfinite(self.bcs["conc"].sel(species=2)))), 
                0., self.bcs["conc"])
        
        if self.assign_init_species:
            self.sconc = self.sconc.expand_dims(species=self.species)
            self.sconc = xr.where(((self.sconc.species == 2) & (self.sconc.sel(species=2)>1.0)), 1.0, self.sconc)
        
        self.bcs["rch_conc"] = xr.DataArray(data=[self.pars["C_f"]]*len(self.species), 
                                 coords=dict(species=self.species), dims=["species"])
    
    def prepare(self, max_sublen=8000):
        """Prepare all the necessary steps to have a working model and
        save the bcs for post-processing.
        
        max_sublen : int
            Maximum amount of years that should not be exceeded by a sub model.
            
            E.g. if set to 4000 years, splits one model of 12000 years into 
            3 consecutive models of 4000 years.
            
        """
        print("...starting preparation...")
        
        if max_sublen>8000:
            raise ValueError("The iMOD calendar cannot handle model extents longer than 8000 years")
        
        self._calculate_dimensionless_numbers(self.pars, self.ncfol)
        self._prepare_time(max_sublen=max_sublen)
        self._save_inputs()
        self._swap_and_reverse_dims()
        self._combine_bcs()
        self._add_species()
        self.prepared=True

#%%Sub model functions

    def _get_timesteps_submod(self, mod_nr):
        ts_mod = self.time_sub["ts_sub"][mod_nr]+self.start_year
        sub_end = self.time_sub["ends"][mod_nr]+self.start_year
        
        timesteps_mod = [cftime.DatetimeProlepticGregorian(t, 1, 1) for t in ts_mod]
        endtime=cftime.DatetimeProlepticGregorian(sub_end, 1, 1)        
        return(timesteps_mod, endtime)

    def __get_kh_submod(self, geo_mod):
        """Select for each submodel the Kh with the least aquitards and confining layer."""
        time_step_min_conf = geo_mod["lith"].where(geo_mod["lith"] == 2).sum(dim=["x", "y", "layer"]).argmin()
        time_step_min_aqtd = geo_mod["lith"].where(geo_mod["lith"] >= 3).sum(dim=["x", "y", "layer"]).argmin()
        
        kh = geo_mod["Kh"].isel(time=time_step_min_aqtd).drop("time")
        lith_min_conf = geo_mod["lith"].isel(time=time_step_min_conf).drop("time")
        lith_min_aqtd = geo_mod["lith"].isel(time=time_step_min_aqtd).drop("time")
            
        kh = xr.where((lith_min_conf != 2) & (lith_min_aqtd == 2), self.pars["Kh_aqf"], kh)
        return(kh)

    def __get_ic_submod(self, mod_nr, geo_mod):
        """Select correct initial conditions"""
        if mod_nr == 0:
            active=(geo_mod["IBOUND"]==1)
            starting_head = xr.where(active,   self.shd, -9999.0)
            starting_conc = xr.where(active, self.sconc, -9999.0)
            
        else:
            year = self.time_sub["ends"][mod_nr-1]+self.start_year
            year_str = cftime.DatetimeProlepticGregorian(
                    year, 1, 1).strftime("%Y%m%d%H%M%S")
            starting_head = "bas/head_{}_l?.idf".format(year_str)
            starting_conc = ["btn/conc_{}_c{}_l?.idf".format(year_str, specie) for specie in self.species]
            starting_conc = xr.DataArray(data=starting_conc, 
                                         coords=dict(species=self.species), 
                                         dims=["species"])
        
        return(starting_head, starting_conc)

    def _correct_wel_time(self, timesteps_mod):
            t_wel = (self.ts[:-1] + self.ts[1:])/2
            self.wel["time"] = self.wel["time"].replace(dict(zip(t_wel, timesteps_mod)))        
            assert(np.all(pd.unique(self.wel["time"])==timesteps_mod))

#%%Write

    def write_model(self, model_fol, mname, write_first_only=False, max_perlen=1000.):
        if self.prepared == False:
            raise ValueError("Model has not been prepared. Please call Synthetic.prepare() first")
        splitted_t = enumerate(zip(self.time_sub["splits"][:-1], self.time_sub["splits"][1:]))
        for mod_nr, (i_start, i_end) in splitted_t:
            print(".........processing model nr. {}..........".format(mod_nr))
            
            timesteps_mod, endtime = self._get_timesteps_submod(mod_nr)

            if self.wel is not None:
                self._correct_wel_time(timesteps_mod)
                wel = self.wel #TODO: Make sure wells can be split up as well
            
            bcs_mod, geo_mod = [da.isel(time=slice(i_start, i_end)).assign_coords(
                    time = timesteps_mod) for da in [self.bcs, self.geo]]
                    
            kh = self.__get_kh_submod(geo_mod)
            
            starting_head, starting_conc = self.__get_ic_submod(mod_nr, geo_mod)
            
            #Create model
            mname_sub = "{}_nr{:02d}".format(mname, mod_nr)
            
            m = imod.wq.SeawatModel(mname_sub, check=None)
            
            bottoms = xr.DataArray(self.topbot[1:], 
                                   {"layer": self.geo.layer}, 
                                   ("layer"))
            
            m["bas"] = imod.wq.BasicFlow(ibound=self.geo["IBOUND"].assign_coords(
                    dx=self.pars["dx"], dy=-self.pars["dy"]), 
                                         top=self.topbot[0], 
                                         bottom=bottoms, 
                                         starting_head=starting_head)
            
            m["lpf"] = imod.wq.LayerPropertyFlow(
                k_horizontal=kh, k_vertical=kh/self.pars["Kh_Kv"], 
                specific_storage=self.pars["S_s"],
                save_budget=True,
            )
            
            m["btn"] = imod.wq.BasicTransport(
                n_species=len(self.species),
                icbund=self.geo["IBOUND"], 
                starting_concentration=starting_conc, 
                porosity=self.pars["n"],
                inactive_concentration = -9999.0
            )
            
            m["adv"] = imod.wq.AdvectionTVD(courant=0.9)
            m["dsp"] = imod.wq.Dispersion(
                    longitudinal=self.pars["a_l"], 
                    ratio_horizontal=self.pars["trpt"],
                    ratio_vertical=self.pars["trpv"],
                    diffusion_coefficient=self.pars["D_m"]
            )
            
            m["vdf"] = imod.wq.VariableDensityFlow(density_concentration_slope=0.7143)
            
            m["ghb"] = imod.wq.GeneralHeadBoundary(head = bcs_mod["heads"],
                                                   conductance=bcs_mod["cond"],
                                                   density=ic.c2dens(bcs_mod["conc"].sel(species=1)), 
                                                   concentration=bcs_mod["conc"])
            
            if self.pars["R"] != 0:
                m["rch"] = imod.wq.RechargeHighestActive(rate=bcs_mod["rch"],
                                                         concentration=self.bcs["rch_conc"])
            if self.wel is not None:
                m["wel"] = imod.wq.Well(wel["name"].to_list(),
                                        wel["x"].to_list(), wel["y"].to_list(),
                                        wel["Q"].to_list(), layer=wel["layer"].to_list(),
                                        time=wel["time"].to_list())
            
            m["pksf"] = imod.wq.ParallelKrylovFlowSolver(
                                                         max_iter=1000, 
                                                         inner_iter=100, 
                                                         hclose=self.hclose, 
                                                         rclose=self.rclose, 
                                                         relax=1.00,
                                                         partition="rcb",
                                                         solver="pcg",
                                                         preconditioner="ilu",
                                                         deflate=False,
                                                         debug=False,
                                                         )
            
            m["pkst"] = imod.wq.ParallelKrylovTransportSolver(
                                                         max_iter=1000, 
                                                         inner_iter=30, 
                                                         cclose=1e-6,
                                                         relax=0.98,
                                                         partition="rcb",
                                                         solver="bicgstab",
                                                         preconditioner="ilu",
                                                         debug=False,
                                                         )
            
            m["oc"] = imod.wq.OutputControl(save_head_idf=True, 
                                            save_concentration_idf=True, 
                                            save_budget_idf=True)
            
            self.model = m
            
            n_timesteps_p1 = 10
            time_util.time_discretization(m, max_perlen, 
                                          endtime=endtime,
                                          n_timesteps_p1=n_timesteps_p1,
                                          timestep_multiplier=7.,
                                          max_n_transport_timestep=999_999,
                                          transport_initial_timestep=100.)
            
            directory = os.path.join(model_fol, mname, mname_sub)
            
            m.write(directory = directory, 
                    result_dir = os.path.join(directory, "results"), 
                    directstop=True)
            
            if (write_first_only==True) & (mod_nr == 0):
                break    