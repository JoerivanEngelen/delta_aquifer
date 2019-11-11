# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:31:17 2019

@author: engelen
"""

import xarray as xr
from scipy import interpolate
from scipy.ndimage.morphology import binary_opening
import numpy as np
import os
if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt

#%%Helper functions
def _griddata(ip, val, targgrid):
    if val.ndim == 1:
        ip.values = np.expand_dims(val, 1)
    else:
        ip.values = val
    return(np.reshape(ip(tuple([i.flatten() for i in targgrid])), targgrid[0].shape))

def _clip_f(arr, sel, a_min=None, a_max=None):
    """Wrapper for clipping with fancy indexing"""
    def check(a_m, sel):
        if a_m is None:
            return(a_m)
        else:
            return(a_m[:, sel])
    
    arr[:, sel] = np.clip(arr[:, sel], check(a_min, sel), check(a_max, sel))
    return(arr)

def _ellipse(phi, D, phis):
    return(np.multiply.outer(D/phi,np.sqrt(phi**2-phis**2)))

def _cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def _pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def _calc_z_onshoreshelf(a, alpha, beta, apex, L_ab, x):
    z = np.zeros(x.shape)
    a_cells = x <= (a*L_ab)
    b_cells = x  > (a*L_ab)
    z[a_cells] = apex - (np.tan(alpha)) * x[a_cells]
    z[b_cells] = 0. - (np.tan(beta)) * x[:np.count_nonzero(b_cells)]
    return(z)

def _calc_z_slope(gamma, z_slope_end, x):
    return(z_slope_end + np.tan(gamma) * x[::-1])

def _calc_apex(a, alpha, L):
    return(L * a * np.tan(alpha))

def _calc_shelf_edge(b, beta, L):
    return(-L * b * np.tan(beta))

def _calc_1D_top_bot(a, alpha, beta, gamma, L, D, dD, x): 
    top = np.zeros(x.shape)
    
    L_ab_old = 0
    L_ab = L
    
    #Iteratively create crosssection as we do not know where the slope and
    #the coastal shelf will intersect
    while np.isclose(L_ab, L_ab_old) == False:
        L_ab_old = L_ab
        apex = _calc_apex(a, alpha, L_ab)
        top_ab = _calc_z_onshoreshelf(a, alpha, beta, apex, L_ab, x)
        
        bot = _calc_bot(apex, a, D, dD, L, x)
        top_slope = _calc_z_slope(gamma, bot[-1], x)
        
        ab_cells = top_ab < top_slope
        
        L_ab = x[ab_cells][-1]
    
    c_cells = ~ab_cells
    
    top[ab_cells] = top_ab[ab_cells]
    top[c_cells] = top_slope[c_cells]
    
    L_a = L_ab * a
    
    return(top, bot, c_cells, L_a)

def _calc_bot(apex, a, D, dD, L, x): 
    D_apex = apex - D*dD
    slope = (D_apex + D)/(L*a)
    return(D_apex - x*slope)

def _get_pizza_slice(d1, length, phis):
    idx_slice = d1["rho"] <= length
    return(_pol2cart(*np.meshgrid(d1["rho"][idx_slice], phis)))

def _cake_cuts(N, phi, rhos, f_offset = 0.):
    """Create N equally spaced lines through the delta.
    Similar to how you would subdivide a slice of cake.
    
    Parameters
    ----------
    N : int
        number of cuts
    
    phi : float
        phi
    
    rhos : array
        array with all the values of rho used in the polar discretization
    
    f_offset : float
        factor from zero to unity that is multiplied with the spacing to get an
        offset
    
    """
    
    assert(0. <= f_offset <= 1.0)
    
    phi_cuts = np.linspace(-phi/2, phi/2, num=N+2)[1:-1]
    offset = phi/(N+1)/2 * f_offset
    phi_cuts += offset
    
    #CHECK: If this leads to wide enough diagonal lines.
    x_cuts, y_cuts = _pol2cart(*np.meshgrid(rhos, phi_cuts))
    return(x_cuts, y_cuts)

def _rho_min_max(interface, rho):
    rho_ct = np.where(~np.isnan(interface), rho, 0)
    rho_max = np.max(rho_ct)
    cutoff = rho_max - 0.1 * rho_max
    rho_min = np.nanmin(np.where(np.max(rho_ct, axis=1) > cutoff, np.max(rho_ct, axis=1), np.nan))
    return(rho_min, rho_max)

#%%Main functions: Geometry
def create_crosssection(a, alpha, beta, gamma, dD, D, L, n_inp):
    d1 = {}
    
    d1["rho"] = np.linspace(0, L, num=n_inp)
     
    d1["top"], d1["bot"], d1["slope"], L_a = _calc_1D_top_bot(a, alpha, beta, gamma, 
                                          L, D, dD, d1["rho"])

    return(d1, L_a)

def create_2D_top_bot(phi, d1, n_inp):
    d2 = {}

    phis = np.linspace(-phi/2, phi/2, num=n_inp-20)
    d2["rho"], d2["phi"] = np.meshgrid(d1["rho"], phis)
    d2["phi_nan"] = np.copy(d2["phi"])
    d2["x"], d2["y"] = _pol2cart(d2["rho"], d2["phi"])
    
    ##Create tops
    d2["tops"], _ = np.meshgrid(d1["top"], phis)

    ##Fit ellipse for bottom
    d2["bots"]= _ellipse(phi/2, d1["bot"], phis).T
    
    ##Clip off weird part where tops and bot intersect
    d2["tops"] = _clip_f(d2["tops"], d1["slope"], a_min = d2["bots"])
    d2["bots"] = _clip_f(d2["bots"], ~d1["slope"], a_max = d2["tops"])
    nan_idx = np.isclose(d2["tops"], d2["bots"])
    nan_idx[:, ~d1["slope"]] = False
    for key, arr in d2.items():
        if key == "phi_nan":
            arr[~nan_idx] = np.nan
        else:
            arr[nan_idx] = np.nan
    return(d2, nan_idx, phis)

def get_targgrid(dx, dy, L, phi_max):
    x_edge, y_edge = _pol2cart(L, phi_max)

    x_max = np.round(L, -3) + 0.5*dx
    x_out = np.arange(0.5 * dx, x_max+1, dx)
    y_max = np.round(y_edge, -3) + 0.5*dy
    y_out = np.arange(-y_max, y_max+1, dy)
    return(np.meshgrid(x_out, y_out))

def pol2griddata(poldata, nan_idx, griddata, key_prep = ""):
    #Remove nans,which broadcasts to 1 dimensional array.
    poldata_no_nan=dict([(key, arr[~nan_idx]) for key, arr in poldata.items()])
    
    vrbls = [i for i in poldata.keys() if not i in ["x", "y", "phi", "phi_nan", "rho"]]
    #Explicitly create NDINterpolator once so QHULL has to be called only once, _griddata then just overrides the values.
    ip = interpolate.LinearNDInterpolator(np.vstack((poldata_no_nan["x"],poldata_no_nan["y"])).T, poldata_no_nan[vrbls[0]])
    for vrbl in vrbls:
        key = key_prep+vrbl
        griddata[key] = _griddata(ip, poldata_no_nan[vrbl], (griddata["x"], griddata["y"]))
    return(griddata)

def get_edges(ibound, top, z_shelf):
    dz = ibound.z[1] - ibound.z[0]
    is_shelf = (top>=z_shelf)
    edges = (
            (
            #Only select top layer of the shelf
            is_shelf & (ibound.z>(top-dz))
            ) | (
            #Select the coastal slope + the part just underneath the shelf break (roll(x=-1))
            ~is_shelf.roll(x=-1, roll_coords=False) | ~is_shelf
                        #Cells need to be active
                    )
            ) & ibound
    
    edges_z = edges.where(edges==1) * edges.z
    edges_x = edges.where(edges==1) * edges.x
    
    edges = xr.where(
        (
                edges.x == edges_x.max(dim="x")
        ) | (
                edges.z == edges_z.max(dim="z")
                ), edges, 0)
        
    edges = edges.transpose("z", "y", "x")
    return(edges)


def pol_to_d2_grid(dx, dy, phi, d2, nan_idx, d2_conf, nan_conf):
    d2_grid = {}
    
    d2_grid["x"], d2_grid["y"] = get_targgrid(dx, dy, 
           np.max(d2["x"][~nan_idx]), phi/2)
    d2_grid["rho"], d2_grid["phi"] = _cart2pol(d2_grid["x"], d2_grid["y"])
    d2_grid = pol2griddata(d2, nan_idx, d2_grid)
    if d2_conf is not None:
        d2_grid = pol2griddata(d2_conf, nan_conf, d2_grid, "conf_")
    return(d2_grid)

def dict_to_ds(dic, dims, coords):
    ds = xr.Dataset(
            dict([(key, (dims, value)) for key, value in dic.items()\
                  if key not in dims]),
            coords = coords)    
    return(ds)

def ds_d2_grid(d2_grid):
    coords = {"x" : d2_grid["x"][0, :],
          "y" : d2_grid["y"][:, 0]}
    ds = dict_to_ds(d2_grid, ["y", "x"], coords)
    ds = ds.dropna("y", "all")
    #For some reason this also clips off the first row where y["all"]!=np.nan??
    return(ds)

def ds_d1(d1):
    return(dict_to_ds(d1, ["rho"], {"rho" : d1["rho"]}))
    
def create_3d_grid(d2_grid, d1, nz):
    
    #TODO: This selection should be actually done when determining d1, 
    #but currently results in 
    to_keep = d1["top"] > d1["bot"]
    
    z = np.linspace(np.min(d1["bot"][to_keep]), 
                    np.max(d1["top"][to_keep]), num=nz)

    d3 = xr.Dataset(coords = dict(d2_grid.coords)).assign_coords(z=z)
    d3["IBOUND"] = xr.where((d3.z<d2_grid["tops"])&(d3.z>d2_grid["bots"]), 1, 0)
    #Brush away small isolated cells.
    d3["IBOUND"] = xr.full_like(d3["IBOUND"],
      #Maybe first allocate larger array to imitate mode="mirror" of convolvend
      binary_opening(d3["IBOUND"], structure=np.ones((2,2,2))).astype(np.int)
      )
    return(d3)

#%%Main functions: Geology
def get_relative_clay_depths(SM, n_clay):
    frac_sand = (1-SM) / (n_clay+1) #+1 because there also is a confining layer
    frac_clay = SM/(n_clay+1)
    
    rel_clay_depths = np.cumsum([frac_clay, frac_sand]*(n_clay+1))
    return(rel_clay_depths)

def create_confining_layer(clay_conf, d2, d1, phis, L_a, frac, n_inp):
    d2_clay = {}
    x_ons, y_ons = _get_pizza_slice(d1, L_a, phis)
    
    x_offset=L_a * (1-clay_conf)
    rho_clay, phi_clay = _cart2pol(x_ons[0][-1] - x_offset, y_ons[0][-1])
    rhos_clay = np.linspace(0, rho_clay, num=n_inp)
    phis_clay = np.linspace(-phi_clay, phi_clay, num=n_inp-20)
  
    topf = interpolate.interp1d(d1["rho"], d1["top"])
    botf = interpolate.interp1d(d1["rho"], d1["bot"])
    d_clay_coast = (botf(rhos_clay[-1]+x_offset) * frac)
    
    depth_conf =  np.linspace(0, d_clay_coast, num=n_inp)

    d2_clay["tops"] = np.zeros((n_inp-20, n_inp)) + topf(rhos_clay+x_offset) #Hack to get 2D tops of clay layer
    d2_clay["bots"] = d2_clay["tops"] - _ellipse(phi_clay, depth_conf, phis_clay).T
    
    d2_clay["x"], d2_clay["y"] = _pol2cart(*np.meshgrid(rhos_clay, phis_clay))
    d2_clay["x"] += x_offset
    
    #deactivate clay_layer past coastline
    nan_clay = np.greater(d2_clay["x"].T, np.nanmax(x_ons, axis = 1)).T
    for key, arr in d2_clay.items():
        arr[nan_clay] = np.nan
    return(d2_clay, nan_clay)

def create_clayer(frac, d1, d2, phis, phi, a_real, bot_offset=0):
    idx = int(d1["top"].shape[0] * a_real)
    base_slope = interpolate.interp1d(np.arange(idx), d1["top"][:idx], fill_value="extrapolate")
    base_slope = base_slope(np.arange(d1["top"].shape[0]))
    
    depth_clay = base_slope  - (base_slope-d1["bot"])*(frac)

    clayer = (_ellipse(phi/2, depth_clay, phis).T + depth_clay)/2.
    in_prism=(clayer<(d2["tops"]+bot_offset)) & (clayer>d2["bots"])
    clayer[~in_prism] = np.nan
    return(clayer)

def create_clayers(fracs, d1, d2, phis, phi, a_real, n_clay):
    rho_min, rho_max = {}, {}
    
    for i in range(n_clay):
        d2["ct%d"%i] = create_clayer(fracs[1::2][i], d1, d2, phis, phi, a_real)
        d2["cb%d"%i] = create_clayer(fracs[2::2][i], d1, d2, phis, phi, a_real)
        
        rho_min["ct%d"%i], rho_max["ct%d"%i]  = _rho_min_max(d2["ct%d"%i], d2["rho"])
        rho_min["cb%d"%i], rho_max["cb%d"%i]  = _rho_min_max(d2["cb%d"%i], d2["rho"])
        
        #Create edges at shoreline, so they extent fully to the sea
        d2["ct%d"%i] = np.where(np.isnan(d2["ct%d"%i]) & (d2["cb%d"%i] < d2["tops"]), d2["tops"], d2["ct%d"%i])
        #Find maximum the clay layer is allowed to be. Used to fix interpolation artifacts in finish_clayer_grid
        d2["ctmax%d"%i] =  np.full_like(d2["ct0"], 1) * np.nanmax(d2["ct0"], axis=0)
        #Create bottom at the side edges
        d2["cb%d"%i] = np.where(np.isnan(d2["cb%d"%i]) & (d2["ct%d"%i] > d2["bots"]), d2["bots"], d2["cb%d"%i])
        #Find minimum the clay layer is allowed to be. Used to fix interpolation artifacts in finish_clayer_grid
        d2["cbmin%d"%i] =  (np.full_like(d2["cb%d"%i], 1).T * np.nanmin(d2["cb%d"%i], axis=1)).T
    
    return(rho_min, rho_max, d2)

def finish_clayer_grid(d2_grid, n_clay, rho_min, rho_max):
    """Final touch to erase errors caused by the interpolation to a grid
    """

    for i in range(n_clay):
        #Make sure empty corners of clay layers are filled
        corner_top = ((d2_grid["rho"] < rho_max['ct%d'%i]) & (d2_grid["rho"] > rho_min['ct%d'%i]) & np.isnan(d2_grid["ct%d"%i]))
        corner_bot = ((d2_grid["rho"] < rho_max['cb%d'%i]) & (d2_grid["rho"] > rho_min['cb%d'%i]) & np.isnan(d2_grid["cb%d"%i]))
        
        d2_grid["ct%d"%i] = np.where(corner_top, np.minimum(d2_grid["tops"], d2_grid["ctmax%d"%i]), d2_grid["ct%d"%i])
        d2_grid["cb%d"%i] = np.where(corner_bot, d2_grid["cbmin%d"%i], d2_grid["cb%d"%i])
    return(d2_grid)

def create_paleo_channels(d2_ds, n_clay, N_pal, s_pal, phi, rhos): #Perhaps also include w_pal
    pal_masks = {}
    
    offset = np.zeros(n_clay)
    offset[::2] = s_pal
    
    for i in range(n_clay):
        pal = _cake_cuts(N_pal, phi, rhos, f_offset = offset[i])
        pal = [p.flatten() for p in pal]

        da = d2_ds.sel(x=xr.DataArray(pal[0], dims="foo"),
                       y=xr.DataArray(pal[1], dims="foo"), 
                       method="nearest")
        
        pal_masks[i] = ((d2_ds.x == da.x) & (d2_ds.y == da.y)).max(dim="foo")

    return(pal_masks)
    
def create_lith(d3, d2_grid, n_clay, pal_mask):
    d3["lith"] = d3["IBOUND"]
    
    #Confining clayer gets nr 2 as nr
    conf_nr = 2
    if "conf_tops" in d2_grid.keys():
        d3["lith"] = xr.where((d3.z<d2_grid["conf_tops"])&(d3.z>d2_grid["conf_bots"]), conf_nr, d3["lith"])
    
    pal_nr = 3
    for i in range(n_clay):
        #The other clay layer are assigned a number exceeding conf_nr
        clay_nr = 1 + pal_nr + i
        in_clayer = (d3.z<d2_grid["ct%d"%i])&(d3.z>d2_grid["cb%d"%i])
        d3["lith"] = xr.where(in_clayer & (~pal_mask[i]), clay_nr, d3["lith"])
        d3["lith"] = xr.where(in_clayer & (pal_mask[i]),  pal_nr , d3["lith"])
    return(d3)

def calc_clay_thicknesses(d2_grid, n_clay):
    if "conf_tops" in d2_grid.keys():
        d2_grid["conf_d"] = d2_grid["conf_tops"] - d2_grid["conf_bots"]
    for i in range(n_clay):
        d2_grid["cd%d"%i] = d2_grid["ct%d"%i] - d2_grid["cb%d"%i]
    return(d2_grid)

def create_Kh(d3, kh=0., kv_mar=0., f_kh_pal=0., ani=0., **kwargs):
    kh_mar = kv_mar*ani
    b = np.log10(kh_mar)
    a = np.log10(kh) - np.log10(kh_mar)
    
    kh_pal = 10**(a*f_kh_pal + b)
    
    d3["Kh"] = xr.zeros_like(d3["lith"])
    d3["Kh"] = xr.where(d3["lith"]==1, kh,      d3["Kh"])
    d3["Kh"] = xr.where(d3["lith"]==2, kh_mar,  d3["Kh"])
    d3["Kh"] = xr.where(d3["lith"]==3, kh_pal,  d3["Kh"])
    d3["Kh"] = xr.where(d3["lith"]>3,  kh_mar,  d3["Kh"])
    return(d3)

#%%Sedimentation/Erosion
def dynamic_confining_layer(d3, sea, t_max):
    sea = sea.max(dim="z")
    d3["lith"] = xr.where((d3["lith"] == 2) & (sea == 1)          , 1, d3["lith"]).astype(np.int64)
    d3["lith"] = xr.where((d3["lith"] == 2) & (d3["time"] > t_max), 1, d3["lith"]).astype(np.int64)
    d3 = d3.transpose("time", "z", "y", "x")
    return(d3)

def erosion_aquitards(d3, is_bc, bcs):
    #Choose maximum over z to keep aquitards active near the coastal slope.
    #Coastal shelf should only have one layer of bc per timestep.
    z_max_z = xr.where(is_bc, bcs.z, np.nan).max(dim="z")
    #Coordinates have to be monotonically increasing for ufunc
    z_max_z = z_max_z.assign_coords(time=z_max_z.time[::-1])
    
    #ignore nans seems not to be supported yet so have to use sentinel values    
    #also erode clay at edges that are just left untouched by our bcs
    sentinel=9999
    ffill = z_max_z.min(dim="y").fillna(sentinel) 
    z_max_z = z_max_z.fillna(ffill)

    ##Take minimum over time
    #This seems to work, I do not fully understand the docs though.
    z_min_t = xr.apply_ufunc(np.minimum.accumulate, z_max_z, 
                           input_core_dims  = [["x", "y"]], 
                           output_core_dims = [["x", "y"]])
    z_min_t = z_min_t.where(z_min_t<sentinel).assign_coords(time=z_min_t.time[::-1])
    
    d3["lith"] = xr.where((d3["lith"] >= 3) & (d3["z"]>z_min_t), 1, d3["lith"]).astype(np.int64)
    d3 = d3.transpose("time", "z", "y", "x")
    return(d3)

#%%Extra information
def add_topsystem(d2_grid):
    d2_grid["topsys"] = (np.nan_to_num(d2_grid["tops"]/d2_grid["tops"]) + 
               (d2_grid["tops"] >= 0)).astype(np.int8)
    
    if "conf_tops" in d2_grid.keys():
        d2_grid["topsys"] += np.nan_to_num(d2_grid["conf_tops"]/d2_grid["conf_tops"]).astype(np.int8)

    return(d2_grid)

#%%Master function
def get_geometry(a=None,  alpha=None,  beta=None,   gamma=None,   L=None, 
                 D=None,  dD=None,    phi=None, 
                 SM=None, n_clay=None,clay_conf=None, N_pal=None, s_pal=None,
                 dx=None, dy=None,    nz=None,  figfol=None, ncfol=None, **kwargs):

    n_inp = 200 #Discretization polar coordinates, not in actual model
    
    #Process central crosssection
    d1, L_a = create_crosssection(a, alpha, beta, gamma, dD, D, L, n_inp)
    
    #Create 3D top & bottom    
    d2, nan_idx, phis = create_2D_top_bot(phi, d1, n_inp)
    
    #Create confining layer
    frac_clay = SM/(n_clay+1)
    
    if clay_conf > 0.:
        d2_conf,nan_conf = create_confining_layer(clay_conf, d2, d1, 
                                                  phis, L_a, frac_clay, n_inp)
    else:
        d2_conf,nan_conf = None, None
    
    #Create clay layers
    rel_clay_depths = get_relative_clay_depths(SM, n_clay)
    rho_min, rho_max, d2 = create_clayers(rel_clay_depths, 
                                          d1, d2, phis, phi, L_a/L, n_clay)
    
    if figfol is not None:
        clayer_plot(d2, d2_conf, n_clay, L_a, L, figfol)
    
    #Convert to Cartesian regular grid
    d2_grid = pol_to_d2_grid(dx, dy, phi, d2, nan_idx, d2_conf, nan_conf)
    
    #Fill up last gaps in clayers
    d2_grid = finish_clayer_grid(d2_grid, n_clay, rho_min, rho_max)
    
    if figfol is not None:
        top_bot_grid_plot(d2_grid,figfol)
    
    #Convert to xarray
    d2_grid = ds_d2_grid(d2_grid)
    d1 = ds_d1(d1)

    #Mask paleo channels
    pal_mask = create_paleo_channels(d2_grid, n_clay, 
                                     N_pal, s_pal, phi, d1["rho"])

    #Add topsystem
    d2_grid = add_topsystem(d2_grid)
    if figfol is not None:
        _plot_imshow(d2_grid["topsys"], os.path.join(figfol, "topsystem_grid.png"))
    
    #Calculate clay thicknesses
    d2_grid = calc_clay_thicknesses(d2_grid, n_clay)

    #Create 3D 
    d3 = create_3d_grid(d2_grid, d1, nz)
    d3 = create_lith(d3, d2_grid, n_clay, pal_mask)
    
    z_shelf_edge = d1["top"][~d1["slope"]][-1]
    d3["edges"] = get_edges(d3["IBOUND"], d2_grid["tops"], z_shelf_edge)
    d3["topsys"], d3["tops"], d3["bots"] = d2_grid["topsys"], d2_grid["tops"], d2_grid["bots"]    

    #Save as netcdf
    if ncfol is not None:
        d3.to_netcdf(os.path.join(ncfol, "geo.nc"))

    #Assign layers after writing nc, as Paraview otherwise struggles with double coordinates
    layers = xr.DataArray(np.arange(len(d3.z))[::-1]+1, coords={"z":d3.z}, dims=["z"])
    d3 = d3.assign_coords(layer = layers)
    
    return(d3, d1, L_a)

#%%Plot functions
def clayer_plot(d2, d2_conf, n_clay, L_a, L, figfol, ext=".png"):
    id_x = int(d2["bots"].shape[1] * L_a /L / 2)
    
    slcs = [np.s_[:, id_x],np.s_[90, :]]
    dims = ["y", "x"]
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4))
    
    for i, slc in enumerate(slcs):
        xtop = xbot = d2[dims[i]][slc]
        axes[i].plot(xbot, d2["bots"][slc], label = "bot")
        axes[i].plot(xtop, d2["tops"][slc], label = "top")
        
        if dims[i] == "y":
            xval = d2["x"][slc][90]
            idx = np.searchsorted(d2["x"][90, :], xval)
            idx = np.s_[:, idx]
        else:
            axes[i].axvline(x=L_a, ls=":", color=".20")
            axes[i].axvline(x=L, color="k")
            idx = slc
        
        if d2_conf is not None:        
            xconf = d2_conf[dims[i]][idx]
            axes[i].fill(np.append(xconf[::-1], xconf), np.append(d2_conf["tops"][idx][::-1], d2_conf["bots"][idx]), label="conf")
        
        for j in range(n_clay):
            xtop = xbot = d2[dims[i]][slc]
            top = d2["ct%d"%j][slc]
            bot = d2["cb%d"%j][slc]
            
            xtop = xtop[~np.isnan(top)]
            xbot = xbot[~np.isnan(bot)]
            
            axes[i].fill(np.append(xtop[::-1], xbot), np.append(top[~np.isnan(top)][::-1], bot[~np.isnan(bot)]), label = "c%d"%j)
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figfol, "clayers%s" % ext))
    plt.close()

def _plot_imshow(da, fp):
    plt.imshow(da)
    plt.savefig(fp)
    plt.close()    

def top_bot_grid_plot(d2_grid,figfol):
    _plot_imshow(d2_grid["tops"], os.path.join(figfol, "top_grid.png"))
    _plot_imshow(d2_grid["bots"], os.path.join(figfol, "bot_grid.png"))
    