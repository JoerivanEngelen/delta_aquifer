# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:31:17 2019

@author: engelen
"""

import xarray as xr
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D #We have to import this to allow 3d plotting, otherwise projection = "3d" not recognized

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

def _calc_top(a, alpha, b, beta, gamma, L,  x): 
    z = np.zeros(x.shape)
    apex =  L * a * np.tan(alpha)
    shelf = -L * b * np.tan(beta)
    
    a_cells = x <= (a*L)
    b_cells = (x > (a*L)) & (x <= ((a+b)*L))
    c_cells = x > ((a+b)*L)
    
    z[a_cells] = apex - (np.tan(alpha)) * x[a_cells] #onshore slope
    z[b_cells] = 0. - (np.tan(beta)) * x[:np.count_nonzero(b_cells)]  #coastal shelf
    z[c_cells] = shelf - np.tan(gamma) * x[:np.count_nonzero(c_cells)]       #coastal slope
    
    return(z, c_cells)

def _calc_bot(apex, a, D, dD, L, x): 
    d1 = apex - D*dD
    slope = (d1 + D)/(L*a)
    return(d1 - x*slope)

#%%Main functions
def create_crosssection(a, alpha, b, beta, gamma, dD, D, L):
    d1 = {}

    gamma *= D/1000
    
    d1["rho"] = np.linspace(0, L, num=ncell)
    d1["top"], d1["slope"] = _calc_top(a, alpha, b, beta, gamma, L, d1["rho"])
    d1["bot"] = _calc_bot(d1["top"][0], a, D, dD, L,  d1["rho"])
    
    val_idx = d1["top"]>=d1["bot"]
    
    for key, arr in d1.items():
        d1[key] = arr[val_idx]

    return(d1)

def create_top_bot(phi, d1):
    d2 = {}

    phis = np.linspace(-phi/2, phi/2, num=ncell-20)
    d2["x"], d2["y"] = _pol2cart(*np.meshgrid(d1["rho"], phis))
    
    ##Fit ellipse for bottom
    d2["bots"]= _ellipse(phi/2, d1["bot"], phis).T
    
    ##Create tops
    d2["tops"], _ = np.meshgrid(d1["top"], phis)
    
    ##Clip off weird part where tops and bot intersect
    d2["tops"] = _clip_f(d2["tops"], d1["slope"], a_min = d2["bots"])
    d2["bots"] = _clip_f(d2["bots"], ~d1["slope"], a_max = d2["tops"])
    nan_idx = np.isclose(d2["tops"], d2["bots"])
    nan_idx[:, ~d1["slope"]] = False
    for i, arr in d2.items():
        arr[nan_idx]=np.nan
    return(d2, nan_idx, phis)

def create_confining_layer(clay_conf, d2, d1, phis, L, a, frac):
    d2_clay = {}
    idx_onshore = d1["rho"] < L*a
    x_ons, y_ons = _pol2cart(*np.meshgrid(d1["rho"][idx_onshore], phis))
    
    x_offset=L * a * (1-clay_conf)
    rho_clay, phi_clay = _cart2pol(x_ons[0][-1] - x_offset, y_ons[0][-1])
    rhos_clay = np.linspace(0, rho_clay, num=ncell)
    phis_clay = np.linspace(-phi_clay, phi_clay, num=ncell-20)
    
    topf = interpolate.interp1d(d1["rho"], d1["top"])
    botf = interpolate.interp1d(d1["rho"], d1["bot"])
    
    depth_conf = topf(rhos_clay+x_offset) + botf(rhos_clay+x_offset) * frac
    d2_clay["bots"] = -_ellipse(phi_clay, depth_conf, phis_clay).T
    d2_clay["tops"] = np.zeros(d2_clay["bots"].shape) + topf(rhos_clay+x_offset) #Hack to get 2D tops of clay layer
    
    d2_clay["x"], d2_clay["y"] = _pol2cart(*np.meshgrid(rhos_clay, phis_clay))
    d2_clay["x"] += x_offset
    
    #deactivate clay_layer past coastline
    nan_clay = np.greater(d2_clay["x"].T, np.nanmax(x_ons, axis = 1)).T
    for key, arr in d2_clay.items():
        arr[nan_clay] = np.nan
#        d2_clay[key] = arr[~nan_clay]
    return(d2_clay, nan_clay)

def get_targgrid(dx, dy, L, phi_max):
    x_edge, y_edge = _pol2cart(L, phi_max)

    x_out = np.arange(0.5 * dx, L-0.5*dx+1, dx)
    y_max = np.round(y_edge, -3) + 0.5*dy
    y_out = np.arange(-y_max, y_max+1, dy)
    return(np.meshgrid(x_out, y_out))

def create_clayer(frac, bot, d2, phis, phi, a_off = 0.):
    
#    depth_clay = d1["top"] + (d1["top"]+d1["bot"])*(frac)
    depth_clay = a_off + bot*frac
    clayer = (_ellipse(phi/2, depth_clay, phis).T + depth_clay)/2.
    in_prism=(clayer<d2["tops"]) & (clayer>d2["bots"])
    clayer[~in_prism] = np.nan
    return(clayer)

def pol2griddata(poldata, nan_idx, griddata, key_prep = ""):
    #Remove nans,which broadcasts to 1 dimensional array.
    for key, arr in poldata.items():
        poldata[key] = arr[~nan_idx]
    
    vrbls = [i for i in poldata.keys() if not i in ["x", "y"]]
    #Explicitly create NDINterpolator once so QHULL has to be called only once, _griddata then just overrides the values.
    ip = interpolate.LinearNDInterpolator(np.vstack((poldata["x"],poldata["y"])).T, poldata[vrbls[0]])
    for vrbl in vrbls:
        griddata[key_prep+vrbl] = _griddata(ip, poldata[vrbl], (griddata["x"], griddata["y"]))
    return(griddata)

#%%Plot functions
def clayer_plot(d2, d2_conf, figfol):
    slcs = [np.s_[:, 50],np.s_[90, :]]
    dims = ["y", "x"]
    fig, axes = plt.subplots(ncols=2, figsize=(8, 4))
    
    for i, slc in enumerate(slcs):
        xtop = xbot = d2[dims[i]][slc]
        axes[i].plot(xbot, d2["bots"][slc], label = "bot")
        axes[i].plot(xtop, d2["tops"][slc], label = "top")
        
        if dims[i] == "y":
            xval = d2["x"][slc][90]
            idx = np.searchsorted(d2_conf["x"][90, :], xval)
            idx = np.s_[:, idx]
        else:
            axes[i].axvline(x=(a+b)*L, ls=":", color=".20")
            axes[i].axvline(x=a*L, ls=":", color=".20")
            axes[i].axvline(x=L, color="k")
            idx = slc
            
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
    plt.savefig(os.path.join(figfol, "clayers.png"))
    plt.close()

#%%Path management
figfol=r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Modelinput\Figures"

#%%Parameters
#Morris parameters
lev = 4

#Model input parameters
ncell = 200

#Geometry input parameters
L = 200000 #FIXED Check all delta lengths, so if this value is representative. Also check leakage factors whether this long enough
b = np.linspace(0.3, 0.6, num=lev) #Check this for deltas
a = np.linspace(0.3, 0.6, num=lev) #Check this for deltas

D = np.array([70, 250, 500, 1000])
dD = np.linspace(0.2, 0.6, num=lev)

alpha  = np.linspace(0.75e-4, 1.5e-4, num=lev) #Check this for deltas
beta = np.linspace(6e-4, 12e-4, num=lev) #Check this for deltas

gamma = 5e-2 #FIXED will be corrected based on thickness.

#alpha = 8e-4 #in rads! For nile roughly: np.arctan(60/75000) = 8e-4 rad * 180/np.pi = 0.046 degrees
#beta = 1e-4 #in rads! For nile roughly: np.arctan(15/150000)

phi = np.linspace(0.125, 0.5, num=lev) * np.pi

clay_conf = np.linspace(0.2, 1., num=lev)

n_clay = np.linspace(0, 3, num=lev, dtype=int)
SM = 0.3 #FIXED FOR NOW ELSE np.linspace(0.1, 0.6, num=4)

#%%For testing
#Domain geometry
a = a[1]
b = b[1]
D = D[2]
dD = dD[2]
alpha = alpha[2]
beta = beta[2]
phi = phi[2]

#Internal geometry
clay_conf = clay_conf[2]
n_clay = n_clay[3]

#%%Process central crosssection
d1 = create_crosssection(a, alpha, b, beta, gamma, dD, D, L)

#%%Create 3D top & bottom
d2, nan_idx, phis = create_top_bot(phi, d1)

##Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(d2["x"], d2["y"], d2["bots"], zorder=1)
ax.plot_surface(d2["x"], d2["y"], d2["tops"], zorder=20)
plt.savefig(os.path.join(figfol, "3D_top_bot.png"))
plt.close()

#%%Create confining layer
frac_sand = (1-SM) / (n_clay + 1)
frac_clay = SM/(n_clay+1)

d2_conf,nan_conf = create_confining_layer(clay_conf, d2, d1, phis, L, a, frac_clay)

#%%Create clay layers
fracs = np.cumsum([frac_sand] + [frac_clay, frac_sand]*n_clay)

for i in range(n_clay):
    d2["ct%d"%i] = create_clayer(fracs[0::2][i], d1["bot"], d2, phis, phi)
    d2["cb%d"%i] = create_clayer(fracs[1::2][i], d1["bot"], d2, phis, phi)

##Plot
clayer_plot(d2, d2_conf, figfol)

#%%Convert to Cartesian regular grid
d2_grid = {}

dx, dy = 1000, 1000

d2_grid["x"], d2_grid["y"] = get_targgrid(dx, dy, L, phi/2)
d2_grid = pol2griddata(d2, nan_idx, d2_grid)
d2_grid = pol2griddata(d2_conf, nan_conf, d2_grid, "conf_")

plt.imshow(d2_grid["tops"])
plt.savefig(os.path.join(figfol, "top_grid.png"))
plt.close()

plt.imshow(d2_grid["bots"])
plt.savefig(os.path.join(figfol, "bot_grid.png"))
plt.close()

#%%Create topsystem
d2_grid["topsys"] = (np.nan_to_num(d2_grid["tops"]/d2_grid["tops"]) + 
           np.nan_to_num(d2_grid["conf_tops"]/d2_grid["conf_tops"]) + 
                                                  (d2_grid["tops"] >= 0)).astype(np.int16)

plt.imshow(d2_grid["topsys"])
plt.savefig(os.path.join(figfol, "topsystem_grid.png"))
plt.close()

#%%Create 3D 
d2_grid = xr.Dataset(dict([(key, (["y", "x"], value)) for key, value in d2_grid.items() if key not in ["x", "y"]]),
                     coords = {"x" : d2_grid["x"][0, :],
                               "y" : d2_grid["y"][:, 0]})

