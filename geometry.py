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

#%%functions
def griddata(x, y, val, targgrid, method="linear"):
    return(np.reshape(interpolate.griddata((x, y), val, tuple([i.flatten() for i in targgrid]), method=method), targgrid[0].shape))

def clip_f(arr, sel, a_min=None, a_max=None):
    """Wrapper for clipping with fancy indexing"""
    def check(a_m, sel):
        if a_m is None:
            return(a_m)
        else:
            return(a_m[:, sel])
    
    arr[:, sel] = np.clip(arr[:, sel], check(a_min, sel), check(a_max, sel))
    return(arr)

def ellipse(phi, D, x):
    return(np.multiply.outer(D/phi,np.sqrt(phi**2-x**2)))

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def calc_top(b, beta, a, alpha, gamma, L,  x): 
    z = np.zeros(x.shape)
    apex =  L * b * np.tan(beta)
    shelf = -L * a * np.tan(alpha)
    
    b_cells = x <= (b*L)
    a_cells = (x > (b*L)) & (x <= ((a+b)*L))
    c_cells = x > ((a+b)*L)
    
    z[b_cells] = apex - (np.tan(beta)) * x[b_cells] #onshore slope
    z[a_cells] = 0. - (np.tan(alpha)) * x[:np.count_nonzero(a_cells)]  #coastal shelf
    z[c_cells] = shelf - np.tan(gamma) * x[:np.count_nonzero(c_cells)]       #coastal slope
    
    return(z, c_cells)

def calc_bot(apex, b, D, dD, L, x): 
    D2 = apex - D*dD
    slope = (D2 + D)/(L*b)
    return(D2 - x*slope)

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

alpha = np.linspace(6e-4, 12e-4, num=lev)
beta  = np.linspace(0.75e-4, 1.5e-4, num=lev)

gamma = 5e-2 #FIXED will be corrected based on thickness.

#alpha = 8e-4 #in rads! For nile roughly: np.arctan(60/75000) = 8e-4 rad * 180/np.pi = 0.046 degrees
#beta = 1e-4 #in rads! For nile roughly: np.arctan(15/150000)

phi = np.linspace(0.125, 0.5, num=lev) * np.pi

clay_conf = np.linspace(0.2, 1., num=lev)

#%%For testing
a = a[1]
b = b[2]
D = D[2]
dD = dD[2]
alpha = alpha[2]
beta = beta[2]
phi = phi[2]

clay_conf = clay_conf[2]

#%%Process central crosssection
d2 = {}

gamma *= D/1000

d2["rho"] = np.linspace(0, L, num=ncell)
d2["top"], d2["slope"] = calc_top(b, beta, a, alpha, gamma, L, d2["rho"])
d2["bot"] = calc_bot(d2["top"][0], b, D, dD, L,  d2["rho"])

val_idx = d2["top"]>=d2["bot"]

for key, arr in d2.items():
    d2[key] = arr[val_idx]

##Plot
plt.plot(d2["rho"], d2["top"])
plt.plot(d2["rho"], d2["bot"])
plt.axvline(x=(a+b)*L, ls=":", color=".20")
plt.axvline(x=b*L, ls=":", color=".20")
plt.axvline(x=L, color="k")
plt.savefig(os.path.join(figfol, "2D_cross-section.png"))
plt.close()

#%%Create 3D top & bottom
d3 = {}

phis = np.linspace(-phi/2, phi/2, num=ncell-20)
d3["x"], d3["y"] = pol2cart(*np.meshgrid(d2["rho"], phis))

##Fit ellipse for bottom
d3["bots"]= ellipse(phi/2, d2["bot"], phis).T

##Create tops
d3["tops"], _ = np.meshgrid(d2["top"], phis)

##Clip off weird part where tops and bot intersect
d3["tops"] = clip_f(d3["tops"], d2["slope"], a_min = d3["bots"])
d3["bots"] = clip_f(d3["bots"], ~d2["slope"], a_max = d3["tops"])
nan_idx = np.isclose(d3["tops"], d3["bots"])
nan_idx[:, ~d2["slope"]] = False
for i, arr in d3.items():
    arr[nan_idx]=np.nan

##Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(d3["x"], d3["y"], d3["bots"], zorder=1)
ax.plot_surface(d3["x"], d3["y"], d3["tops"], zorder=20)
plt.savefig(os.path.join(figfol, "3D_top_bot.png"))
plt.close()

#%%Create confining layer
d3_clay = {}
idx_onshore = d2["rho"] < L*a
x_ons, y_ons = pol2cart(*np.meshgrid(d2["rho"][idx_onshore], phis))

rech=L * a * (1-clay_conf)
rho_clay, phi_clay = cart2pol(x_ons[0][-1] - rech, y_ons[0][-1])

rhos_clay = np.linspace(0, rho_clay, num=ncell)
phis_clay = np.linspace(-phi_clay, phi_clay, num=ncell-20)
d3_clay["x"], d3_clay["y"] = pol2cart(*np.meshgrid(rhos_clay, phis_clay))
d3_clay["x"] += rech

#deactivate clay_layer past coastline
nan_clay = np.greater(d3_clay["x"].T, np.nanmax(x_ons, axis = 1)).T
for key, arr in d3_clay.items():
    d3_clay[key] = arr[~nan_clay]

#%%Convert to Cartesian regular grid
dx, dy = 1000, 1000

x_edge, y_edge = pol2cart(L, phis[-1])

x_out = np.arange(0.5 * dx, L-0.5*dx+1, dx)
y_max = np.round(y_edge, -3) + 0.5*dy
y_out = np.arange(-y_max, y_max+1, dy)
targgrid=np.meshgrid(x_out, y_out)

for key, arr in d3.items():
    d3[key] = arr[~nan_idx]

topgrid = griddata(d3["x"], d3["y"], d3["tops"], targgrid)
botgrid = griddata(d3["x"], d3["y"], d3["bots"], targgrid)
claygrid = griddata(d3_clay["x"], d3_clay["y"], np.ones(d3_clay["y"].shape), targgrid, method="linear")

plt.imshow(topgrid)
plt.savefig(os.path.join(figfol, "top_grid.png"))
plt.close()

plt.imshow(botgrid)
plt.savefig(os.path.join(figfol, "bot_grid.png"))
plt.close()

#%%Create topsystem

topsys = np.nan_to_num(topgrid/topgrid) + np.nan_to_num(claygrid) + (topgrid >= 0)
plt.imshow(topsys)
plt.savefig(os.path.join(figfol, "topsystem_grid.png"))
plt.close()