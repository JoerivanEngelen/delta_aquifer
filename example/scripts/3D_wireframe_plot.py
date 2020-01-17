# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:01:36 2019

@author: engelen

Plot 3D wireframe

"""

import matplotlib.pyplot as plt
#We have to import this to allow 3d plotting, otherwise 
#projection = "3d" not recognized
from mpl_toolkits.mplot3d import Axes3D 

from delta_aquifer import geometry
from pkg_resources import resource_filename

import os

#%%Path management

figfol = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "example", "scratch_figures")))
os.makedirs(figfol, exist_ok=True)

#%%Specify parameters
#Geometry
l_a = 0.4
H_b = 412
f_H = 0.2
alpha = 0.00015
beta = 0.0016
gamma = 0.08
L = 200_000
phi_f = 1.57

f_aqt = 0.3
l_conf = 0.8
N_aqt = 3

n_inp = 200 #Discretization polar coordinates, not in actual model

#%%Process
#Process central crosssection
d1, L_a = geometry.create_crosssection(l_a, alpha, beta, gamma, H_b, f_H, L, n_inp)

#Create 3D top & bottom    
d2, nan_idx, phis = geometry.create_2D_top_bot(phi_f, d1, n_inp)

#Create confining layer
frac_clay = f_aqt/(N_aqt+1)

if l_conf > 0.:
    d2_conf,nan_conf = geometry.create_confining_layer(l_conf, d2, d1, 
                                              phis, L_a, frac_clay, n_inp)
else:
    d2_conf,nan_conf = None, None

#Create clay layers
rel_clay_depths = geometry.get_relative_clay_depths(f_aqt, N_aqt)
rho_min, rho_max, d2 = geometry.create_clayers(rel_clay_depths, 
                                      d1, d2, phis, phi_f, L_a/L, N_aqt)

#%%Interpolate coastal slope
#You need to manually fiddle around with this depending on parameters
idx_strt_slope = 160

refine = slice(idx_strt_slope, 184)
rest = slice(None, idx_strt_slope)
plt_vars = ["x", "y", "tops"]

tp_slope = dict([(var, d2[var][:, refine]) for var in plt_vars])
tp_nope = dict([(var, d2[var][:, rest]) for var in plt_vars])

#%%Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ccount=50

ax.plot_wireframe(tp_slope["x"], tp_slope["y"], tp_slope["tops"], zorder=20, rcount=1, ccount=ccount, alpha=0.7, color="k")
ax.plot_wireframe(tp_nope["x"], tp_nope["y"], tp_nope["tops"], zorder=20, rcount=1, ccount=ccount, alpha=0.7, color="k")
ax.plot_wireframe(d2["x"], d2["y"], d2["bots"], zorder=1, rcount=1, alpha=0.3, ccount=ccount, color="k")
ax.set_axis_off()

plt.savefig(os.path.join(figfol, "3D_top_bot.pdf"))
plt.close()

#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ccount=50

ax.plot_surface(d2["x"], d2["y"], d2["cb2"], zorder=1, alpha=1,linewidth=0, rcount=1,  ccount=ccount, color="#512c14ff", shade=True)
ax.plot_surface(d2["x"], d2["y"], d2["ct2"], zorder=1, alpha=1, linewidth=0, rcount=1, ccount=ccount, color="#784421ff", shade=True)

ax.plot_surface(d2["x"], d2["y"], d2["cb1"], zorder=1, alpha=1, linewidth=0, rcount=1, ccount=ccount, color="#512c14ff", shade=True)
ax.plot_surface(d2["x"], d2["y"], d2["ct1"], zorder=1, alpha=1, linewidth=0, rcount=1, ccount=ccount, color="#784421ff", shade=True)

ax.plot_surface(d2["x"], d2["y"], d2["cb0"], zorder=2, linewidth=0, alpha=1, rcount=1, ccount=ccount, color="#512c14ff", shade=True)
ax.plot_surface(d2["x"], d2["y"], d2["ct0"], zorder=2, linewidth=0, alpha=1, rcount=1, ccount=ccount, color="#784421ff", shade=True)

ax.plot_surface(d2_conf["x"], d2_conf["y"], d2_conf["bots"], zorder=3, alpha=1, rcount=1, ccount=ccount, linewidth=0, color="#512c14ff", shade=True)
ax.plot_surface(d2_conf["x"], d2_conf["y"], d2_conf["tops"], zorder=3, alpha=1, rcount=1, ccount=ccount, linewidth=0, color="#c87137ff", shade=True)
                
ax.set_axis_off()

plt.savefig(os.path.join(figfol, "3D_clayers.pdf"))
plt.close()