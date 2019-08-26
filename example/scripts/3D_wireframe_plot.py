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
a = 0.4
b = 0.4
D = 412
dD= 0.2
alpha = 0.00015
beta = 0.0016
gamma = 0.08
L = 200_000
phi = 1.57

n_inp = 200 #Discretization polar coordinates, not in actual model

#%%Process
#Process central crosssection
d1 = geometry.create_crosssection(a, alpha, b, beta, gamma, dD, D, L, n_inp)

#Create 3D top & bottom    
d2, nan_idx, phis = geometry.create_top_bot(phi, d1, n_inp)

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