# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 10:41:38 2019

@author: engelen
"""

import os
from delta_aquifer import geometry
from pkg_resources import resource_filename
#%%Path management

figfol = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "example", "scratch_figures")))
os.makedirs(figfol, exist_ok=True)

#%%Specify parameters
pars = {}

a = 0.4
b = 0.4
D = 412
dD= 0.2
alpha = 0.00015
beta = 0.0016
gamma = 0.08
L = 200_000
phi = 1.57

SM = 0.3
n_clay = 3
clay_conf = 0.8

n_inp = 200 #Discretization polar coordinates, not in actual model

#%%

#Process central crosssection
d1 = geometry.create_crosssection(a, alpha, b, beta, gamma, dD, D, L, n_inp)

#Create 3D top & bottom    
d2, nan_idx, phis = geometry.create_top_bot(phi, d1, n_inp)

#Create confining layer
frac_clay = SM/(n_clay+1)    
d2_conf,nan_conf = geometry.create_confining_layer(clay_conf, d2, d1, 
                                          phis, L, a, frac_clay, n_inp)

#Create clay layers
rel_clay_depths = geometry.get_relative_clay_depths(SM, n_clay)
rho_min, rho_max, d2 = geometry.create_clayers(rel_clay_depths, 
                                      d1, d2, phis, phi, a, n_clay)

geometry.clayer_plot(d2, d2_conf, n_clay, a, b, L, figfol, ext=".pdf")
