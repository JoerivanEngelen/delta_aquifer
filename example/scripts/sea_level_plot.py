# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 10:57:19 2019

@author: engelen
"""

import os
from delta_aquifer import boundary_conditions as bc
from pkg_resources import resource_filename
import numpy as np

#%%Path management
datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))
sl_curve = os.path.join(datafol, "spratt2016.txt")

figfol = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "example", "scratch_figures")))
os.makedirs(figfol, exist_ok=True)

#%%

# Time discretization
ts = np.array([30000, 25000, 20000, 15000, 13000, 12000, 11000, 10000, 9000,
               8000, 7000, 6000, 5000,  4000, 3000, 2000, 1000, 0])

ts = np.concatenate((np.arange(125000, 30000, -8000), ts))  #Add Pleistocene

#Transform to ka
ts = ts / 1000.

sea_level = bc.get_sea_level(sl_curve, ts, qt="50%", figfol=figfol, ext=".pdf")