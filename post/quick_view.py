# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:31:11 2019

@author: engelen
"""

import imod
import xarray as xr
import matplotlib.pyplot as plt


concs = imod.idf.open(r"c:\Users\engelen\test_imodpython\synth_delta_test\test_small\results\conc\conc_45000101120000_l*.IDF")
#concs.where(concs==concs.max, drop=True)

concs = concs.sel(y=500)
concplot = concs.plot()
plt.gca().invert_yaxis()