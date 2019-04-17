# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:31:11 2019

@author: engelen
"""

import imod
import xarray as xr

concs = imod.idf.open(r"c:\Users\engelen\test_imodpython\synth_delta_test\test_small\results\conc\conc_44520101000000_l*.IDF")
concs.where(concs==concs.max, drop=True)