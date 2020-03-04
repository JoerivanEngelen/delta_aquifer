# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:27:20 2020

@author: engelen
"""

from pkg_resources import resource_filename
import seaborn as sns
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import os

#%%Path management

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas")))

#%%Read inputs
par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)

#%%Select and add
mk = par.loc[153:161, ["Kh_aqf", "Kv_aqt"]]

mk["rating"] = [1, 1, 1, 1, 3, 2, 1, 3, 2]

#%%

class ScalarFormatterForceFormat(ticker.ScalarFormatter):
    def _set_format(self):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here

tick = ScalarFormatterForceFormat(useOffset=False, useMathText=True)
tick.set_powerlimits((0,0))

tx = [u"${}$".format(tick.format_data(x)) for x in np.sort(pd.unique(mk["Kv_aqt"]))]
ty = [u"${}$".format(tick.format_data(x)) for x in np.sort(pd.unique(mk["Kh_aqf"]))]

ax = sns.heatmap(mk.pivot("Kh_aqf", "Kv_aqt", "rating"), xticklabels = tx, yticklabels = ty)