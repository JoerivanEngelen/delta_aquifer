# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:55:29 2020

@author: engelen
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

#%%Path management
path_Nile = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Nile_Delta\Data\Well_Data\Extractions_All.xlsx"
path_out = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\Data\depth_abstractions\Nile.png"
#%%Read
abs_Nile = pd.read_excel(path_Nile)

#%%Process data
abs_Nile["Screen_Mid"] = (abs_Nile["Screen_Top"] + abs_Nile["Screen_Bottom"])/2

abs_Nile = abs_Nile[abs_Nile["Screen_Mid"] < 1000.]
abs_Nile = abs_Nile.dropna(subset=["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"])

abs_Nile["depth (m)"] = pd.cut(abs_Nile["Screen_Mid"], np.arange(0, 321, 20), duplicates="drop")

Q_sum_depth = abs_Nile.groupby("depth (m)").sum()[["Irrigation", "Drinking", "Industrial", "Conjunctive_Use"]].reset_index()
Q_sum_depth = Q_sum_depth.set_index("depth (m)").sort_index(ascending=False)

mid_depth = [(a.left + a.right)/2 for a in Q_sum_depth.index]

#Scale
Q_sum_depth /= Q_sum_depth.values.sum()

Q_sum_depth = Q_sum_depth.reindex(mid_depth)
Q_sum_depth = Q_sum_depth.reindex(mid_depth)


#%%Exponential distribution
x = np.arange(0, 321, 20)
scale = 1/0.30
y = stats.expon.pdf(x/15, scale=scale)
#y = stats.gamma.pdf(x/20., 1., scale=scale)

#%%Plot
sns.set()
ax = Q_sum_depth.plot(kind="barh", stacked=True)
#ax = Q_sum_depth.plot(kind="bar", stacked=True)
#ax.set_ylabel("rel GW abstracted (-)")
plt.plot(y, np.arange(len(x))[::-1]-1, color="k")
ax.set_xlabel("rel GW abstracted (-)")

plt.tight_layout()

#%%Save
plt.savefig(path_out, dpi=300)