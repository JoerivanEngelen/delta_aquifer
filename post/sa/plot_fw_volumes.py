# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:55:25 2019

@author: engelen
"""
import os, sys
import xarray as xr

if os.name == "posix":
    import matplotlib 
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib.pyplot as plt

import seaborn as sns

#%%Path management
if len(sys.argv) > 1:
    fol  = sys.argv[1]
else:
    #Local testing on my own windows laptop
    fol = r"g:\synthdelta\results\test_output\synth_SD_i029_m24_7065039"

#%%Load data
frac_vols = xr.open_dataset(os.path.join(fol, "frac_vols.nc"))

#%%Plot
sns.set_style("darkgrid", {"axes.facecolor": ".4", "xtick.bottom": True})

lines = []
labels = ["fw", "bw", "sw", "ow"]
colors = [[0.5703125, 0.7734375, 0.87109375],
          [0.99609375, 0.87109375, 0.6015625],
          [0.9609375, 0.5625, 0.32421875],
          [0.1, 0.1, 0.1]]

for i, lab in enumerate(labels):
    lines.extend(frac_vols[lab].plot())
    lines[i].set_label(lab)
    lines[i].set_color(colors[i])

ax = plt.gca()

ax.invert_xaxis()
ax.xaxis.grid(False)

plt.title(None)
plt.ylabel("frac. volume (-)")
plt.ylim((0, 1))
plt.xlabel("time (ka)")
plt.legend()

plt.savefig(os.path.join(fol, "fw_volume.png"), dpi=300)

plt.close()