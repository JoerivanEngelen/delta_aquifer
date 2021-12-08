# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:27:58 2020

@author: engelen
"""

import xarray as xr
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

from itertools import product

def prepare_for_plot(ds):
    ds = ds.where(ds > -9998.)
    ds = ds.assign_coords(x = ds["x"]/100000)
    ds = ds.assign_coords(z = ds["z"]/10)
    ds = ds.reset_coords(drop=True)
    ds = ds.clip(min=0., max=35.)
    return(ds)

#%%Path management
outf = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results\2D_cross"

folders = [#r"g:\synthdelta\results_SA_zip\synth_SD_i018_m24_7080926",
           r"g:\synthdelta\results_SA_zip\synth_SD_i147_m24_7114854",
           r"g:\synthdelta\results_SA_zip\synth_SD_i043_m24_7065053",
           r"g:\synthdelta\results_SA_zip\synth_SD_i049_m24_7087160",
           r"g:\synthdelta\results_SA_zip\synth_SD_i123_m24_7114121",
           r"g:\synthdelta\results_SA_zip\synth_SD_i140_m24_7114849"]

files = ["results_006.nc",
         "results_006.nc",
         "results_006.nc",
         "results_001.nc",
         "results_006.nc"]

files = [os.path.join(folder, f) for folder, f in zip(folders, files)]

#%%Classification
types = [3, 5, 4, 1, 2]

#%%Reorder
files = [x for _,x in sorted(zip(types,files))]

#%%Read and select
ds_ls = [xr.open_dataset(f)["conc1"].isel(time=-1) for f in files]
ds_2d = [ds.isel(y=-1) for ds in ds_ls]

#%%Prepare for plot
ds_2d = [prepare_for_plot(ds) for ds in ds_2d]

#%%Plot
colors = [
 [0.171875, 0.48046875, 0.7109375], 
 [0.5703125, 0.7734375, 0.87109375], 
 [0.8671875, 0.9375, 0.8125], 
 [0.99609375, 0.87109375, 0.6015625], 
 [0.9609375, 0.5625, 0.32421875], 
 [0.83984375, 0.09765625, 0.109375], 
 [0.62890625, 0.09765625, 0.109375], 
 # [0.421875, 0.09765625, 0.109375], 
]

# levels = [0, 1, 5, 10, 15, 20, 25, 30, 35]
levels = [0, 1, 5, 10, 15, 20, 25, 35]

#http://xarray.pydata.org/en/stable/generated/xarray.DataArray.plot.contourf.html
#ds_2d[0].plot.contourf(levels = levels, colors = colors)

#%%Setup Plot
widths = [10, 10, 1]
widths = [10, 10]

dark = ".1"
light = ".95"

style_dict = {"axes.facecolor": dark, 
               "axes.labelcolor" : dark,
               'text.color': light, 
               'figure.facecolor': dark, 
               'xtick.color': light,
               'ytick.color': light}

with sns.axes_style("dark", style_dict):
    with sns.plotting_context("paper", font_scale=0.7):
        fig = plt.figure(constrained_layout=True, figsize = (19/2.54/3, 11.5/2.54))
        gs = fig.add_gridspec(3,2, width_ratios=widths)
        
        
        plot_locs = list(product([0, 1, 2], [0, 1]))[:-1]
        axes = [fig.add_subplot(gs[loc]) for loc in plot_locs]
        cbar_ax = fig.add_subplot(gs[2, 1])
        cbar_axin = inset_axes(cbar_ax,
                            width="30%",  # width = 50% of parent_bbox width
                            height="90%",
                            loc='center left')
        cbar_ax.axis("off")
 
        for spine in cbar_axin.spines.values():
           spine.set_edgecolor(dark)       
 
        plot_cbars = [False, False, False, False, True]
        cbar_axes = [None, None, None, None, cbar_axin]
        plot_yaxes = [True, False] * 3
        plot_xaxes = [False, False, False, False, True]
        cbar_kwargs = {}
        cbar_kwargs = [None, None, None, None, cbar_kwargs]
        
        for ds, ax, plot_cbar, cbar_ax, plot_yaxis, plot_xaxis, cbar_kwarg, typ in zip(
                ds_2d, axes, plot_cbars, cbar_axes, plot_yaxes, plot_xaxes, cbar_kwargs, 
                range(1,6)
                ):
            ds.plot.contourf(levels=levels, colors=colors, ax=ax, 
                             add_colorbar = plot_cbar, cbar_ax=cbar_ax, cbar_kwargs=cbar_kwarg)
            
            ax.xaxis.label._color = light
            ax.yaxis.label._color = light
        
            if plot_yaxis == False:
                ax.set_ylabel("")
            else:
                ax.set_ylabel("z [10 m]")
            if plot_xaxis == False:
                ax.set_xlabel("")
                ax.set_xticklabels("")
            else:
                ax.set_xlabel("x [100 km]")
            
            for spine in ax.spines.values():
                spine.set_edgecolor(dark)

            ax.text(0.1, 0.5, str(typ), horizontalalignment='center',
                 verticalalignment='center', transform=ax.transAxes, color=light, 
                 bbox=dict(edgecolor=light, facecolor=dark))

        cbar_axin.set_ylabel("salinity (TDS g/l)", rotation=270, labelpad=10.)
        cbar_axin.xaxis.label._color = light
        cbar_axin.yaxis.label._color = light
        plt.savefig(outf + ".png", facecolor=fig.get_facecolor(), edgecolor='none', dpi=300)
        plt.savefig(outf + ".svg", facecolor=fig.get_facecolor(), edgecolor='none')