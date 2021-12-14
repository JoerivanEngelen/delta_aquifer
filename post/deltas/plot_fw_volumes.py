# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 09:43:24 2020

@author: engelen
"""

import xarray as xr
import cftime
import pandas as pd
import os
from pkg_resources import resource_filename
from glob import glob
import numpy as np
import re

from pathlib import Path

from analyze_input.ci_range import pointplotrange

from itertools import product

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#%%Functions to parse
def get_model_id(file):
    pattern = r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return number


#%%Function to process data
def time_in_ka(ds):
    return [t.year / 1000 for t in ds.time.values[::-1]]


def read_and_process_ds(files, delta_len, append, start, stop):
    files_mod = [i for i in files if get_model_id(i) in range(start, stop)]
    #    assert(len(files_mod)==delta_len)

    ds_ls = [xr.open_dataset(f).sortby("time", ascending=True) for f in files_mod]
    ds_ls = [ds.assign_coords(time=time_in_ka(ds)) for ds in ds_ls]
    ds_ls = [append_varname(ds, append=append) for ds in ds_ls]
    return ds_ls


def append_varname(ds, append):
    keys = ds.keys()
    renamed = dict([(key, key + append) for key in keys])
    return ds.rename(**renamed)


def calc_ddim(da, dim):
    shifts = {dim: 1}
    ddim = da[dim] - da[dim].shift(**shifts)
    ddim[0] = ddim[1]
    return ddim


#%%Function to adapt plot
def label_offset(ax, axis="y", label=None, unit=None):
    if axis == "y":
        fmt = ax.yaxis.get_major_formatter()
        ax.yaxis.offsetText.set_visible(False)
        set_label = ax.set_ylabel
    #        label = ax.get_ylabel()

    elif axis == "x":
        fmt = ax.xaxis.get_major_formatter()
        ax.xaxis.offsetText.set_visible(False)
        set_label = ax.set_xlabel
        label = ax.get_xlabel()

    def update_label(event_axes):
        offset = fmt.get_offset()
        if offset == "":
            set_label("{} ({})".format(label, unit))
        else:
            minus_sign = u"\u2212"
            expo = np.float(offset.replace(minus_sign, "-").split("e")[-1])
            set_label("%s ($\mathregular{10^{%d}}$ %s)" % (label, expo, unit))

    ax.callbacks.connect("ylim_changed", update_label)
    ax.callbacks.connect("xlim_changed", update_label)
    ax.figure.canvas.draw()
    update_label(None)


def logmean(arr):
    return 10 ** np.mean(np.log10(arr + 1))


#%%Path management
datafol = os.path.abspath(
    resource_filename("delta_aquifer", os.path.join("..", "data", "30_deltas"))
)
val_fol = os.path.join(datafol, "..", "validation")
res_fol = (
    r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results_30deltas"
)

mas_f = glob(os.path.join(res_fol, "mas", "mas_i*.nc"))
mas_pump_f = glob(os.path.join(res_fol, "mas", "mas_pump_i*.nc"))
vol_f = glob(os.path.join(res_fol, "vols", "vol_i*.nc"))
vol_pump_f = glob(os.path.join(res_fol, "vols", "vol_pump_i*.nc"))
val_f = glob(os.path.join(val_fol, "FW_vol*.csv"))

recharge_path = os.path.join(res_fol, "recharges", "recharges.csv")

par = pd.read_csv(os.path.join(datafol, r"model_inputs.csv"), index_col=0)
abstractions = glob(os.path.join(datafol, "abstractions", "*.nc"))

#%%Create start and endpoints
delta_len = 9

# to_finish = np.array([63, 198, 207])

starts = np.array([0, 9, 18, 27, 36, 45, 54, 63, 81, 90, 153, 171, 198, 207, 216])
# starts = starts[~np.isin(starts, to_finish)]
stops = starts + delta_len

k_values = list(product(["min ", "max ", "mean "], ["max", "min", "mean"]))
k_values = pd.DataFrame(
    [[i, j] for (i, j) in list(k_values)], columns=["Kh_aqf", "Kv_aqt"]
)

#%%Process to long format for seaborn
not_porous = ["dx", "dy", "Kh_aqf", "Kv_aqt", "time", "x_sal", "x_sal_pump", "delta"]

ds_deltas = {}
for start, stop in zip(starts, stops):
    delta = par.loc[start, "Delta"]
    por = par.loc[start, "n"]

    file_append = zip(
        [mas_f, vol_f, mas_pump_f, vol_pump_f], ["", "", "_pump", "_pump"]
    )
    ds_ls = [
        read_and_process_ds(f, delta_len, append, start, stop)
        for (f, append) in file_append
    ]

    ds_deltas[delta] = [xr.merge(ds_ls_sim) for ds_ls_sim in zip(*ds_ls)]

    ds_deltas[delta] = [
        ds.to_dataframe().reset_index().assign(**k_values.loc[i].to_dict())
        for i, ds in enumerate(ds_deltas[delta])
    ]

    ds_deltas[delta] = pd.concat(ds_deltas[delta])
    porous = [var for var in ds_deltas[delta].columns if var not in not_porous]
    ds_deltas[delta].loc[:, porous] = ds_deltas[delta].loc[:, porous] * por

df_deltas = pd.concat([ds.assign(delta=delta) for delta, ds in ds_deltas.items()])
df_deltas = df_deltas.rename(
    columns={"time": "time (ka)", "Kh_aqf": "$K_{h,aqf}$", "Kv_aqt": "$K_{v,aqt}$"}
)
#%%Validations
val_df = dict(
    [
        (
            os.path.splitext(os.path.basename(path))[0].split("FW_vol_")[1],
            # Reindex to pad with NaNs
            pd.read_csv(path).reindex(range(delta_len)),
        )
        for path in val_f
    ]
)

# Add Mekong data
# Gunnink, J.L., Hung, V.P., Oude Essink, G.H.P., Bierkens, M.F.P. 2020.
# Geostatistical estimation of the 3D groundwater salinity distribution in the
# Mekong Delta, Vietnam. Submitted to HESSD

val_df["Mekong"] = pd.DataFrame(data={"FW_vols": [830e9, 876e9, 890e9]}).reindex(
    range(delta_len)
)

#%% Estimations
est_df = pd.read_csv(Path(val_fol) / "FW_estimations.csv", index_col=0)

#%%Plot settings
plot_df = False

agu_whole = (19 / 2.54, 23 / 2.54)
col_wrap = 4

height = agu_whole[1] / col_wrap
aspect = agu_whole[0] / agu_whole[1]

sns.set(style="whitegrid")

plt_kwargs = dict(
    x="time (ka)",
    hue="$K_{h,aqf}$",
    style="$K_{v,aqt}$",
    col="delta",
    col_wrap=col_wrap,
    data=df_deltas,
    kind="line",
    hue_order=["max ", "mean ", "min "],
    style_order=["max", "mean", "min"],
    height=height,
    aspect=aspect,
    facet_kws=dict(sharey=False, legend_out=False, margin_titles=False, xlim=(125, 1)),
)

#%%Plot
to_plot = [
    "fw_onshore",
    "fw_onshore_pump",
    "fw_offshore",
    "fw_offshore_pump",
    "sal_onshore",
    "sal_onshore_pump",
    "ol_sal_onshore",
    "ol_sal_onshore_pump",
    "x_sal",
    "x_sal_pump",
]

palettes = [
    "Blues_r",
    "Blues_r",
    "Blues_r",
    "Blues_r",
    "Reds_r",
    "Reds_r",
    "Greys_r",
    "Greys_r",
    "Oranges_r",
    "Oranges_r",
]

ylabels = [
    "$V_{fw,on}$",
    "$V_{fw,on}$",
    "$V_{fw,off}$",
    "$V_{fw,off}$",
    "$S_{on}$",
    "$S_{on}$",
    "$S_{init}$",
    "$S_{init}$",
    "$x_{toe}$",
    "$x_{toe}$",
]

units = ["m3", "m3", "m3", "m3", "kg", "kg", "kg", "kg", "m", "m"]

if plot_df == True:
    for var, palette, ylabel, unit in zip(to_plot, palettes, ylabels, units):
        g = sns.relplot(y=var, palette=palette, **plt_kwargs)
        for ax in g.axes:
            ax.set_ylim(ymin=0)
            label_offset(ax, axis="y", label=ylabel, unit=unit)
            ax.set_xscale("log")
            ax.set_xticks([100, 10, 1])

        g.axes[0].legend(
            loc="lower right",
            bbox_to_anchor=(0.975, 0.025),
            bbox_transform=g.fig.transFigure,
        )

        plt.subplots_adjust(hspace=0.25, wspace=0.8)
        g.set_titles(template="{col_name}")
        g.savefig(os.path.join(res_fol, "{}.png".format(var)), dpi=300)
        g.savefig(os.path.join(res_fol, "{}.pdf".format(var)), dpi=300)
        plt.close()

#%%Process abstractions
deltas = pd.unique(df_deltas["delta"])
name = "__xarray_dataarray_variable__"

abs_ls = [
    xr.open_dataset(os.path.join(datafol, "abstractions", "{}.nc".format(delta)))[name]
    for delta in deltas
]
cellsizes = [calc_ddim(da, "x") * calc_ddim(da, "y") for da in abs_ls]
abs_ls = [da * cellsizes[i] for i, da in enumerate(abs_ls)]

end_abs = [da.isel(time=-1).sum().values for da in abs_ls]
end_abs = pd.DataFrame(data=end_abs, index=deltas, columns=["Q"])
end_abs = end_abs.reset_index().rename(columns={"index": "delta"})

##Update deltas for which we have data
# end_abs.loc[end_abs["delta"]=="Nile",   "Q"] = 7.162345e6 * 365.25
# end_abs.loc[end_abs["delta"]=="Mekong", "Q"] = 2.4e6 * 365.25

#%%Process recharges
rch = pd.read_csv(recharge_path, index_col=0)
rch = rch.sort_index()
rch = rch.join(par["Delta"]).rename(columns={"Delta": "delta"})
# Recharge constant through time for paleo reconstruction,
# so dropping duplicates results in the present-day recharge.
rch = rch.drop_duplicates()

#%%Process historical abstractions
hist_abs = [da.sel(time=slice(1958, 2013)).sum(dim=["x", "y"]) for da in abs_ls]
hist_abs = [
    da.assign_coords(delta=delta).to_dataframe(name="Q_hist")
    for da, delta in zip(hist_abs, deltas)
]
hist_abs = pd.concat(hist_abs).set_index(["delta"], append=True)

# Compute historical depletion
hist_depletion = pd.DataFrame()
hist_depletion["Q_hist"] = hist_abs["Q_hist"] - rch.set_index("delta")["rch"]
hist_depletion["Q_hist_4"] = hist_abs["Q_hist"] * 4 - rch.set_index("delta")["rch"]
hist_depletion["Q_hist_0.4"] = hist_abs["Q_hist"] * 0.41 - rch.set_index("delta")["rch"]

# Assume all excess recharge water is drained,
# so we can clip negative values
hist_depletion = hist_depletion.clip(lower=0.0)

# Sum across time
hist_depletion = hist_depletion.groupby("delta").sum()

#%%Create formatted table
col_select = ["fw_onshore_pump", "fw_offshore", "ol_sal_onshore"]
to_pointplot = col_select + ["t_depleted"]
to_log = ["fw_onshore_pump", "fw_offshore", "t_depleted"]
xlabels = [
    "$V_{fw, on}$ (m3)",
    "$V_{fw, off}$ (m3)",
    "$S_{init}$ (%)",
    "$t_{d}$ (year)",
]
colors = ["navy", "royalblue", "firebrick", "darkslategray"]

# Process
df_fin = df_deltas.loc[
    df_deltas["time (ka)"] == 1.0, (col_select + ["delta", "sal_onshore"])
]
df_fin = df_fin.merge(end_abs, on="delta")
df_fin = df_fin.merge(rch, on="delta")
df_fin = df_fin.merge(hist_depletion, on="delta")

df_fin["t_depleted"] = (df_fin["fw_onshore_pump"] - df_fin["Q_hist"]).clip(
    lower=0.0
) / (df_fin["Q"] - df_fin["rch"])
df_fin["t_depleted_4"] = (df_fin["fw_onshore_pump"] - df_fin["Q_hist_4"]).clip(
    lower=0.0
) / ((df_fin["Q"] * 4) - df_fin["rch"])
df_fin["t_depleted_0.4"] = (df_fin["fw_onshore_pump"] - df_fin["Q_hist_0.4"]).clip(
    lower=0.0
) / ((df_fin["Q"] * 0.41) - df_fin["rch"])

df_fin["ol_sal_onshore"] = df_fin["ol_sal_onshore"] / df_fin["sal_onshore"]

# df_fin = df_fin.loc[df_fin["fw_onshore_pump"]!=0] #Saloum has a few simulations where V_fw = 0., these cannot be True.
df_fin = df_fin.sort_values("delta")

# Add observed FW vols
obs_name = "fw_onshore_observed"
df_fin[obs_name] = np.nan
for delta in val_df.keys():
    df_fin.loc[df_fin["delta"] == delta, obs_name] = val_df[delta]["FW_vols"].values

# Add estimated FW vols
est_name = "fw_onshore_estimated"
est_df = est_df.set_index(est_df.index.str.replace("Red_River", "Red River"))

df_fin[est_name] = np.nan
for delta in est_df.index.values:
    df_fin.loc[df_fin["delta"] == delta, est_name] = est_df.loc[delta, "V_fw,on"]


depletion_columns = ["t_depleted_0.4", "t_depleted", "t_depleted_4"]
# Negative depletion times are infinite
for colname in depletion_columns:
    df_fin.loc[df_fin[colname] < 0, [colname]] = 1e8

# Set everything between 0 and 1 to
for colname in depletion_columns:
    df_fin.loc[df_fin[colname] < 1, [colname]] = 1.0


#%% Plot
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=agu_whole)
axes = axes.flatten()

opts = dict(ci="range", join=False, errwidth=1.5, scale=0.4)
for i, var in enumerate(to_pointplot):
    if var in to_log:
        estimator = logmean
    else:
        estimator = np.mean

    if var == "t_depleted":
        pointplotrange(
            y="delta",
            x="t_depleted_4",
            data=df_fin,
            estimator=estimator,
            ax=axes[i],
            color="peru",
            **opts
        )
        pointplotrange(
            y="delta",
            x="t_depleted_0.4",
            data=df_fin,
            estimator=estimator,
            ax=axes[i],
            color="burlywood",
            **opts
        )
        axes[i].axvline(200, linestyle=":", color="darkgray")
        customlines = [
            Line2D([0], [0], color=c) for c in ["burlywood", colors[i], "peru"]
        ]
        axes[i].legend(
            customlines,
            ["$0.4 \; Q$", "$1.0 \; Q$", "$4.0 \; Q$"],
            loc="upper right",
            bbox_to_anchor=[0.94, 1.01],
        )

    pointplotrange(
        y="delta",
        x=var,
        data=df_fin,
        estimator=estimator,
        ax=axes[i],
        color=colors[i],
        **opts
    )

    if var == "fw_onshore_pump":
        col_obs = "limegreen"
        col_est = "tomato"
        m_obs = "d"
        m_est = "^"

        pointplotrange(
            y="delta",
            x="fw_onshore_observed",
            data=df_fin,
            estimator=estimator,
            ax=axes[i],
            color=col_obs,
            markers=m_obs,
            **opts
        )
        pointplotrange(
            y="delta",
            x="fw_onshore_estimated",
            data=df_fin,
            estimator=estimator,
            ax=axes[i],
            color=col_est,
            markers=m_est,
            **opts
        )

        colors_vol = [colors[i], col_obs, col_est]
        markers = ["o", m_obs, m_est]

        customlines = [
            Line2D([0], [0], color=c, marker=m, markerfacecolor=c)
            for c, m in zip(colors_vol, markers)
        ]
        axes[i].legend(
            customlines,
            ["simulated", "observed", "estimated"],
            loc="upper left",
            bbox_to_anchor=[-0.01, 1.01],
        )

    labels = dict(xlabel=xlabels[i])
    if (i % 2) == 1:
        labels["ylabel"] = ""
    axes[i].set(**labels)
    if var in to_log:
        axes[i].set_xscale("log")

    if var in ["fw_onshore_pump", "fw_offshore"]:
        axes[i].set_xticks([1e4, 1e6, 1e8, 1e10, 1e12])
    elif var == "t_depleted":
        axes[i].set_xticks([1e0, 1e2, 1e4, 1e6, 1e8])
        axes[i].set_xticklabels(
            [r"$\mathregular{10^{%d}}$" % i for i in [0, 2, 4, 6]] + [r"$\infty$"]
        )

plt.tight_layout()
plt.savefig(os.path.join(res_fol, "end-states.png"), dpi=300)
plt.savefig(os.path.join(res_fol, "end-states.pdf"))
plt.close()

# %%
