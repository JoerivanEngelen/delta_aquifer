# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 14:59:38 2019

@author: engelen
"""

import numpy as np
import pandas as pd
import seaborn.utils as utils
from seaborn.utils import remove_na
from seaborn.categorical import _PointPlotter
from seaborn.algorithms import bootstrap
import matplotlib.pyplot as plt

"""Override the following function to allow for error bars to present the data range.
"""


def pointplotrange(
    x=None,
    y=None,
    hue=None,
    data=None,
    order=None,
    hue_order=None,
    estimator=np.mean,
    ci=95,
    n_boot=1000,
    units=None,
    seed=None,
    markers="o",
    linestyles="-",
    dodge=False,
    join=True,
    scale=1,
    orient=None,
    color=None,
    palette=None,
    errwidth=None,
    capsize=None,
    ax=None,
    **kwargs
):

    plotter = _PointPlotterRange(
        x,
        y,
        hue,
        data,
        order,
        hue_order,
        estimator,
        ci,
        n_boot,
        units,
        seed,
        markers,
        linestyles,
        dodge,
        join,
        scale,
        orient,
        color,
        palette,
        errwidth,
        capsize,
    )

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax)
    return ax


class _PointPlotterRange(_PointPlotter):
    def estimate_statistic(self, estimator, ci, n_boot, seed):

        if self.hue_names is None:
            statistic = []
            confint = []
        else:
            statistic = [[] for _ in self.plot_data]
            confint = [[] for _ in self.plot_data]

        for i, group_data in enumerate(self.plot_data):

            # Option 1: we have a single layer of grouping
            # --------------------------------------------

            if self.plot_hues is None:

                if self.plot_units is None:
                    stat_data = remove_na(group_data)
                    unit_data = None
                else:
                    unit_data = self.plot_units[i]
                    have = pd.notnull(np.c_[group_data, unit_data]).all(axis=1)
                    stat_data = group_data[have]
                    unit_data = unit_data[have]

                # Estimate a statistic from the vector of data
                if not stat_data.size:
                    statistic.append(np.nan)
                else:
                    statistic.append(estimator(stat_data))

                # Get a confidence interval for this estimate
                if ci is not None:

                    if stat_data.size < 2:
                        confint.append([np.nan, np.nan])
                        continue

                    if ci == "sd":

                        estimate = estimator(stat_data)
                        sd = np.std(stat_data)
                        confint.append((estimate - sd, estimate + sd))

                    elif ci == "range":
                        confint.append((np.min(stat_data), np.max(stat_data)))

                    else:

                        boots = bootstrap(
                            stat_data, func=estimator, n_boot=n_boot, units=unit_data
                        )
                        confint.append(utils.ci(boots, ci))

            # Option 2: we are grouping by a hue layer
            # ----------------------------------------

            else:
                for j, hue_level in enumerate(self.hue_names):

                    if not self.plot_hues[i].size:
                        statistic[i].append(np.nan)
                        if ci is not None:
                            confint[i].append((np.nan, np.nan))
                        continue

                    hue_mask = self.plot_hues[i] == hue_level
                    if self.plot_units is None:
                        stat_data = remove_na(group_data[hue_mask])
                        unit_data = None
                    else:
                        group_units = self.plot_units[i]
                        have = pd.notnull(np.c_[group_data, group_units]).all(axis=1)
                        stat_data = group_data[hue_mask & have]
                        unit_data = group_units[hue_mask & have]

                    # Estimate a statistic from the vector of data
                    if not stat_data.size:
                        statistic[i].append(np.nan)
                    else:
                        statistic[i].append(estimator(stat_data))

                    # Get a confidence interval for this estimate
                    if ci is not None:

                        if stat_data.size < 2:
                            confint[i].append([np.nan, np.nan])
                            continue

                        if ci == "sd":

                            estimate = estimator(stat_data)
                            sd = np.std(stat_data)
                            confint[i].append((estimate - sd, estimate + sd))

                        elif ci == "range":
                            confint[i].append((np.min(stat_data), np.max(stat_data)))

                        else:

                            boots = bootstrap(
                                stat_data,
                                func=estimator,
                                n_boot=n_boot,
                                units=unit_data,
                                seed=seed,
                            )
                            confint[i].append(utils.ci(boots, ci))

        # Save the resulting values for plotting
        self.statistic = np.array(statistic)
        self.confint = np.array(confint)
