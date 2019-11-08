# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:19:07 2019

@author: engelen

Script to convert the trajectories to a textfile which can be used for the gifmatrix

"""

from pkg_resources import resource_filename
import pandas as pd
import os
import numpy as np

datafol= os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))
traj_id = os.path.join(datafol, "traj_id.csv")

traj_len=24

df = pd.read_csv(traj_id, index_col=0)
diff = df.rolling(window=2).apply(lambda x: x[1] - x[0], raw=True)
diff.loc[0::traj_len] = np.nan
base_id = list(diff.loc[0::traj_len].index)

changed_parameters = diff[diff != 0].stack().index.tolist()
changed_parameters.extend([(idx, "base") for idx in base_id])
changed_parameters.sort()

names=list(zip(*changed_parameters))[1]