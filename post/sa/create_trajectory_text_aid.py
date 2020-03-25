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


#%%Path management
datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))
traj_id  = os.path.join(datafol, "traj_id.csv")
out_path = os.path.join(datafol, "text_aid.csv")

#%%Process
traj_len=24

df = pd.read_csv(traj_id, index_col=0)
diff = df.rolling(window=2).apply(lambda x: x[1] - x[0], raw=True)
diff.loc[0::traj_len] = np.nan

base_id = list(diff.loc[0::traj_len].index)
sign = np.sum(np.sign(diff), axis=1) #Direction parameters changed
convert = {1.0 : "+",
           -1.0: "-",
           0.0:  ""}

changed_parameters = diff[diff != 0].stack().index.tolist()
changed_parameters.extend([(idx, "base") for idx in base_id])
changed_parameters.sort()

changed_parameters = pd.DataFrame(data=changed_parameters).set_index(0)
changed_parameters.columns = ["par"]
changed_parameters["par"] += sign.map(convert)

changed_parameters.to_csv(out_path)