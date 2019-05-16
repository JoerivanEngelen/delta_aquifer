# -*- coding: utf-8 -*-
"""
Calculate modflow timestepping with tsmult

By: Joeri van Engelen

"""

import numpy as np

def dt1(perlen, tsmult, nstp):
    return(perlen*(tsmult-1)/(tsmult**nstp-1))
    
def get_dt(perlen, tsmult, nstp):
    dt=[]
    dt.append(dt1(perlen, tsmult, nstp))
    for i in range(0, nstp-1):
        dt.append(dt[i]*tsmult)
    return(np.array(dt))


dt = get_dt(1000.*365.25, 5., 8)


print(dt)
print(np.cumsum(dt))