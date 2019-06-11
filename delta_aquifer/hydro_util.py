# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:46:42 2019

@author: engelen
"""

def rayleigh(rho_f, rho_s, D, diff, Kv):
    return((rho_s-rho_f)*Kv*D/(rho_f * diff))
