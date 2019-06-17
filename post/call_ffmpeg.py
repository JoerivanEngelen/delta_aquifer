# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:54:43 2019

@author: engelen
"""

import subprocess
import os
from glob import glob


#%%Path management
ffmpeg = r"c:\Users\engelen\OneDrive - Stichting Deltares\ffmpeg\ffmpeg-20160412-git-196cfc2-win64-static\bin\ffmpeg.exe"

#video_fol = r"g:\synthdelta\test_output\SD_i202\synth_SD_i202_m24_6378797\screenshots"
video_fol = r"g:\synthdelta\test_output\*\*\screenshots"
models = glob(video_fol)

#%%
#args = [ffmpeg, "-framerate", "5", "-i", None, "-q", "1", "-r", "20", None]
args = [ffmpeg, "-framerate", "5", "-i", None, "-q", "1", "-r", "25", None, "-y"]

inp, outp = [i for i, arg in enumerate(args) if arg is None]
#%%
for modfol in models:
    args[inp] = os.path.join(modfol, "results" + ".%04d.png")
    args[outp]= os.path.join(modfol, "results" + ".avi")
    proc = subprocess.run(args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print(proc.stderr)