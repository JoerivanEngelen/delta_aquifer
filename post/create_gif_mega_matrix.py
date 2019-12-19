# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 13:20:38 2019

Script to convert multiple multi-panel .gifs to one big one.

@author: engelen
"""

from PIL import Image, ImageSequence
import os
from glob import glob
import numpy as np


#%%Path management
out_fol = r"g:\synthdelta\results"
out_file = os.path.join(out_fol, "mega_matrix.gif")

files=glob(os.path.join(out_fol, "SD*.gif"))
files.sort()

#%%Process

nrows, ncols = 4, 3

dcell = 10

ims = [Image.open(f) for f in files]
n_frames = [im.n_frames for im in ims]
#max_n = max(n_frames) #For some reason there are some gifs with 47 frames
max_n = 46
iterators=[ImageSequence.Iterator(im) for im in ims]

mag = 1
#mag = 6

size = ims[0]._size

panel_size = int(size[0]/ncols)*mag, int(size[0]/nrows)*mag
output_size = panel_size[0] * ncols + (ncols-1) * dcell, panel_size[1] * nrows + (nrows-1) * dcell

idxs = [(int(i/ncols), i%ncols) for i in range(len(ims))]

out_frames = []
for i in range(max_n):
    stitched_image = Image.new("RGBA", output_size)

    frames = [iterator[i] for j, iterator in enumerate(iterators)]
    resized_frames = [frame.resize(panel_size).convert("RGBA") for frame in frames]
    
    for j, resized_frame in enumerate(resized_frames):
        r, c = idxs[j]
        sub_extent=np.array([c*(panel_size[0]), r*(panel_size[1]), (c+1)*panel_size[0], (r+1)*panel_size[1]])
        sub_extent[::2]  += dcell * c
        sub_extent[1::2] += dcell * r
        stitched_image.paste(resized_frame, sub_extent)

    stitched_image = stitched_image.convert("P", palette=Image.ADAPTIVE)
    out_frames.append(stitched_image)

out_frames[0].save(out_file, format='GIF', 
          append_images=out_frames[1:], save_all=True, duration=130, loop=0,
          include_color_table=True)