# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:21:38 2019

@author: engelen
"""

from glob import glob
from PIL import Image, ImageSequence
import sys, os

#%%TODO
#Crop images
#Hoe kan ik een color lookup table geven?

#%%Path management
if len(sys.argv) > 1:
    globpath  = sys.argv[1]
    out_fol = sys.argv[2]
    start = sys.argv[3]
    stop = sys.argv[4]
    mod_idx = slice(int(start), int(stop))
else:
    #Local testing on my own windows laptop
    globpath = r"g:\synthdelta\results\gifs\*.gif"
    out_fol = r"g:\synthdelta\results"
    mod_idx = slice(0, 24)

files = glob(globpath)
files.sort()

#%%Specify output image
nrows, ncols = 4, 6

#%%Load and peek
ims = [Image.open(f) for f in files[mod_idx]]
w_original, h_original = ims[0].size

#%%Prepare
new_size = int(w_original/ncols), int(h_original/nrows)
output_size = new_size[0] * ncols, new_size[1] * nrows
ids = [(int(i/ncols), i%ncols) for i in range(len(ims))]

out_frames = []
iterators=[ImageSequence.Iterator(im) for im in ims]

#%%Process
for frames in zip(*iterators):
    stitched_image = Image.new("RGB", output_size)
    resized_frames = [frame.resize(new_size) for frame in frames]
    for i, resized_frame in enumerate(resized_frames):
        r, c=ids[i]
        sub_extent=(c*new_size[0], r*new_size[1], (c+1)*new_size[0], (r+1)*new_size[1])
        stitched_image.paste(resized_frame, sub_extent)
    out_frames.append(stitched_image)

#%%Save
out_frames[0].save(os.path.join(out_fol, "total.gif"), format='GIF', 
          append_images=out_frames[1:], save_all=True, duration=100, loop=0,
          include_color_table=True)