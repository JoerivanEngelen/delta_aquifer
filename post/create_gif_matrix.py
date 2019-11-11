# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:21:38 2019

@author: engelen
"""

from glob import glob
from PIL import Image, ImageSequence, ImageDraw, ImageFont
from pkg_resources import resource_filename
from matplotlib import font_manager
import sys, os
import pandas as pd

#%%TODO
#Time + (Colorbar?)

#%%Path management
if len(sys.argv) > 1:
    globpath  = sys.argv[1]
    out_fol = sys.argv[2]
    start = sys.argv[3]
    stop = sys.argv[4]
else:
    #Local testing on my own windows laptop
    globpath = r"g:\synthdelta\results\gifs\*.gif"
    out_fol = r"g:\synthdelta\results"
    start = 0
    stop = 23

mod_idx = slice(start, stop+1)

files = glob(globpath)
files.sort()

#Get default font path of Matplotlib
font_path = font_manager.findfont(None)

#Path with text aid
text_path = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "text_aid.csv")))

out_gif_path = os.path.join(out_fol, "SD_i%s-%s.gif" % (start, stop))

#%%Specify output image
nrows, ncols = 3, 8

#%%Get texts
texts = list(pd.read_csv(text_path, index_col=0)["par"][mod_idx])

#%%Load and peek
ims = [Image.open(f) for f in files[mod_idx]]
#w_original, h_original = ims[0].size
#Original image 1000, 850
cropbox = (220, 40, 850, 850)
w_crop_original = cropbox[2]-cropbox[0]
h_crop_original = cropbox[3]-cropbox[1]

#%%Prepare
font = ImageFont.truetype(font_path, 16)

new_size = int(w_crop_original/ncols)*2, int(h_crop_original/nrows)
output_size = new_size[0] * ncols, new_size[1] * nrows
ids = [(int(i/ncols), i%ncols) for i in range(len(ims))]

out_frames = []
iterators=[ImageSequence.Iterator(im) for im in ims]

#%%Process
for i, frames in enumerate(zip(*iterators)):
    #We convert everything to RGB, as each .gif has a different color lookup table
    #In the end convert to Palette mode
    stitched_image = Image.new("RGB", output_size)
    draw = ImageDraw.Draw(stitched_image)
    resized_frames = [frame.crop(cropbox).resize(new_size).convert("RGB") for frame in frames]
    for i, resized_frame in enumerate(resized_frames):
        r, c=ids[i]
        sub_extent=(c*new_size[0], r*new_size[1], (c+1)*new_size[0], (r+1)*new_size[1])
        stitched_image.paste(resized_frame, sub_extent)
        draw.text((c*new_size[0], r*new_size[1]),
                  texts[i],(255,255,255), font=font)
    del draw
    stitched_image = stitched_image.convert("P", palette=Image.ADAPTIVE)
    out_frames.append(stitched_image)

#%%Save
out_frames[0].save(out_gif_path, format='GIF', 
          append_images=out_frames[1:], save_all=True, duration=100, loop=0,
          include_color_table=True)