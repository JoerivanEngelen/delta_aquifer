# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:21:38 2019

@author: engelen
"""

from glob import glob
from PIL import Image, ImageSequence
from pkg_resources import resource_filename
import matplotlib.pyplot as plt
import sys, os
import pandas as pd
import re
import numpy as np

#%%Functions to parse 
def get_model_id(file):     
    pattern=r"i([0-9][0-9][0-9])"
    number = int(re.search(pattern, file).group(1))
    return(number)

def str_to_mathtxt(s):
    s, sign = s[:-1], s[-1]
    s_ls = s.split("_")
    var = s_ls[0]
    #Create subscript
    if len(s_ls) > 1:
        sub = s_ls[1:]
        sub = ",".join(sub)
        sub = "_{%s}" % sub
    else:
        sub=""
    s = "${}{}{}$".format(var, sub, sign)
    return(s)

#%%Functions to convert matplotlib figs to pillow Images
def _get_text_imagearray(text, size):
    dpi=300
    fontsize=5
    #dpi=1
    #fontsize=10000
    
    figsize=tuple([s/dpi for s in size])
    
    fig = plt.figure(figsize=figsize, dpi=dpi)
    plt.text(0.0, 1.0, text, transform=fig.transFigure, fontsize=fontsize, ha="left", va="top")
    plt.axis('off')
    fig.canvas.draw()
    s, (width, height) = fig.canvas.print_to_buffer()
    arr = np.frombuffer(s, dtype='uint8').reshape(height, width, 4).copy()
    plt.close()
    return(arr)

def _black2white(arr, background=110):
    """Convert black to white and set black transparent. 
    Only works for grayscale images.
    
    Parameters
    ----------
    arr : np.ndarray(N, M, 4)
        RGBA array with grayscale array
    
    """
    arr[:, :, 3] = 255 - arr[:, :, 0]  
    arr[:, :, :3] = (255-background) * (arr[:, :, :3]/255)
    arr[:, :, :3] = 255 - arr[:, :, :3]
    return(arr)

def get_text_image(text, size):
    arr = _get_text_imagearray(text, size)
    arr = _black2white(arr)
    image = Image.fromarray(arr)
    return(image)

def get_texts(text_path, start, stop):
    tex = {"alpha" : r"\alpha",
           "beta"  : r"\beta", 
           "phi_f"   : r"\phi_f",
           "Kh_Kv" : r"Kh/Kv"}
    
    mod_idx = slice(start, stop+1)
    texts = pd.read_csv(text_path, index_col=0)["par"][mod_idx]
    texts = list(texts)
    texts = [tex[t[:-1]]+t[-1] if t[:-1] in tex.keys() else t for t in texts]
    texts = [str_to_mathtxt(s) for s in texts]
    return(texts)

#%%TODO
#Time + (Colorbar?)

#%%Path management
if len(sys.argv) > 1:
    globpath  = sys.argv[1]
    out_fol = sys.argv[2]
else:
    #Local testing on my own windows laptop
    globpath = r"g:\synthdelta\results\gifs\*.gif"
    out_fol = r"g:\synthdelta\results"

files = glob(globpath)
files.sort()

#Path with text aid
text_path = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data", "text_aid.csv")))

#%%
starts = np.array([0, 24, 48, 72, 120, 144, 168, 192, 240, 264])
stops = starts+23

for start, stop in zip(starts, stops):
    files_mod = [i for i in files if get_model_id(i) in range(start, stop+1)]
    
    out_gif_path = os.path.join(out_fol, "SD_i%s-%s.gif" % (start, stop))

    #%%Specify output image size
    cropbox = (220, 40, 850, 850) #Original image 1000, 850
    w_crop_original = cropbox[2]-cropbox[0]
    h_crop_original = cropbox[3]-cropbox[1]
    
    nrows, ncols = 3, 8 

    panel_size = int(w_crop_original/ncols)*2, int(h_crop_original/nrows)
    output_size = panel_size[0] * ncols, panel_size[1] * nrows
    
    #%%Load
    ims = [Image.open(f) for f in files_mod]
    n_frames = [im.n_frames for im in ims]
    assert(len(files_mod)==23)

    #%%Prepare
    texts = get_texts(text_path, start, stop) #Get texts
    
    idxs = [(int(i/ncols), i%ncols) for i in range(len(ims))]
    
    max_n = max(n_frames)
    start_frames = [max_n-n_frame for n_frame in n_frames]
    iterators=[ImageSequence.Iterator(im) for im in ims]
    
    #%%Create text overlay
    text_overlay = Image.new("RGBA", output_size)
    
    for j, (r, c) in enumerate(idxs):
        sub_extent=(c*panel_size[0], r*panel_size[1], (c+1)*panel_size[0], (r+1)*panel_size[1])
        text_im = get_text_image(texts[j], panel_size)
        text_overlay.paste(text_im, sub_extent)
    
    #%%Process
    out_frames = []
    for i in range(max_n):
        #We convert everything to RGB, as each .gif has a different color lookup table
        #In the end convert to Palette mode
        stitched_image = Image.new("RGBA", output_size)
        #Select frame if available, else show first frame.
        frames = [iterator[i-start_frames[j]] if i >= start_frames[j] else iterator[0] for j, iterator in enumerate(iterators)]
        
        resized_frames = [frame.crop(cropbox).resize(panel_size).convert("RGBA") for frame in frames]
        for j, resized_frame in enumerate(resized_frames):
            r, c = idxs[j]
            sub_extent=(c*panel_size[0], r*panel_size[1], (c+1)*panel_size[0], (r+1)*panel_size[1])
            stitched_image.paste(resized_frame, sub_extent)
        
        stitched_image = Image.alpha_composite(stitched_image, text_overlay)
        stitched_image = stitched_image.convert("P", palette=Image.ADAPTIVE)
        out_frames.append(stitched_image)
    
    #%%Save
    out_frames[0].save(out_gif_path, format='GIF', 
              append_images=out_frames[1:], save_all=True, duration=100, loop=0,
              include_color_table=True)