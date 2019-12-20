# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:21:38 2019

@author: engelen
"""

from glob import glob
from PIL import Image, ImageSequence, ImageEnhance
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

#%%Functions to create texts
def get_changes(traj_id, traj_len):
    df = pd.read_csv(traj_id, index_col=0)
    diff = df.rolling(window=2).apply(lambda x: x[1] - x[0], raw=True)
    diff.loc[0::traj_len] = np.nan
    return(diff)

def create_text_aid(diff, traj_len, text_path=None):
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
    
    if text_path is not None:
        changed_parameters.to_csv(text_path)
    
    return(changed_parameters)

def create_texts(diff, traj_len):
    tex = {"alpha" : r"\alpha",
           "beta"  : r"\beta", 
           "phi_f"   : r"\phi_f",
           "Kh_Kv" : r"Kh/Kv"}
    
    text_aid = create_text_aid(diff, traj_len)
    texts = text_aid["par"]
    texts = list(texts)
    texts = [tex[t[:-1]]+t[-1] if t[:-1] in tex.keys() else t for t in texts]
    texts = [str_to_mathtxt(s) for s in texts]
    return(texts)

#%%Functions to convert matplotlib figs to pillow Images
def _get_text_imagearray(text, size):
    dpi=300
    
    lp = (270/300) #Appropriate scaling relationship for fontsize
    mag = size[1]/dpi /lp
    fontsize=5*mag
    
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


#%%Functions to edit pillow images
def get_figure_sizes(cropbox, nrows, ncols, mag=1):
    w_crop_original = cropbox[2]-cropbox[0]
    h_crop_original = cropbox[3]-cropbox[1]    
    panel_size = int(w_crop_original/ncols)*2*mag, int(h_crop_original/nrows)*mag
    output_size = panel_size[0] * ncols, panel_size[1] * nrows
    return(panel_size, output_size)

def create_text_overlay(texts, output_size, panel_size, idxs):
    text_overlay = Image.new("RGBA", output_size)
    
    for j, (r, c) in enumerate(idxs):
        sub_extent=(c*panel_size[0], r*panel_size[1], (c+1)*panel_size[0], (r+1)*panel_size[1])
        text_im = get_text_image(texts[j], panel_size)
        text_overlay.paste(text_im, sub_extent)
    return(text_overlay)

def create_frames(files, texts, nrows, ncols, enhance=None):
    #specify output image size
    cropbox = (220, 40, 850, 850) #Original image 1000, 850    
    panel_size, output_size = get_figure_sizes(cropbox, nrows, ncols)

    #Load
    ims = [Image.open(f) for f in files]

    #Prepare
    idxs = [(int(i/ncols), i%ncols) for i in range(len(ims))]
    
    n_frames = [im.n_frames for im in ims]
    max_n = max(n_frames)
    start_frames = [max_n-n_frame for n_frame in n_frames]
    
    iterators=[ImageSequence.Iterator(im) for im in ims]
    
    #Create text overlay
    text_overlay = create_text_overlay(texts, output_size, panel_size, idxs)
    
    #Process
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
        if enhance is not None:
            #enhance = 1.4 is a good value
            stitched_image = ImageEnhance.Contrast(stitched_image).enhance(enhance)
            stitched_image = ImageEnhance.Brightness(stitched_image).enhance(enhance)
        stitched_image = stitched_image.convert("P", palette=Image.ADAPTIVE)
        out_frames.append(stitched_image)
    return(out_frames)

def save_frames(out_frames, out_path):
    out_gif_path = out_path % "gif"
    out_png_path = out_path % "png"
    
    out_frames[0].save(out_gif_path, format='GIF', 
              append_images=out_frames[1:], save_all=True, duration=130, loop=0,
              include_color_table=True)
    
    out_frames[-1].save(out_png_path, include_color_table=True)

#%%TODO
#Time + (Colorbar?)

#%%Settings
plot_trajectories=False
plot_inputs=True

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
datafol  = os.path.abspath(resource_filename("delta_aquifer", os.path.join("..", "data")))
traj_id  = os.path.join(datafol, "traj_id.csv")
text_path = os.path.join(datafol, "text_aid.csv")

#%%Get text aids
traj_len = 24

diff = get_changes(traj_id, traj_len)
texts = create_texts(diff, traj_len)

#%%Plot all trajectories
skip = [96, 216]

starts = np.arange(0, 12) * traj_len
starts[~np.isin(starts, skip)]
stops = starts+traj_len

if plot_trajectories:
    for start, stop in zip(starts, stops):
        files_mod = [i for i in files if get_model_id(i) in range(start, stop)]
        assert(len(files_mod)==traj_len)
        
        nrows, ncols = 3, 8 
        mod_idx = slice(start, stop)
        out_frames = create_frames(files_mod, texts[mod_idx], nrows, ncols)
        
        #Save
        out_path = os.path.join(out_fol, "SD_i{:03d}-{:03d}.%s".format(start, stop-1))
        save_frames(out_frames, out_path)

#%%Plot change per parameter
bad_chars = r"{}\\/"
pattern = re.compile('[%s]' % bad_chars)

if plot_inputs:
    df_text = pd.DataFrame(texts, columns=["par"])
    invalid_ids = np.concatenate([np.arange(s, s+traj_len) for s in skip])
    
    inps = list(np.unique(df_text["par"].str[1:-2])) #Strip of sign (dollar signs)
    inps.remove("bas") #Remove base run
    
    for inp in inps:
        print(inp)
        text_inp = df_text[df_text["par"].str.contains(inp, regex=False)]
        
        idx = text_inp.index
        idx = np.array(list(set(idx) - set(invalid_ids))) #Remove the values from idx that are in invalid_ids
        idx.sort()
        
        text_inp = list(text_inp.loc[idx]["par"])
        
        idx = np.concatenate([idx-1, idx])
        
        files_mod = [i for i in files if get_model_id(i) in idx]
        files_mod = files_mod[::2] + files_mod[1::2]
        
        nrows = 2
        ncols = int(len(files_mod)/nrows)
        
        text_inp = [""] * ncols + text_inp
        
        out_frames = create_frames(files_mod, text_inp, nrows, ncols)
        
        out_path = os.path.join(out_fol, "per_parameter", pattern.sub("", inp) + ".%s")
        save_frames(out_frames, out_path)
    
