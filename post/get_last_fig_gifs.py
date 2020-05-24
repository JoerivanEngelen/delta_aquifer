# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:25:56 2020

@author: engelen
"""

from glob import glob
from PIL import Image, ImageSequence, ImageEnhance, ImageDraw
import os

#Source functions: https://gist.github.com/BigglesZX/4016539
def analyseImage(path):
    '''
    Pre-process pass over the image to determine the mode (full or additive).
    Necessary as assessing single frames isn't reliable. Need to know the mode 
    before processing all frames.
    '''
    im = Image.open(path)
    results = {
        'size': im.size,
        'mode': 'full',
    }
    try:
        while True:
            if im.tile:
                tile = im.tile[0]
                update_region = tile[1]
                update_region_dimensions = update_region[2:]
                if update_region_dimensions != im.size:
                    results['mode'] = 'partial'
                    break
            im.seek(im.tell() + 1)
    except EOFError:
        pass
    return results

def processImage(path):
    '''
    Iterate the GIF, extracting each frame.
    '''
    frames = []
    mode = analyseImage(path)['mode']
    
    im = Image.open(path)

    i = 0
    p = im.getpalette()
    last_frame = im.convert('RGBA')
    
    try:
        while True:
            #print("saving %s (%s) frame %d, %s %s" % (path, mode, i, im.size, im.tile))
            
            '''
            If the GIF uses local colour tables, each frame will have its own palette.
            If not, we need to apply the global palette to the new frame.
            '''
            if not im.getpalette():
                im.putpalette(p)
            
            new_frame = Image.new('RGBA', im.size)
            
            '''
            Is this file a "partial"-mode GIF where frames update a region of a different size to the entire image?
            If so, we need to construct the new frame by pasting it on top of the preceding frames.
            '''
            if mode == 'partial':
                new_frame.paste(last_frame)
            
            new_frame.paste(im, (0,0), im.convert('RGBA'))
#            new_frame.save('%s-%d.png' % (''.join(os.path.basename(path).split('.')[:-1]), i), 'PNG')
            frames.append(new_frame)

            i += 1
            last_frame = new_frame
            im.seek(im.tell() + 1)
    except EOFError:
        pass
    return(frames)

#%%Path management
gif_path = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Synth_Delta\results_30deltas\gifs\Natural\*.gif"
files = glob(gif_path)
out_folder = os.path.join(gif_path, "..", "..", "last_frame")
os.makedirs(out_folder, exist_ok=True)

fnames = [os.path.splitext(os.path.basename(f))[0] for f in files]
out_png_paths = [os.path.join(out_folder, fname+".png") for fname in fnames]
#%%Read
ims = [Image.open(f) for f in files]

#%%{Process}
cropbox = (220, 40, 850, 850)

#last_frames = [[frame for frame in iterator][-1] for iterator in iterators]
last_frames = [processImage(f)[-1] for f in files]

resized_frames = [frame.crop(cropbox) for frame in last_frames]


[frame.save(path, include_color_table=True) for path, frame in zip(out_png_paths, resized_frames)]