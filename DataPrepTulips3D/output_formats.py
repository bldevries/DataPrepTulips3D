import mesaPlot as mp
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import numpy as np
import copy
from time import time
from PIL import Image
import os
import openexr_numpy as exr
from . import blackbody
from . import colormodels


def save_texture(array, output_path, do_exr=True):

    if array.ndim != 3:
        raise ValueError(f"Expected 3D array, got {array.ndim}D array")

    if np.max(array) > 1.:
        raise ValueError(f"Array not normed to one")

    if do_exr:
        float32_array = np.copy(array).astype(np.float32)
        exr.imwrite(output_path+".exr", float32_array)
    else:
        uint8_array = np.clip(array * 255, 0, 255).astype(np.uint8)

        # 4. Save as Image
        # Mode 'RGBA' handles the 4 channels
        img = Image.fromarray(uint8_array)#, mode='RGBA')
        # Save with optimization
        img.save(output_path+".png", optimize=True, compress_level=6)
    
def set_2D_data_to_R_channel(array):
    if array.ndim != 2:
        raise ValueError(f"Expected 2D array, got {array.ndim}D array")

    d_2d = np.zeros((array.shape[0], array.shape[1], 4))
    d_2d[:, :, 0] = array[:, :]/np.max(array)
    d_2d[:, :, 3] = 1.

    return d_2d, np.max(array)


def set_1D_data_to_R_channel(array):
    if array.ndim != 1:
        raise ValueError(f"Expected 1D array, got {array.ndim}D array")
    
    d_2d = np.zeros((1, len(array), 4))
    d_2d[0, :, 0] = array/np.max(array)
    d_2d[0, :, 3] = 1.

    return d_2d, np.max(array)

def set_color_map_chem_abun(array):
    CMAP_RGBA = [(0.25098039215686274, 0.7686274509803922, 1.0, 1.0), (0.39215686274509803, 0.7098039215686275, 0.9647058823529412, 1.0), (0.1607843137254902, 0.3843137254901961, 1.0, 1.0), \
            (0.7764705882352941, 1.0, 0.0, 1.0), (0.09411764705882353, 1.0, 1.0, 1.0), (0.0, 0.7843137254901961, 0.3254901960784314, 1.0), (0.4117647058823529, 0.9411764705882353, 0.6823529411764706, 1.0), \
            (0.5529411764705883, 0.43137254901960786, 0.38823529411764707, 1.0), (0.36470588235294116, 0.25098039215686274, 0.21568627450980393, 1.0), (0.7411764705882353, 0.7411764705882353, 0.7411764705882353, 1.0), \
            (0.1843137254901961, 0.0196078431372549, 0.5882352941176471, 1.0), (0.28627450980392155, 0.011764705882352941, 0.6274509803921569, 1.0), (0.3803921568627451, 0.0, 0.6549019607843137, 1.0), \
            (0.5294117647058824, 0.027450980392156862, 0.6509803921568628, 1.0), (0.7294117647058823, 0.2, 0.5333333333333333, 1.0), (0.8705882352941177, 0.3803921568627451, 0.39215686274509803, 1.0), \
            (0.9019607843137255, 0.4235294117647059, 0.3607843137254902, 1.0), (0.9294117647058824, 0.4745098039215686, 0.3254901960784314, 1.0), (0.996078431372549, 0.7176470588235294, 0.17647058823529413, 1.0)]

    colored_array = np.zeros((*array.shape,4))
    for r in range(array.shape[0]):
        for th in range(array.shape[1]):
            colored_array[r, th, :] = CMAP_RGBA[int(array[r, th])]#np.array([CMAP_RGBA[array[r, th]] for r, th in np.ndindex(a.shape)])
    
    return colored_array

def set_color_map_tulips(array, cmap_name='viridis', vmin=-10, vmax=10):

    if array.ndim != 1:
        raise ValueError(f"Expected 1D array, got {array.ndim}D array")
    
    # 1. Setup Colormap and Normalization
    # cmap = plt.get_cmap(cmap_name)
    p = mp.plot(rcparams_fixed=False)
    cmap = p.mergeCmaps([plt.cm.Purples_r, plt.cm.hot_r], [[0.0, 0.5], [0.5, 1.0]])
    
    # Create the Normalize object. 
    # If you pass specific vmin/vmax, it clamps values outside that range.
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # 2. Apply Normalization and Colormap
    # norm(array) scales data to [0, 1]
    # cmap(...) converts those scaled values to RGBA (floats 0.0-1.0)
    rgba_array = cmap(norm(array))
    
    # The result from cmap is shape (H, W, 4) with dtype float64
    # We need to check the shape just in case the input was masked or weird, 
    # but usually for a simple 2D array, this works directly.

    return rgba_array

def set_color_map_blackbody(array):
    ''' Turn the N-length 1d array into blackbody colors of shape Nx4'''

    if array.ndim != 1:
        raise ValueError(f"Expected 1D array, got {array.ndim}D array")

    colored = []

    for T in array:
        _c = [i/255 for i in colormodels.irgb_from_xyz(blackbody.blackbody_color(10**T))]
        # Add the alpha channel
        colored.append(_c + [1.])

    return np.array(colored)
