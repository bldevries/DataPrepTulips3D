# import bpy
import numpy as np
import sys, os
import pickle

# from time import time
# import matplotlib.pyplot as plt
# from matplotlib.colors import Normalize, LogNorm
# from scipy.interpolate import interp1d

import mesaPlot as mp
from scipy.interpolate import interp1d

from . import blackbody
from . import colormodels
from . import roche_lobe

# import Data1D.Data1D
from DataPrepTulips3D.output_formats import *

# import DataPrepTulips3D as DP

# from .blender_keys_prep import *
from . import blender_keys_prep as key

def save_to_texture(d, directory):
    # directory = os.path.join(directory, "VDB")
    # os.makedirs(directory, exist_ok=True)
    # We have several dimensions to save into textures. A texture has 2 dim (called 
    # Texture dimension below) and using a sequence of images we have a 3rd 
    # dimension (this dimension we call the file dimension).

    # First save all the data that depends on time and radius
    # Texture dimensions: nr_radial_points x 1
    # File dimension: time index

    for label in d[key.prof_labels]:
        print(f"Texture save of {label=}")
        # Make a directory for every profile
        prof_directory = os.path.join(directory, label)
        os.makedirs(prof_directory, exist_ok=True)
        # Get the index of the profile
        label_index = d[key.prof_labels].index(label)
        # Get the data
        d_label = d[key.data_prof_t_r][label_index, :, :]
        # Loop over time index
        for time_index in range(d_label.shape[0]):
            # Set the filename. For blender the numbering with dots seem important
            # so you can load it in as an image texture sequence using a # for image
            # number.
            filename_full = f"data.{str(time_index)}" #.zfill(4)
            # And add the filename convention to the dict
            d[key.dir_structure].update(\
                {key.data_t_r_filename: "data.0.exr"})

            # Calculate colors
            d_2d = set_color_map_tulips(d_label[time_index, :])
            # Add a dimension to the 1d array otherwise it will not save to a texture
            # -1 will inherit the dimension of the given array (d_label[time_index, :])
            # which is nr_radial_points
            d_2d = np.reshape(d_2d, (1, -1, 4))
            # Savind the texture:
            save_texture(d_2d, os.path.join(prof_directory, filename_full))#save_colormapped_texture

    # Now save the chemical profiles, which depend on time, radius and theta
    # Texture dimensions: nr_radial_points x theta_points
    # File dimension: time index
    # chem_abun_id[time, radius, theta]
    #
    shape = d[key.data_chem_abun_id].shape
    prof_directory = os.path.join(directory, "chem_abun")
    prof_directory_col = os.path.join(directory, "chem_abun_color")
    os.makedirs(prof_directory, exist_ok=True)
    os.makedirs(prof_directory_col, exist_ok=True)
    for t in range(shape[0]):
        d_2, d_max = set_2D_data_to_R_channel(d[key.data_chem_abun_id][t,:,:])
        filename_full = f"chem_data_{round(d_max,4)}_nrTh{shape[2]}_.{str(t)}" #.zfill(4)
        save_texture(d_2, os.path.join(prof_directory, filename_full))

        d_3 = set_color_map_chem_abun(d[key.data_chem_abun_id][t,:,:])
        filename_full = f"color_chem_data_nrTh{shape[2]}_.{str(t)}" #.zfill(4)
        save_texture(d_3, os.path.join(prof_directory_col, filename_full))
        d[key.dir_structure].update(\
                {key.chem_abun_filename: f"color_chem_data_nrTh{shape[2]}_.0.exr"})

    # Now save the Roche lobe shape grids (theta x phi, per frame), one per
    # star - see roche_lobe.py / blender_keys_prep.py. Unlike the chem_abun
    # block above (which normalizes each frame independently, fine for a
    # categorical abundance-species ID), this data is a physical radius that
    # must convert back to Rsun *consistently across the whole animation* -
    # so normalization uses ONE max value across ALL frames+directions
    # (like the plain data_t scalar path uses one max across the whole time
    # series), not a fresh per-frame max, and that single max_value is
    # recorded once in dir_structure rather than embedded per-frame in the
    # filename.
    for shape_key, filename_key, subdir_name in [
        (key.data_rochelobe_shape_1, key.rochelobe_shape_filename_1, "rochelobe_1"),
        (key.data_rochelobe_shape_2, key.rochelobe_shape_filename_2, "rochelobe_2"),
    ]:
        if shape_key not in d:
            continue
        rl_array = d[shape_key]  # (n_frames, n_theta, n_phi), Rsun
        rl_max = float(np.max(rl_array))
        n_theta_rl, n_phi_rl = rl_array.shape[1], rl_array.shape[2]
        rl_directory = os.path.join(directory, subdir_name)
        os.makedirs(rl_directory, exist_ok=True)
        for t in range(rl_array.shape[0]):
            d_2d = np.zeros((n_theta_rl, n_phi_rl, 4))
            d_2d[:, :, 0] = rl_array[t, :, :] / rl_max
            d_2d[:, :, 3] = 1.
            filename_full = f"rochelobe_data_nrTh{n_theta_rl}_nrPhi{n_phi_rl}_.{str(t)}"
            save_texture(d_2d, os.path.join(rl_directory, filename_full))
        d[key.dir_structure].update({
            filename_key: {
                # filename includes the subdir_name/ prefix (unlike the
                # bare-filename convention chem_abun uses, where the
                # consuming code hardcodes "chem_abun_color/" separately) -
                # keeps the dir_structure entry self-contained, one join
                # away from a full path via os.path.join(texture_dir, ...).
                "filename": f"{subdir_name}/rochelobe_data_nrTh{n_theta_rl}_nrPhi{n_phi_rl}_.0.exr",
                "max_value": rl_max,
                "n_theta": n_theta_rl,
                "n_phi": n_phi_rl,
                "n_frames": rl_array.shape[0],
            }
        })
        print(f"  Baked Roche lobe shape sequence: {subdir_name} (max={rl_max:.4f} Rsun, "
              f"{rl_array.shape[0]} frames, {n_theta_rl}x{n_phi_rl} grid)")

    # Now save the full Roche potential "landscape" (height-field grid over
    # the whole orbital plane) - see roche_lobe.roche_potential_surface_grid.
    # Values are dimensionless and uniformly negative (a potential, not a
    # physical distance like the shape grids above) - normalize by
    # max(abs(value)) across the WHOLE animation (one consistent divisor,
    # same reasoning as the Roche lobe shape grids and L1_dist_from_CM: the
    # Blender side needs to reconstruct a real height consistently frame to
    # frame, not have the normalization silently shift per frame).
    if key.data_rochelobe_potential_grid in d:
        pot_array = d[key.data_rochelobe_potential_grid]  # (n_frames, n_y, n_x), dimensionless Phi
        pot_absmax = float(np.max(np.abs(pot_array)))
        n_y_pot, n_x_pot = pot_array.shape[1], pot_array.shape[2]
        pot_directory = os.path.join(directory, "rochelobe_potential")
        os.makedirs(pot_directory, exist_ok=True)
        for t in range(pot_array.shape[0]):
            d_2d = np.zeros((n_y_pot, n_x_pot, 4))
            d_2d[:, :, 0] = pot_array[t, :, :] / pot_absmax
            d_2d[:, :, 3] = 1.
            filename_full = f"potential_data_nrX{n_x_pot}_nrY{n_y_pot}_.{str(t)}"
            save_texture(d_2d, os.path.join(pot_directory, filename_full))
        d[key.dir_structure].update({
            key.rochelobe_potential_filename: {
                "filename": f"rochelobe_potential/potential_data_nrX{n_x_pot}_nrY{n_y_pot}_.0.exr",
                "max_value": pot_absmax,
                "n_x": n_x_pot,
                "n_y": n_y_pot,
                "n_frames": pot_array.shape[0],
                "x_range": d[key.potential_x_range],
                "y_range": d[key.potential_y_range],
            }
        })
        print(f"  Baked Roche potential surface grid (absmax={pot_absmax:.4f}, "
              f"{pot_array.shape[0]} frames, {n_x_pot}x{n_y_pot} grid)")

    # Now save the L2/L3 "outer" equipotential rings (see
    # roche_lobe.equipotential_boundary_xy) - each is a (n_frames, n_theta, 2)
    # array of (x,y) Rsun positions (NOT a scalar radius-from-a-fixed-center
    # like the L1 lobe shape grids above - the merged shape isn't star-shaped
    # from any single fixed center for realistic mass ratios, so both x AND y
    # are baked directly per ring point). Packed R=x, G=y (each independently
    # abs-max normalized - x and y have different natural scales) in a single
    # (1, n_theta, 4) image per frame, same per-frame-EXR-sequence convention
    # as the shape grids/potential grid above.
    for ring_key, subdir_name, filename_key in (
        (key.data_rochelobe_l2_ring, "rochelobe_l2_ring", key.rochelobe_l2_ring_filename),
        (key.data_rochelobe_l3_ring, "rochelobe_l3_ring", key.rochelobe_l3_ring_filename),
    ):
        if ring_key not in d:
            continue
        ring_array = d[ring_key]  # (n_frames, n_theta, 2), Rsun
        x_absmax = float(np.max(np.abs(ring_array[:, :, 0])))
        y_absmax = float(np.max(np.abs(ring_array[:, :, 1])))
        n_theta_ring = ring_array.shape[1]
        ring_directory = os.path.join(directory, subdir_name)
        os.makedirs(ring_directory, exist_ok=True)
        for t in range(ring_array.shape[0]):
            d_2d = np.zeros((1, n_theta_ring, 4))
            d_2d[0, :, 0] = ring_array[t, :, 0] / x_absmax
            d_2d[0, :, 1] = ring_array[t, :, 1] / y_absmax
            d_2d[0, :, 3] = 1.
            filename_full = f"ring_data_nrTheta{n_theta_ring}_.{str(t)}"
            save_texture(d_2d, os.path.join(ring_directory, filename_full))
        d[key.dir_structure].update({
            filename_key: {
                "filename": f"{subdir_name}/ring_data_nrTheta{n_theta_ring}_.0.exr",
                "x_max_value": x_absmax,
                "y_max_value": y_absmax,
                "n_theta": n_theta_ring,
                "n_frames": ring_array.shape[0],
            }
        })
        print(f"  Baked equipotential ring: {subdir_name} (x_absmax={x_absmax:.4f}, "
              f"y_absmax={y_absmax:.4f} Rsun, {ring_array.shape[0]} frames, {n_theta_ring} points)")

    # Now save all data only dependent on time
    data_dir_dict = {}
    for index, (_key, data) in enumerate(d[key.data_t].items()):
        if _key == "logTeff":
            filename_full = "colored_"+_key
            # d_2d = np.reshape(data, (1, -1))
            d_2d = set_color_map_blackbody(data)
            d_2d = np.reshape(d_2d, (1, -1, 4))
            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({"colored_"+_key: {"filename":filename_full+".exr", "max_value":d_max}})

            d_2d, d_max = set_1D_data_to_R_channel(data)
            filename_full = f"data_{round(d_max,4)}_"+_key
            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({_key: {"filename":filename_full+".exr", "max_value":d_max}})
        elif _key == "logTeff_color":
            print("Colored already done")
        elif _key in ["log_abs_mdot", "lg_mstar_dot_1", "lg_mstar_dot_2"]:
            if not np.all(data < 0.):
                raise ValueError(f"Not all massloss rates have a negative exponent! Check your code.")

            # We make all values positive before normalizing, not a nice thing to do!
            d_2d, d_max = set_1D_data_to_R_channel(abs(data))
            d_max = -1* d_max
            filename_full = f"data_{round(d_max,4)}_"+_key

            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({_key: {"filename":filename_full+".exr", "max_value":d_max}})
        elif _key in ("phi_l1_softened", "phi_l2_softened", "phi_l3_softened"):
            # Uniformly negative (like the main potential grid) - NOT
            # sign-swinging like L1_dist_from_CM, but the generic branch's
            # plain max(data) normalization is still wrong here: max(data)
            # for all-negative data is the value CLOSEST to zero, so
            # data/max(data) can land well outside [-1,1] and even flip
            # sign for the more-negative frames. abs-max keeps every frame
            # in [-1,1] with sign preserved, same fix as the main potential
            # grid's own baking block above.
            d_absmax = float(np.max(np.abs(data)))
            d_2d = np.zeros((1, len(data), 4))
            d_2d[0, :, 0] = data / d_absmax
            d_2d[0, :, 3] = 1.
            filename_full = f"data_{round(d_absmax,4)}_"+_key
            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({_key: {"filename":filename_full+".exr", "max_value":d_absmax}})
        elif _key in ("L1_dist_from_CM", "L2_dist_from_CM", "L3_dist_from_CM"):
            # Unlike log_abs_mdot (always negative) or everything in the
            # generic branch below (always positive), this genuinely swings
            # sign (positive on star 1's side of the CM, negative on star
            # 2's) - normalizing by plain max(data) (the generic branch's
            # approach) can drive most frames to a large NEGATIVE normalized
            # value (e.g. -24x, if one frame happens to sit near +0 while
            # most others are strongly negative) - safe as a linear
            # transform in principle, but empirically caused badly wrong
            # sampled values in Blender (well outside plain reconstruction
            # error, likely the image texture pipeline not tolerating such
            # a wide/lopsided value range gracefully). Normalizing by
            # max(abs(data)) instead keeps every frame's value within
            # [-1, 1] - still a signed value, sign preserved exactly,
            # reconstructed the same way (normalized * max_value) - just a
            # safer divisor.
            d_absmax = float(np.max(np.abs(data)))
            d_2d = np.zeros((1, len(data), 4))
            d_2d[0, :, 0] = data / d_absmax
            d_2d[0, :, 3] = 1.
            filename_full = f"data_{round(d_absmax,4)}_"+_key
            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({_key: {"filename":filename_full+".exr", "max_value":d_absmax}})
        else:
            d_2d, d_max = set_1D_data_to_R_channel(data)
            filename_full = f"data_{round(d_max,4)}_"+_key
            save_texture(d_2d, os.path.join(directory, filename_full))
            data_dir_dict.update({_key: {"filename":filename_full+".exr", "max_value":d_max}})


        # if key == "logTeff":
        #     filename_full = f"data_{key}"+ext
        #     # d_2d = np.reshape(data, (1, -1))/np.max(data)
        #     d_2d = np.zeros((1, len(data), 4))
        #     d_2d[0, :, 0] = data/np.max(data)
        #     d_2d[0, :, 3] = 1.
        
        #     save_texture(d_2d, os.path.join(directory, filename_full))
        # elif key == "logTeff_color":
        #     filename_full = f"data_{key}"+ext
        #     # print(data)
        #     d_2d = np.reshape(data, (1, -1, 4))/np.max(data)
        #     save_texture(d_2d, os.path.join(directory, filename_full))
    d[key.dir_structure].update(\
            {key.data_t_filename_and_max: data_dir_dict})


    return d

def save_to_pickle(data_dict, filename_path_no_ext, add_ext=None):
    '''Saves a dict into a pickle file'''
    # dbfile = open(filepath, 'ab')
    if add_ext == None:
        filename_path = filename_path_no_ext
    else:
        filename_path = filename_path_no_ext+add_ext

    if not os.path.isfile(filename_path):
        print(f"Saving pickle file to: ", filename_path)
        dbfile = open(filename_path, 'wb')
        pickle.dump(data_dict, dbfile)
        dbfile.close()
    else:
        print("File not found: ", filename_path)

def load_from_pickle(filepath):
    '''Reads a dict from a pickle file'''
    dbfile = open(filepath, 'rb')    
    data_dict = pickle.load(dbfile)
    return data_dict

def save(data_dict, directory, filename_base):
    save_to_pickle(data_dict, directory, filename_base)
    save_to_texture(d, directory, filename_base)



def convertMesaData(mesa_LOGS_directory, t_resolution, r_resolution,\
                    save_to_dir, pickle_filename = "MESA_data_dict.pkl",\
                    time_scale_type="log_to_end",\
                    filename_history = None, verbose_timing = False, \
                    profiles=[], r_grid_name="mass", is_binary=False, binary_nr=1,\
                    n_theta_rochelobe=128, n_phi_rochelobe=24,\
                    n_x_potential=48, n_y_potential=48, potential_softening=0.1):

    # binary_nr = 1 or 2

    pickle_file = os.path.join(save_to_dir, pickle_filename)

    if os.path.isfile(pickle_file):
        print("Dir exists, doing nothing")
        return load_from_pickle(pickle_file)
    else:
        d = loadMesaData(mesa_LOGS_directory, t_resolution, r_resolution,\
                        time_scale_type,\
                        filename_history, verbose_timing, \
                        profiles, r_grid_name,\
                        n_theta_rochelobe, n_phi_rochelobe,\
                        n_x_potential, n_y_potential, potential_softening)
        d.update({key.dir_structure:{"pickle_filename": pickle_filename}})
        d.update({key.is_binary: is_binary})
        d.update({key.binary_nr: binary_nr})

        text_dir = "textures"
        texture_dir = os.path.join(save_to_dir, text_dir)
        d[key.dir_structure].update({key.texture_dir: text_dir})
        d = save_to_texture(d, texture_dir)

        save_to_pickle(d, pickle_file)

        return d

def loadMesaData(mesa_LOGS_directory, t_resolution, r_resolution,\
                time_scale_type="log_to_end",\
                filename_history = None, verbose_timing = False, \
                profiles=[], r_grid_name="mass",\
                n_theta_rochelobe=24, n_phi_rochelobe=24,\
                n_x_potential=48, n_y_potential=48, potential_softening=0.1):
    """
    Loads MESA data, sets it to a specified resolution and saves it 
    into a dictionairy. All values are set to the same grid (given by
    t_resolution, r_resolution). The returned dict. contains profile data 
    that is a function of radius and time (data_prof_t_r) and data that 
    only depends on time (data_t). It also contains r_max, which gives
    the maximum radius of profile data (depending on the choosen r_grid_name)
    and can be used with r_resolution to recontruct the radial values.
    Note that all data is forced to have a value at the radial value at 0.

    Parameters:
    mesa_LOGS_directory (str): directory of the LOGS file
    t_resolution (int): the amount of time steps you want the output to have
    r_resolution (int): the amount of radial steps you want the output to have
    filename_history (str): if the history file has a non-standard filename, indicate it here
    profiles (list of str): optionally you can list the profiles you want yourself
    r_grid_name (str): the grid profile you want to use
    verbose_timing (bool): if you want more output, set to true

    Returns:
    mesa_data: dict containing the MESA data
    """

    if len(profiles) == 0:
        profiles = ['mass', 'logT', 'logRho', 'he4'] + ["en"]

    if verbose_timing: _ = time()

    # Create MESA object
    m = mp.MESA() 

    # Figure out the history file to use
    if filename_history:
        f = os.path.join(mesa_LOGS_directory, filename_history)
    else:
        f = os.path.join(mesa_LOGS_directory, "history.data")

    # Load the mesa history data
    m.loadHistory(filename_in = f)

    _age_indices = np.round(np.linspace(0, len(m.hist.star_age) - 1, t_resolution)).astype(int)
    _age_indices = rescale_time(_age_indices, m, time_scale_type=time_scale_type)
    new_age_grid = m.hist.star_age[_age_indices]

    if verbose_timing: print("Timing load mesa file: ", time()-_)

    t_resolution_orig = len(m.hist.star_age)
    print("Original/new time resolution: ", t_resolution_orig, len(new_age_grid))

    # Load the Teff data of the star
    print("Loading Teff")
    logTeff_values, logTeff_colors = loadMesaTeffData(m, age_indices=_age_indices)

    print("Loading Profiles")
    data_array_prof, chem_array, R_star_from_grid, R_star_from_grid_Rsun, elem_list, chem_abun_id \
                                         = loadMesaProfile(m , mesa_LOGS_directory, \
                                                    r_resolution = r_resolution,\
                                                    profile_names=profiles, \
                                                    r_grid_name=r_grid_name, \
                                                    age_indices = _age_indices)

    print("Loading History")
    _data_t = loadMesaHistory(m, age_indices=_age_indices)

    _data_t.update({"logTeff": logTeff_values, "logTeff_color": logTeff_colors, "Rmax": R_star_from_grid[0,:], "Rmax_Rsun": R_star_from_grid_Rsun[0,:], "age": new_age_grid})

    # Roche lobe geometry (physically-accurate equipotential surface, not
    # just the Eggleton-formula scalar MESA already gives us in rl_1/rl_2) -
    # only possible/meaningful for a binary with both masses and the orbital
    # separation available. See roche_lobe.py for the actual math. Computed
    # for BOTH stars in every pass (regardless of whether mesa_LOGS_directory
    # is star 1's or star 2's own LOGS dir) since MESA's binary history
    # columns (masses, separation, r1, r2) are already replicated in both
    # stars' own history files - same reasoning as why r1 AND r2 already end
    # up in _data_t above regardless of which star this is.
    data_rochelobe_shape_1 = None
    data_rochelobe_shape_2 = None
    if all(k in _data_t for k in ("star_1_mass", "star_2_mass", "binary_separation", "r1")):
        m1 = _data_t["star_1_mass"]
        m2 = _data_t["star_2_mass"]
        sep = _data_t["binary_separation"]  # Rsun
        q1 = m2 / m1  # star 1's own lobe: companion(2)/self(1)
        q2 = m1 / m2  # star 2's own lobe: companion(1)/self(2)

        r_grid_1, x_l1 = roche_lobe.roche_lobe_radius_grid(q1, n_theta_rochelobe, n_phi_rochelobe)
        r_grid_2, _ = roche_lobe.roche_lobe_radius_grid(q2, n_theta_rochelobe, n_phi_rochelobe)
        data_rochelobe_shape_1 = r_grid_1 * sep[:, None, None]  # Rsun
        data_rochelobe_shape_2 = r_grid_2 * sep[:, None, None]  # Rsun

        # L1's position relative to the center of mass, in Rsun, signed along
        # the same X axis convention the two stars already use (star 1 sits
        # at -r1, star 2 at +r2 - see add_geo_nodes_profile's distance_sign
        # in the TULIPS-3D addon): x_l1 is measured from star 1 toward star
        # 2, so L1's absolute position is star 1's own position (-r1) plus
        # that offset.
        L1_dist_from_1 = x_l1 * sep  # Rsun, from star 1 toward star 2
        _data_t.update({"L1_dist_from_CM": -_data_t["r1"] + L1_dist_from_1})
        print("  Added Roche lobe shape grids and L1_dist_from_CM")

        # Potential AT L1, using the same softening as the potential-surface
        # height field (roche_potential_surface_grid below) so a curve drawn
        # at this exact height lands precisely on that surface's own L1
        # elevation - not the true unsoftened critical potential, which
        # would be a hair off from where the softened landscape actually
        # sits. Dimensionless, uniformly negative (like the potential grid
        # itself) - baked via the same abs-max convention as L1_dist_from_CM
        # is NOT appropriate here (that one swings sign; this one doesn't),
        # see the "phi_l1_softened" branch in save_to_texture below instead.
        _data_t.update({
            "phi_l1_softened": roche_lobe.roche_potential_softened(x_l1, 0., 0., q1, potential_softening)
        })
        print("  Added Phi(L1) (softened)")

        # L2/L3 (outer Lagrange points, beyond star 2 / star 1 respectively
        # along the line joining the stars) - exact critical points via
        # root-finding (unlike L4/L5, there's no closed form). Computed for
        # ALL frames up front so the potential-surface domain (x_range/
        # y_range below) can be widened to comfortably contain them (Ben:
        # "increase the potential wireframe further in the +/-x directions",
        # later "make the grid cover the L2 and L3 curves") - unlike L4/L5's
        # fixed position, x_l2/x_l3 depend on q1, which varies over the run
        # as mass transfer proceeds, so the margin has to be based on the
        # actual worst-case extent seen in THIS run, not a universal
        # constant.
        x_l2 = roche_lobe.find_L2(q1)
        x_l3 = roche_lobe.find_L3(q1)
        phi_l1_true = roche_lobe.roche_potential(x_l1, 0., 0., q1)  # unsoftened - for the equipotential backoff below
        phi_l2_true = roche_lobe.roche_potential(x_l2, 0., 0., q1)
        phi_l3_true = roche_lobe.roche_potential(x_l3, 0., 0., q1)
        _data_t.update({
            "phi_l2_softened": roche_lobe.roche_potential_softened(x_l2, 0., 0., q1, potential_softening),
            "phi_l3_softened": roche_lobe.roche_potential_softened(x_l3, 0., 0., q1, potential_softening),
        })
        print("  Added Phi(L2)/Phi(L3) (softened)")

        # L2/L3's position relative to the CM, in Rsun - same signed
        # convention as L1_dist_from_CM above, for a small marker object at
        # each (Ben: "make a marker for the L2 and L3 locations just as for
        # the L1"). Both are always on the same side of the CM regardless
        # of q (x_l2>1 always -> always beyond star 2's side; x_l3<0 always
        # -> always beyond star 1's side) so, unlike L1_dist_from_CM, these
        # never swing sign - still baked via the same abs-max branch below
        # for consistency/safety, not because it's strictly required here.
        _data_t.update({
            "L2_dist_from_CM": x_l2 * sep - _data_t["r1"],
            "L3_dist_from_CM": x_l3 * sep - _data_t["r1"],
        })
        print("  Added L2_dist_from_CM/L3_dist_from_CM")

        # L2/L3 "outer" equipotential rings - a closed contour just inside
        # (Phi a small fraction below) each critical value, since the EXACT
        # critical potential is itself a degenerate/pinched contour (same
        # reasoning as the theta=0 pole special-case in
        # roche_lobe_radius_grid, just pinching at the L2/L3 point instead
        # of L1). Uses roche_lobe.equipotential_boundary_xy per frame (not
        # vectorizable across frames like the rest of this module - it does
        # a flood-fill + boundary trace per frame - but fast, a few ms per
        # frame even at n_grid=140). Falls back to the previous frame's ring
        # on a rare failure (e.g. an unexpectedly extreme mass ratio pushing
        # the region to touch the search-grid edge) rather than crashing the
        # whole data-prep run over one frame.
        #
        # Computed HERE (before x_range/y_range below) using a generous,
        # independent search domain - not the final display domain, which
        # is only known once the rings' own extent is known (a chicken-and-
        # egg problem otherwise). n_theta_equipot=96 (up from an initial 48,
        # Ben: "increase the amount of points used for the L1/2/3 curves so
        # that they are smoother") - the underlying boundary trace already
        # has much finer (grid-resolution, n_grid=140) detail than either
        # point count resamples to, so this is a pure smoothness win, not
        # limited by the trace itself.
        n_theta_equipot = 200
        backoff_frac = 0.03  # small: stays close to the true critical value while avoiding its exact (degenerate) pinch
        ring_search_pad = 1.0  # generous margin beyond x_l2/x_l3 for the search grid itself (independent of the final display x_pad below)
        ring_domain_x = (float(x_l3.min() - ring_search_pad), float(x_l2.max() + ring_search_pad))
        ring_domain_y = (-3.0, 3.0)  # generous fixed range - verified safe across q=0.1-10 (this project's actual MESA range)
        (_, l4_y), (l5_x, l5_y) = roche_lobe.l4_l5_positions()
        seed_xy = (0.5, l4_y)  # L4 - the domain's own potential max away from the two singularities, always "inside"
        n_frames_binary = len(q1)
        ring_l2 = np.zeros((n_frames_binary, n_theta_equipot, 2))
        ring_l3 = np.zeros((n_frames_binary, n_theta_equipot, 2))
        for t in range(n_frames_binary):
            phi_target_l2 = phi_l2_true[t] - backoff_frac * abs(phi_l2_true[t] - phi_l1_true[t])
            phi_target_l3 = phi_l3_true[t] - backoff_frac * abs(phi_l3_true[t] - phi_l1_true[t])
            try:
                ring_l2[t] = roche_lobe.equipotential_boundary_xy(
                    q1[t], phi_target_l2, ring_domain_x, ring_domain_y, seed_xy, n_theta=n_theta_equipot)
            except ValueError as e:
                print(f"  X L2 ring failed at frame {t} ({e}) - reusing previous frame")
                ring_l2[t] = ring_l2[t-1] if t > 0 else 0.
            try:
                ring_l3[t] = roche_lobe.equipotential_boundary_xy(
                    q1[t], phi_target_l3, ring_domain_x, ring_domain_y, seed_xy, n_theta=n_theta_equipot)
            except ValueError as e:
                print(f"  X L3 ring failed at frame {t} ({e}) - reusing previous frame")
                ring_l3[t] = ring_l3[t-1] if t > 0 else 0.
        data_rochelobe_l2_ring = ring_l2 * sep[:, None, None]  # Rsun (scales both x and y)
        data_rochelobe_l3_ring = ring_l3 * sep[:, None, None]
        print(f"  Added L2/L3 equipotential rings ({n_frames_binary} frames, {n_theta_equipot} points each)")

        # Full Roche potential "landscape" (a height-field surface plot of
        # the whole potential, not just the single Roche-lobe equipotential
        # contour above) - a visualization aid, prototype stage (Ben: "using
        # the MESA data, reproduce the 3d plotted full roche potential in
        # Blender"). Shared between both stars (like the Roche lobe shape
        # grids), computed once here using q1 (star 1 at the origin, the
        # convention roche_potential_surface_grid/l4_l5_positions assume).
        # x_range/y_range cover the ACTUAL computed extent of the L2/L3
        # rings above (which, thanks to equipotential_boundary_xy's small
        # backoff, can reach marginally past x_l2/x_l3 themselves) unioned
        # with the L4/L5-based extent, so the grid is guaranteed - by
        # construction, not by a hopeful margin guess - to fully contain
        # both rings (Ben: "make the grid cover the L2 and L3 curves").
        margin_y = 1.3  # L4/L5's own margin, per Ben's original request - kept as a floor even where the rings don't reach this far
        pad = 0.3  # a bit further than whichever (ring or L4/L5) extent ends up larger, in both x and y
        x_min = min(float(ring_l2[:, :, 0].min()), float(ring_l3[:, :, 0].min()), float(x_l3.min()))
        x_max = max(float(ring_l2[:, :, 0].max()), float(ring_l3[:, :, 0].max()), float(x_l2.max()))
        y_min = min(float(ring_l2[:, :, 1].min()), float(ring_l3[:, :, 1].min()), l5_y * margin_y)
        y_max = max(float(ring_l2[:, :, 1].max()), float(ring_l3[:, :, 1].max()), l4_y * margin_y)
        x_range = (x_min - pad, x_max + pad)
        y_range = (y_min - pad, y_max + pad)
        potential_grid = roche_lobe.roche_potential_surface_grid(
            q1, x_range, y_range, n_x=n_x_potential, n_y=n_y_potential, softening=potential_softening)[0]
        data_rochelobe_potential_grid = potential_grid  # dimensionless (G(M1+M2)=1, a=1) - NOT yet scaled to physical units
        print(f"  Added Roche potential surface grid (x_range={x_range}, y_range={y_range})")
    else:
        print("  X Skipped Roche lobe geometry (not enough binary data present)")
        data_rochelobe_potential_grid = None
        data_rochelobe_l2_ring = None
        data_rochelobe_l3_ring = None
        x_range, y_range = None, None

    print(_age_indices)
    print(type(_age_indices))
    print(type(_age_indices[0]))
    print()

    result = {\
            "info":"",\
            "MESA_file": mesa_LOGS_directory,\
            "filename_history":filename_history,\
            "r_label": r_grid_name,\
            "age": new_age_grid,\
            "age_indices": list(_age_indices),\
            key.prof_labels: profiles,\
            key.chem_elem_labels: elem_list,\
            key.t_resolution: t_resolution, \
            key.t_resolution_orig: t_resolution_orig, \
            key.time_scale_type: time_scale_type, \
            key.r_resolution: r_resolution,\
            key.data_prof_t_r: data_array_prof,\
            key.data_chem_t_r: chem_array,\
            key.data_chem_abun_id: chem_abun_id,\
            key.data_r_max: R_star_from_grid,\
            key.data_r_max_Rsun: R_star_from_grid_Rsun,\
            key.data_t: _data_t\
            }
    if data_rochelobe_shape_1 is not None:
        result.update({
            key.data_rochelobe_shape_1: data_rochelobe_shape_1,
            key.data_rochelobe_shape_2: data_rochelobe_shape_2,
            key.nr_theta_points_rochelobe: n_theta_rochelobe,
            key.nr_phi_points_rochelobe: n_phi_rochelobe,
        })
    if data_rochelobe_potential_grid is not None:
        result.update({
            key.data_rochelobe_potential_grid: data_rochelobe_potential_grid,
            key.nr_x_points_potential: n_x_potential,
            key.nr_y_points_potential: n_y_potential,
            key.potential_x_range: x_range,
            key.potential_y_range: y_range,
            key.data_rochelobe_l2_ring: data_rochelobe_l2_ring,
            key.data_rochelobe_l3_ring: data_rochelobe_l3_ring,
            key.nr_theta_points_equipotential: n_theta_equipot,
        })
    return result

def loadMesaProfile(m, mesa_LOGS_directory, profile_names, \
                    r_resolution, 
                    r_grid_name, age_indices):
    """Reads in a profile data from a mesa file"""

    def find_profile(m, mesa_LOGS_directory, time_ind=0):
        model_number = m.hist.model_number[time_ind]
        m.loadProfile(num=model_number, f=mesa_LOGS_directory, silent=True)
        return m.prof
    
    _prof = find_profile(m, mesa_LOGS_directory)
    r_grid = _prof.data[r_grid_name][:]
    print("Profile data keys:", _prof.data.dtype.names)
    print()

    chem_prof = True
    if chem_prof:
        p = mp.plot(rcparams_fixed=False)
        elem_list = p._listAbun(_prof)
        # Remove "neut", "prot", and ionized hydrogen from the list of isotopes
        if "neut" in elem_list:
            elem_list.remove("neut")
        if "prot" in elem_list:
            elem_list.remove("prot")
        if "h1_1" in elem_list:
            elem_list.remove("h1_1")


    # DOING THE PROFILES
    print("Working on profile data")

    # This will contain the MESA data in the new resolutions
    data_array = np.zeros((len(profile_names), len(age_indices), r_resolution))

    # This will contain the max Radius at each time index
    R_star_from_grid = np.zeros((len(profile_names), len(age_indices)))
    R_star_from_grid_Rsun = np.zeros((len(profile_names), len(age_indices)))

    for i_prof_name, pname in enumerate(profile_names):
        for i, t in enumerate(age_indices):
            if pname == "en":
                _R_max, _prop = loadMesaEnergyData(m, t, r_resolution)
            else:
                prof = find_profile(m, mesa_LOGS_directory, time_ind=t)
                _r = prof.data[r_grid_name][:]
                _r_Rsun = 10**prof.data["logR"][:] # log_R ! log10 radius in Rsun units

                _prop = prof.data[pname][:]

                # If the order is descending, flip the arrays
                if _r[0] > _r[-1]: 
                    _r = np.flip(_r)
                    _r_Rsun = np.flip(_r_Rsun)
                    _prop = np.flip(_prop)
                # If there is no element at r=0.0, add it
                if _r[0] != 0.: 
                    _r = np.concatenate([[0.],_r])
                    _r_Rsun = np.concatenate([[0.],_r_Rsun])
                    _prop = np.concatenate([[_prop[0]],_prop])

                _new_r = np.linspace(min(_r), max(_r), num=r_resolution)

                # print(i_prof_name, prop.shape)
                f = interp1d(_r, _prop, bounds_error=False, fill_value=np.nan)
                _prop = f(_new_r)
                _R_max = max(_r)

                f = interp1d(_r, _r_Rsun, bounds_error=False, fill_value=np.nan)
                _r_Rsun = f(_new_r)
                _R_max_Rsun = max(_r_Rsun)

            R_star_from_grid[i_prof_name,i] = _R_max
            R_star_from_grid_Rsun[i_prof_name,i] = _R_max_Rsun
            data_array[i_prof_name,i,:] = _prop

    # DOING THE CHEM
    print("Working on Chem. abundance data")

    # This will contain the MESA data in the new resolutions
    chem_array = np.zeros((len(elem_list), len(age_indices), r_resolution))

    # This will contain the max Radius at each time index
    # R_star_from_grid = np.zeros((len(elem_list), len(age_indices)))

    for i_prof_name, pname in enumerate(elem_list):
        for i, t in enumerate(age_indices):
            prof = find_profile(m, mesa_LOGS_directory, time_ind=t)
            _r = prof.data[r_grid_name][:]
            _prop = prof.data[pname][:]

            # If the order is descending, flip the arrays
            if _r[0] > _r[-1]: 
                _r = np.flip(_r)
                _prop = np.flip(_prop)
            # If there is no element at r=0.0, add it
            if _r[0] != 0.: 
                _r = np.concatenate([[0.],_r])
                _prop = np.concatenate([[_prop[0]],_prop])

            _new_r = np.linspace(min(_r), max(_r), num=r_resolution)

            # print(i_prof_name, prop.shape)
            f = interp1d(_r, _prop, bounds_error=False, fill_value=np.nan)
            _prop = f(_new_r)
            _R_max = max(_r)

            # R_star_from_grid[i_prof_name,i] = _R_max
            chem_array[i_prof_name,i,:] = _prop

    
    # We will sum away the different chemical profile labels. So we need an array that has indices 
    # abundances_cummulative[chem.prof., time index, radial index].
    abundances_cummulative = np.zeros(chem_array.shape)#(_chem_data.shape[1], _chem_data.shape[2]))
    # We need to know the number of theta indices
    nr_Th = 50#ob[key.nr_theta_points]
    # Iterate over the time indices
    for t in range(chem_array.shape[1]):
        for r in range(chem_array.shape[2]):
            abun = chem_array[:, t, r] * nr_Th # The abundances at t, r scaled by the nr of theta points
            # abundances = np.array([i for i in v[:, r_index]])*nr_Th
            # abundances = v*nr_Th # We multiply with nr_Th since it will be in ratio to the theta indices
            abundances_cummulative[:, t, r] = np.array([np.sum(abun[0:i+1]) for i in range(len(abun))]) 

            # abundances_cummulative[t, :] = np.array([
            #     np.array([np.sum(abundances[0:i+1, r_index]) 
            #         for i in range(len(abundances[:, r_index]))]) 
            #         for r_index in range(len(abundances[0, :]))])
    

    # list_color = []
    # cmap = CMAP_BASE
    # CMAP_DEFAULT = cmr.get_sub_cmap("cmr.pride", 0, 0.8)
    # cmap = plt.get_cmap(CMAP_DEFAULT, len(labels))

    # We make a new array with indices chem_abun_id[time, radius, theta]
    # And it contains the chemical species ID of the element filling up the (r, th)
    # of the star pie
    chem_abun_id = np.zeros((abundances_cummulative.shape[1], abundances_cummulative.shape[2], nr_Th))
    for t in range(chem_array.shape[1]):
        for r_index in range(chem_array.shape[2]):
            for th_index in range(nr_Th):
                v = abundances_cummulative[:, t, :]
                # r_index = mesh.attributes['vert_col_radial_index'].data[vert_i_mesh].value
                # th_index = mesh.attributes['vert_col_th_index'].data[vert_i_mesh].value
                chem_abun_id[t, r_index, th_index] = np.searchsorted(v[:, r_index], th_index)
                # vert_color = CMAP_RGBA[idx]
                # list_color.append(vert_color)
        # abundances_cummulative

    return data_array, chem_array, R_star_from_grid, R_star_from_grid_Rsun, elem_list, chem_abun_id



def loadMesaEnergyData(m, time_index, r_resolution, verbose_timing=False):
    '''Reads in Energy production data from a mesa file'''

    # We need these indices for the burning data
    qtop = "burn_qtop_"
    qtype = "burn_type_"

    # A check if data is avaliable
    try:
        m.hist.data[qtop + "1"]
    except ValueError:
        raise KeyError(
            "No field " + qtop + "* found, add mixing_regions 40 and burning_regions 40 to your history_columns.list")

    sm = m.hist.star_mass

    time_indices = len(sm)

    if verbose_timing: _ = time()

    num_burn_zones = int([xx.split('_')[2] for xx in m.hist.data.dtype.names if qtop in xx][-1])
    # Per time stamp we will have radii and values, which we list in these variables:
    list_r = [np.abs(m.hist.data[qtop + str(region)][time_index] * sm[time_index]) for region in range(1, num_burn_zones + 1)]
    list_E = [m.hist.data[qtype + str(region)][time_index] for region in range(1, num_burn_zones + 1)]

    # We make one data array
    _d = np.array([list_r, list_E])
    # We need to remove duplicates at the end where value==-9999
    _mask = _d[1] != -9999
    # Our cleaned up arrays containing the radii and burning values
    _r, _prop = _d[0][_mask], _d[1][_mask] # Now you have data you can interpolate f = interp1d(r, v, kind='cubic')

    # If the order is descending, flip the arrays
    if _r[0] > _r[-1]: 
        _r = np.flip(_r)
        _prop = np.flip(_prop)
    # If there is no element at r=0.0, add it
    if _r[0] != 0.: 
        _r = np.concatenate([[0.],_r])
        _prop = np.concatenate([[_prop[0]],_prop])

    _new_r = np.linspace(min(_r), max(_r), num=r_resolution)
    f = interp1d(_r, _prop, bounds_error=False, fill_value=np.nan)
    _prop = f(_new_r)
    _R_max = max(_r)

    return _R_max, _prop #v, r# E_value, E_grid


def loadMesaTeffData(mesa_object, age_indices, verbose_timing=False):
    '''Reads in Teff data from a mesa file'''

    sm = mesa_object.hist.star_mass
    time_indices = len(sm)
    logTeff_values = []
    logTeff_colors = []

    for t in age_indices:
        T = mesa_object.hist.log_Teff[t]
        logTeff_values.append(T)

        # And make a colour of the logTeff
        _c = [i/255 for i in colormodels.irgb_from_xyz(blackbody.blackbody_color(10**T))]

        # Add the alpha channel
        logTeff_colors.append(_c + [1.])

    return np.array(logTeff_values), np.array(logTeff_colors)

def loadMesaHistory(mesa_object, age_indices):
    labels = ["omega", "v_rot", "j_tot"]
    labels += ["period_days", "binary_separation", "eccentricity", "v_orb_1", "v_orb_2", "rl_1", "rl_2", "log_abs_mdot", "lg_mstar_dot_1", "lg_mstar_dot_2"]#, "star_1_radius", "star_2_radius"]
    labels += ["star_1_mass", "star_2_mass"]  # needed for the Roche-lobe mass ratio (roche_lobe.py)

    sm = mesa_object.hist.star_mass
    time_indices = len(sm)
    data_t = {}
    logTeff_values = []
    logTeff_colors = []

    # print("Keys history: ", mesa_object.hist.keys())

    # print("History keys available: ", mesa_object.hist.keys())
    for lab in labels:
        if lab in mesa_object.hist.keys():
            print(f"  Loading {lab}")
            _d = []
            for t in age_indices:
                _d.append(mesa_object.hist[lab][t])
            data_t.update({lab: np.array(_d)})
        else:
            print(f"  X did not find {lab}")

    if all(k in data_t for k in ("binary_separation","v_orb_1", "v_orb_2")):
        # Using r1/r2 = v1/v2 and r1+r2=r
        v1 = data_t["v_orb_1"]
        v2 = data_t["v_orb_2"]
        r = data_t["binary_separation"]
        r1 = v1/(v1+v2) * r
        r2 = v2/(v2+v1) * r

        data_t.update({"r1":r1, "r2":r2})

        print("  Added r1 and r2")

    return data_t






def rescale_time(indices, m, time_scale_type="model_number"):
    """Rescale the time.
    
    Rescale time indices depending on the time_type.
    
    Parameters
    ----------    
    indices : np.array or list of int
        Containing selected indices.
    m : mesa Object
    time_scale_type : str
        One of `model_number`, `linear`, or `log_to_end`. For `model_number`, the time follows the moment when a new MESA model was saved. For `linear`, the time follows linear steps in star_age. For `log_to_end`, the time axis is tau = log10(t_final - t), where t_final is the final star_age of the model.
    
    Returns
    -------
    ind_select : list
        New list of indices that reflect the rescaling in time.
    """

    def find_closest(ary, value):
        return int(np.abs(ary - value).argmin())

    age = m.hist.star_age
    if time_scale_type == "model_number":
        return [int(i) for i in np.array(indices)]
    elif time_scale_type == "linear":
        val_select = np.linspace(age[indices[0]], age[indices[-1]], len(indices))
        ind_select = [find_closest(val, age) for val in val_select]
        return ind_select
    elif time_scale_type == "log_to_end":
        time_diff = (age[-1] - age)
        # Avoid invalid values for log
        time_diff[time_diff <= 0] = 1e-5
        logtime = np.log10(time_diff)
        # Find indices
        val_select = np.linspace(logtime[indices[0]], logtime[indices[-1]], len(indices))
        ind_select = [find_closest(val, logtime) for val in val_select]
        return ind_select
    else:
        raise ValueError('Invalid time_type. Choose one of "model_number", "linear", or "log_to_end"')


