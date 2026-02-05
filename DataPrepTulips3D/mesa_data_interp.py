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

# import Data1D.Data1D

import DataPrepTulips3D as DP

key_DataPrepTulips3D_prof_labels = "prof_labels"
key_DataPrepTulips3D_r_resolution = "r_resolution"
key_DataPrepTulips3D_t_resolution = "t_resolution"
key_DataPrepTulips3D_data_prof_t_r = "data_prof_t_r"
key_DataPrepTulips3D_data_r_max = "data_r_max"

def save_to_pickle(data_dict, filepath):
    '''Saves a dict into a pickle file'''
    dbfile = open(filepath, 'ab')
    pickle.dump(data_dict, dbfile)   
    dbfile.close()

def load_from_pickle(filepath):
    '''Reads a dict from a pickle file'''
    dbfile = open(filepath, 'rb')    
    data_dict = pickle.load(dbfile)
    return data_dict

def loadMesaData(mesa_LOGS_directory, t_resolution, r_resolution,\
                t_range = -1,\
                filename_history = None, verbose_timing = False, \
                profiles=['mass', 'logT', 'logRho', 'he4'], \
                r_grid_name="mass"):
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
    t_resolution (int): 
    r_resolution (int):
    filename_history (str):
    profiles (list of str):
    r_grid_name (str):
    verbose_timing (bool):

    Returns:
    mesa_data: dict containing the MESA data
    """
    # Besides the given profiles we also want to read in the energy production and Teff
    profiles = profiles + ["en"]

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

    # We want to reduce the resolution of the t_grid to t_resolution
    # Get t_resolution elements from t_grid using these indices:
    _t_grid = m.hist.star_age

    if t_range == -1: _t_indices = np.round(np.linspace(0, len(_t_grid) - 1, t_resolution)).astype(int)
    else: _t_indices = np.round(np.linspace(t_range[0], t_range[1] - 1, t_resolution)).astype(int)
    new_t_grid = _t_grid[_t_indices]

    if verbose_timing: print("Timing load mesa file: ", time()-_)

    print("Original/new time resolution: ", len(_t_grid), len(new_t_grid))

    # Load the Teff data of the star
    logTeff_values = loadMesaTeffData(m, new_t_grid)

    data_array_prof, R_star_from_grid = loadMesaProfile(m , mesa_LOGS_directory, \
                                                    r_resolution = r_resolution,\
                                                    profile_names=profiles, \
                                                    r_grid_name=r_grid_name, t_grid = new_t_grid)
    
    return  {\
            "info":"",\
            "MESA_file": mesa_LOGS_directory,\
            "filename_history":filename_history,\
            "r_label": r_grid_name,\
            key_DataPrepTulips3D_prof_labels: profiles,\
            key_DataPrepTulips3D_t_resolution: t_resolution, \
            key_DataPrepTulips3D_r_resolution: r_resolution,\
            key_DataPrepTulips3D_data_prof_t_r: data_array_prof,\
            key_DataPrepTulips3D_data_r_max: R_star_from_grid,\
            "data_t": {"logTeff": logTeff_values}\
            }


def loadMesaProfile(m, mesa_LOGS_directory, profile_names, \
                    r_resolution, 
                    r_grid_name, t_grid):
    
    """
    Reads in a profile data from a mesa file

    Parameters:
    
    Returns:
    
    """


    def find_profile(m, mesa_LOGS_directory, time_ind=0):
        model_number = m.hist.model_number[time_ind]
        m.loadProfile(num=model_number, f=mesa_LOGS_directory, silent=True)
        return m.prof
    
    _prof = find_profile(m, mesa_LOGS_directory)
    r_grid = _prof.data[r_grid_name][:]





    # The two output arrays.
    # This will contain the MESA data in the new resolutions
    data_array = np.zeros((len(profile_names), t_grid.shape[0], r_resolution))
    # This will contain the max Radius at each time index
    R_star_from_grid = np.zeros((len(profile_names), t_grid.shape[0]))

    for i_prof_name, pname in enumerate(profile_names):
        for t, age in enumerate(t_grid):
            if pname == "en":
                _R_max, _prop = loadMesaEnergyData(m, t, r_resolution)
            else:
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

            R_star_from_grid[i_prof_name,t] = _R_max
            data_array[i_prof_name,t,:] = _prop

    return data_array, R_star_from_grid



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

    #print("EN:", len(list_r), len(list_E), list_r[0], list_r[-1])
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


def loadMesaTeffData(mesa_object, t_grid, verbose_timing=False):
    '''Reads in Teff data from a mesa file'''

    sm = mesa_object.hist.star_mass
    time_indices = len(sm)
    logTeff_values = []

    for t, age in enumerate(t_grid):
        logTeff_values.append(mesa_object.hist.log_Teff[t])

    #data.set_total_property("logTeff", logTeff_values)
    return logTeff_values

# def loadMesaTeffData(mesa_object, t_grid, verbose_timing=False):
#     '''Reads in Teff data from a mesa file'''

#     sm = mesa_object.hist.star_mass
#     time_indices = len(sm)
#     logTeff_values = []

#     start_ind = 0
#     time_index_step=1
#     while start_ind < time_indices:
#         logTeff_values.append(mesa_object.hist.log_Teff[start_ind])
#         start_ind += time_index_step

#     #data.set_total_property("logTeff", logTeff_values)
#     return logTeff_values

# def loadMesaEnergyData_old(mesa_object, verbose_timing=False):
#     '''Reads in Energy production data from a mesa file'''

#     # We need these indices for the burning data
#     qtop = "burn_qtop_"
#     qtype = "burn_type_"

#     # A check if data is avaliable
#     try:
#         mesa_object.hist.data[qtop + "1"]
#     except ValueError:
#         raise KeyError(
#             "No field " + qtop + "* found, add mixing_regions 40 and burning_regions 40 to your history_columns.list")

#     sm = mesa_object.hist.star_mass

#     time_indices = len(sm)

#     if verbose_timing: _ = time()
#     R_star = []
#     t_index = []
#     E_value = []
#     E_grid = []

#     start_ind = 0
#     time_index_step=1
#     while start_ind < time_indices:
#         # print(start_ind)
#         num_burn_zones = int([xx.split('_')[2] for xx in mesa_object.hist.data.dtype.names if qtop in xx][-1])
#         # Per time stamp we will have radii and values, which we list in these variables:
#         list_r = [np.abs(mesa_object.hist.data[qtop + str(region)][start_ind] * sm[start_ind]) for region in range(1, num_burn_zones + 1)]
#         list_E = [mesa_object.hist.data[qtype + str(region)][start_ind] for region in range(1, num_burn_zones + 1)]

#         print("EN:", len(list_r), len(list_E), list_r[0], list_r[-1])
#         # We make one data array
#         _d = np.array([list_r, list_E])
#         # We need to remove duplicates at the end where value==-9999
#         _mask = _d[1] != -9999
#         # Our cleaned up arrays containing the radii and burning values
#         r, v = _d[0][_mask], _d[1][_mask] # Now you have data you can interpolate f = interp1d(r, v, kind='cubic')
#         E_value.append(v)
#         E_grid.append(r)
#         # We save the times where we sample and the stellar radius at that time        
#         t_index.append(start_ind)
#         R_star.append(r[-1]) 
#         # And we make vertex colors
#         #tulips3dGeometry.make_vertex_colors(r, v, settings, vertex_colors_name_base="energy_ver_col_"+str(start_ind))
#         # And increment the time index
#         start_ind += time_index_step

#     # data.set_total_property("stellar_mass", sm)
#     # data.set_total_property("stellar_radius", R_star)
#     # data.set_grid_property("energy", E_value, E_grid)
#     return E_value, E_grid
