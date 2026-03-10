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

# import DataPrepTulips3D as DP

key_DataPrepTulips3D_prof_labels = "prof_labels"
key_DataPrepTulips3D_r_resolution = "r_resolution"
key_DataPrepTulips3D_t_resolution = "t_resolution"
key_DataPrepTulips3D_data_prof_t_r = "data_prof_t_r"
key_DataPrepTulips3D_data_r_max = "data_r_max"

def save_to_pickle(data_dict, filepath):
    '''Saves a dict into a pickle file'''
    # dbfile = open(filepath, 'ab')
    dbfile = open(filepath, 'wb')
    pickle.dump(data_dict, dbfile)   
    dbfile.close()

def load_from_pickle(filepath):
    '''Reads a dict from a pickle file'''
    dbfile = open(filepath, 'rb')    
    data_dict = pickle.load(dbfile)
    return data_dict

def loadMesaData(mesa_LOGS_directory, t_resolution, r_resolution,\
                time_scale_type="log_to_end",\
                filename_history = None, verbose_timing = False, \
                profiles=[], r_grid_name="mass"):
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

    print("Original/new time resolution: ", len(m.hist.star_age), len(new_age_grid))

    # Load the Teff data of the star
    logTeff_values = loadMesaTeffData(m, age_indices=_age_indices)


    data_array_prof, R_star_from_grid = loadMesaProfile(m , mesa_LOGS_directory, \
                                                    r_resolution = r_resolution,\
                                                    profile_names=profiles, \
                                                    r_grid_name=r_grid_name, \
                                                    age_indices = _age_indices)
    
    return  {\
            "info":"",\
            "MESA_file": mesa_LOGS_directory,\
            "filename_history":filename_history,\
            "r_label": r_grid_name,\
            "age": new_age_grid,\
            "age_indices": list(_age_indices),\
            key_DataPrepTulips3D_prof_labels: profiles,\
            key_DataPrepTulips3D_t_resolution: t_resolution, \
            key_DataPrepTulips3D_r_resolution: r_resolution,\
            key_DataPrepTulips3D_data_prof_t_r: data_array_prof,\
            key_DataPrepTulips3D_data_r_max: R_star_from_grid,\
            "data_t": {"logTeff": logTeff_values}\
            }

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

    # This will contain the MESA data in the new resolutions
    data_array = np.zeros((len(profile_names), len(age_indices), r_resolution))

    # This will contain the max Radius at each time index
    R_star_from_grid = np.zeros((len(profile_names), len(age_indices)))

    for i_prof_name, pname in enumerate(profile_names):
        for i, t in enumerate(age_indices):
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

            R_star_from_grid[i_prof_name,i] = _R_max
            data_array[i_prof_name,i,:] = _prop

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

    for t in age_indices:
        logTeff_values.append(mesa_object.hist.log_Teff[t])

    return logTeff_values


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
        return indices
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


