# import bpy
import numpy as np
import sys, os
# from time import time
# import matplotlib.pyplot as plt
# from matplotlib.colors import Normalize, LogNorm
# from scipy.interpolate import interp1d

import mesaPlot as mp

# import Data1D.Data1D

import DataPrepTulips3D as DP

def loadMesaData(mesa_LOGS_directory, filename_history = None, verbose_timing = False, profiles=['mass', 'logT', 'logRho', 'he4']):
    """
    Loads MESA data into a data1D object

    Parameters:
    mesa_LOGS_directory (str): directory of the LOGS file
    publication_year (int): The year the book was published.
    filename_history (str): optional, default=history.data
    
    Returns:
    Data1D: data1d object containing mesa data
    """

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
        
    if verbose_timing: print("Timing load mesa file: ", time()-_)

    # Create the Data1D object. Use the file name as identifier to
    # figure out later what this data file is about
    data = DP.Data1D(time_grid=m.hist.star_age, identifier=f)

    # Load the Energy generated into the data object
    loadMesaEnergyData(m, data)
    # Load the Teff data of the star
    loadMesaTeffData(m, data)

    # Load the radial grid a.f.o. time    
    property_name_r_axis = 'logR'
    loadMesaProfile(m , data, mesa_LOGS_directory, property_name=property_name_r_axis)
        
    # Load the different profiles
    property_name_values = profiles #, 'zone', 'logP', 'h1', 'he3', ]#, , 'energy', 'x_mass_fraction_H', 'y_mass_fraction_He', 'z_mass_fraction_metals',  'c12', 'n14', 'o16', 'ne20', 'mg24', 'si28', 's32', 'ar36', 'ca40', 'ti44', 'cr48', 'fe52', 'ni56', 'fe54', 'fe56', 'cr56', 'opacity', 'luminosity', 'mlt_mixing_length', 'mlt_mixing_type', 'conv_vel', 'mixing_type', 'log_D_mix', 'log_D_mix_non_rotation', 'log_D_conv', 'log_D_semi', 'log_D_ovr', 'log_D_thrm', 'tau', 'omega', 'j_rot', 'fp_rot', 'ft_rot', 'r_polar', 'r_equatorial', 'am_log_D_visc', 'am_log_D_DSI', 'am_log_D_SH', 'am_log_D_SSI', 'am_log_D_ES', 'am_log_D_GSF', 'am_log_D_ST', 'am_log_nu_ST', 'dynamo_log_B_r', 'dynamo_log_B_phi']
    for p in property_name_values:
        loadMesaProfile(m , data, mesa_LOGS_directory, property_name=p, r_grid = data.GridPropertiesList["logR"].value_afo_time)
    return data


def loadMesaProfile(m, data, mesa_LOGS_directory, property_name="logRho", r_grid=None):#, raxis="mass"):
    '''Reads in a profile data from a mesa file'''

    def find_profile(m, mesa_LOGS_directory, time_ind=0):
        model_number = m.hist.model_number[time_ind]
        m.loadProfile(num=model_number, f=mesa_LOGS_directory, silent=True)
        return m.prof
    
    # prof = find_profile(m, mesa_LOGS_directory, time_ind=0)
    # print("sdfsdf", prof.data.dtype.names)
    
    # list_r = []
    list_prop = []
    for i, age in enumerate(m.hist.star_age):
        prof = find_profile(m, mesa_LOGS_directory, time_ind=i)
        # r = prof.data[raxis]
        prop = prof.data[property_name][:]
        # print(i, age, r.shape, prop.shape)
        # if r.shape != prop.shape:
            # raise ValueError(f"Radius and property array length do not match for model number {m.hist.model_number[i]}")

        # list_r.append(r)
        list_prop.append(prop)

    if r_grid:
        data.set_grid_property(property_name, list_prop, r_afo_time = r_grid)#, list_r)
    else:
        data.set_grid_property(property_name, list_prop, r_afo_time = list_prop)#, list_r)

def loadMesaTeffData(mesa_object , data, verbose_timing=False):
    '''Reads in Teff data from a mesa file'''

    sm = mesa_object.hist.star_mass
    time_indices = len(sm)
    logTeff_values = []

    start_ind = 0
    time_index_step=1
    while start_ind < time_indices:
        logTeff_values.append(mesa_object.hist.log_Teff[start_ind])
        start_ind += time_index_step

    data.set_total_property("logTeff", logTeff_values)



def loadMesaEnergyData(mesa_object, data, verbose_timing=False):
    '''Reads in Energy production data from a mesa file'''

    # We need these indices for the burning data
    qtop = "burn_qtop_"
    qtype = "burn_type_"

    # A check if data is avaliable
    try:
        mesa_object.hist.data[qtop + "1"]
    except ValueError:
        raise KeyError(
            "No field " + qtop + "* found, add mixing_regions 40 and burning_regions 40 to your history_columns.list")

    sm = mesa_object.hist.star_mass

    time_indices = len(sm)

    if verbose_timing: _ = time()
    R_star = []
    t_index = []
    E_value = []
    E_grid = []

    start_ind = 0
    time_index_step=1
    while start_ind < time_indices:
        # print(start_ind)
        num_burn_zones = int([xx.split('_')[2] for xx in mesa_object.hist.data.dtype.names if qtop in xx][-1])
        # Per time stamp we will have radii and values, which we list in these variables:
        list_r = [np.abs(mesa_object.hist.data[qtop + str(region)][start_ind] * sm[start_ind]) for region in range(1, num_burn_zones + 1)]
        list_E = [mesa_object.hist.data[qtype + str(region)][start_ind] for region in range(1, num_burn_zones + 1)]

        # We make one data array
        _d = np.array([list_r, list_E])
        # We need to remove duplicates at the end where value==-9999
        _mask = _d[1] != -9999
        # Our cleaned up arrays containing the radii and burning values
        r, v = _d[0][_mask], _d[1][_mask] # Now you have data you can interpolate f = interp1d(r, v, kind='cubic')
        E_value.append(v)
        E_grid.append(r)
        # We save the times where we sample and the stellar radius at that time        
        t_index.append(start_ind)
        R_star.append(r[-1]) 
        # And we make vertex colors
        #tulips3dGeometry.make_vertex_colors(r, v, settings, vertex_colors_name_base="energy_ver_col_"+str(start_ind))
        # And increment the time index
        start_ind += time_index_step

    data.set_total_property("stellar_mass", sm)
    data.set_total_property("stellar_radius", R_star)
    data.set_grid_property("energy", E_value, E_grid)
