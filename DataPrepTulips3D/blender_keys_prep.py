
# OBJECT KEYS TO TRACK
key_ob_DataPrepTulips3D_data = 'tulips_grid_data' # This holds the DataPrepTulips3D dict
key_ob_active_data_label = 'tulips_active_profile_key' # Label of the currently shown data 
key_ob_active_time_index = 'active_time_index' # Currently shown time index

key_blender_ob_data_dict = "MESA_data_dict"


prof_labels = "prof_labels"
chem_elem_labels = "chem_elem_labels"
r_resolution = "r_resolution"
t_resolution = "t_resolution"
t_resolution_orig = "t_resolution_orig" # Number of timesteps in the original (un-resampled) MESA history
time_scale_type = "time_scale_type" # How the resampled time steps were spaced ("model_number", "linear", or "log_to_end")

key_ob_cumm_abun_label = "abundances_cummulative"
nr_theta_points = "nr_Th"

# The data_array (given by DP.load_from_pickle) will contain 
# the "data_prof_t_r" key. It contains the MESA data in the 
# resolution (prof_labels, t_resolution, r_resolution)
data_prof_t_r = "data_prof_t_r"
data_chem_t_r = "data_chem_t_r"
data_chem_abun_id = "data_chem_abun"
data_t = "data_t"

data_t_Teff = "logTeff"

data_r_max = "data_r_max"
data_r_max_Rsun = "data_r_max_r_sun"

# Roche lobe equipotential-surface radius grids (theta x phi, per frame), one
# per star in a binary - computed from the mass ratio via roche_lobe.py.
# Kept as their own top-level dict keys (like data_chem_abun_id above), not
# inside data_t, since they're 3D (frame x theta x phi) arrays and data_t's
# generic per-key baking in save_to_texture assumes a plain 1D
# (one-scalar-per-frame) array.
data_rochelobe_shape_1 = "data_rochelobe_shape_1"
data_rochelobe_shape_2 = "data_rochelobe_shape_2"
nr_theta_points_rochelobe = "nr_Th_RL"
nr_phi_points_rochelobe = "nr_Phi_RL"

# Full Roche potential "landscape" (a height-field grid over the whole
# orbital plane, not just the single Roche-lobe equipotential contour above)
# - a prototype 3D surface-plot visualization of the potential itself. One
# shared grid per binary (not per-star), (frame x n_y x n_x). Same reasoning
# as the Roche lobe shape grids for living outside data_t.
data_rochelobe_potential_grid = "data_rochelobe_potential_grid"
nr_x_points_potential = "nr_x_potential"
nr_y_points_potential = "nr_y_potential"
potential_x_range = "potential_x_range"  # (x_min, x_max), dimensionless (separation units), fixed per run
potential_y_range = "potential_y_range"  # (y_min, y_max), dimensionless (separation units), fixed per run

# L2/L3 "outer" equipotential rings (closed contours just inside each
# critical value, enclosing both stars) - see roche_lobe.equipotential_boundary_xy.
# Unlike the L1 lobe shape grids (a scalar radius per theta, from each
# star's own center), these store full (x,y) per ring point - the merged
# shape isn't star-shaped from any single fixed center for realistic mass
# ratios - so each is (n_frames, n_theta, 2), Rsun.
data_rochelobe_l2_ring = "data_rochelobe_l2_ring"
data_rochelobe_l3_ring = "data_rochelobe_l3_ring"
nr_theta_points_equipotential = "nr_Th_equipot"

dir_structure = "dir_structure"
texture_dir = "texture_dir"
data_t_r_filename = "dir_structure_data_t_r_filename"
chem_abun_filename = "dir_structure_chem_abun_filename"
data_t_filename_and_max = "dir_structure_data_t_filename_and_max"
rochelobe_shape_filename_1 = "dir_structure_rochelobe_shape_filename_1"
rochelobe_shape_filename_2 = "dir_structure_rochelobe_shape_filename_2"
rochelobe_potential_filename = "dir_structure_rochelobe_potential_filename"
rochelobe_l2_ring_filename = "dir_structure_rochelobe_l2_ring_filename"
rochelobe_l3_ring_filename = "dir_structure_rochelobe_l3_ring_filename"

is_binary = "is_binary"
binary_nr = "binary_nr"