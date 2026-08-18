import DataPrepTulips3D as DP
import numpy as np
import matplotlib.pyplot as plt

# d_log = DP.loadMesaData(mesa_LOGS_directory = "../../../example_MESA_data/single_11Msun/LOGS", \
#                     t_resolution = 500, r_resolution = 100,\
#                     time_scale_type="log_to_end",\
#                     filename_history = "history.data")
# d_log.keys(), 
# DP.save_to_pickle(d_log, "../../../example_MESA_data/DataDictFormat/", "single_11Msun_log_to_end")


d = DP.load_from_pickle("../../../example_MESA_data/DataDictFormat/"+"single_11Msun_log_to_end.pkl")
	#"../../../example_MESA_data/DataDictFormat/binary.pkl")


DP.save_to_texture(d, "test_textures")

