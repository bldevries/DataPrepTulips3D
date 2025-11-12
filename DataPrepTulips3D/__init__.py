
import importlib.metadata
__version__ = importlib.metadata.version('DataPrepTulips3D')

from DataPrepTulips3D.mesa_data import loadMesaData
from DataPrepTulips3D.Data1D import Data1D, save_Data1D_to_pickle, load_Data1D_from_pickle

