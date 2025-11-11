from scipy.interpolate import interp1d, CubicSpline
import pickle


def save_Data1D_to_pickle(data1d, filepath):
    '''Saves a Data1D object into a pickle file'''
    dbfile = open(filepath, 'ab')
    pickle.dump(data1d, dbfile)   
    dbfile.close()

def load_Data1D_from_pickle(filepath):
    '''Reads a Data1D object from a pickle file'''
    dbfile = open(filepath, 'rb')    
    d_recov = pickle.load(dbfile)
    return d_recov


        
class Data1D():
    ''' Data object that holds multiple properties. These can be a function of time only (TotalProperty) or of time and radius (GridProperty)'''
    def __init__(self, time_grid, identifier):
        self.identifier = identifier
        self.time_grid = time_grid
        self.TotalPropertiesList = {}
        self.GridPropertiesList = {}

        # self._original_filename = None

    class TotalProperty():
        """A property of the whole object a.f.o. time"""
        def __init__(self, name, value_afo_time):
            self.name = name
            self.value_afo_time = value_afo_time

    class GridProperty():
        """ A property that varies over the radial grid  and a.f.o. time"""
        def __init__(self, name, value_afo_time, r_afo_time=[]):
            self.name = name
            self.value_afo_time = value_afo_time
            self.r_afo_time = r_afo_time
            # print(name, len(self.r_afo_time), len(self.value_afo_time ))

        def getValue(self):
            return self.value_afo_time

        def getGrid(self):
            return self.r_afo_time

    def print_summary(self):
        print("--")
        print(f"** Content of Data1D object **")
        print(f" * id = {self.identifier}")
        print(f" * Time grid length: {len(self.time_grid)}")
        print()
        print(f"Properties of the total object: ")
        for i, (k, v) in enumerate(self.TotalPropertiesList.items()):
            print(f" - {k} (Tot): len={len(v.value_afo_time)}")
        print(f"Properties as a function of radius: ")
        for i, (k, v) in enumerate(self.GridPropertiesList.items()):
            print(f" - {k} (Grid): len={len(v.value_afo_time)}")
        # for k in self.TotalPropertiesList.keys():
            # print(f"{k}: {self.TotalPropertiesList[k].shape}")
        print("--")

    def set_grid_property(self, name, value_afo_time, r_afo_time=[]):
        '''Sets a property that is a function of radius, and time'''
        if len(value_afo_time) != len(self.time_grid):
            raise ValueError("Value not the same length as time grid") 
        if len(r_afo_time) != 0:
            if len(r_afo_time) != len(self.time_grid): 
                raise ValueError("Grid not the same length as time grid") 
            for i in range(len(value_afo_time)):
                if(len(value_afo_time[i]) != len(r_afo_time[i])): 
                    raise ValueError(f"Value entry {i} and grid entry not the same length") 

        self.GridPropertiesList.update(\
            {name: self.GridProperty(name, value_afo_time, r_afo_time)}\
            )


    def set_total_property(self, name, value_afo_time):
        '''Sets a property that is a function time'''
        self.TotalPropertiesList.update(\
            {name: self.TotalProperty(name, value_afo_time)}\
            )

    def getTime(self):
        return self.time_grid

    def get_grid_property(self, name):
        return self.GridPropertiesList[name]

    def get_radial_interpolated_property(self, name, time_index=None):
        '''Get an interpolated property that is a function of radius'''
        prop = self.GridPropertiesList[name]

        if len(prop.r_afo_time) == len(prop.value_afo_time):
            if time_index:
                value_interpol = interp1d(prop.r_afo_time[i], prop.value_afo_time[i], kind='cubic')
            else:
                value_interpol = [interp1d(prop.r_afo_time[i], prop.value_afo_time[i], kind='cubic') for i in range(len(prop.value_afo_time))]
        else:
            value_interpol = -1

        return value_interpol
