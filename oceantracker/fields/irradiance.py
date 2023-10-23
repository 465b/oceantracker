from oceantracker.particle_properties._base_properties import ParticleProperty
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from oceantracker.fields._base_field import UserFieldBase
from oceantracker.util.parameter_checking import ParamValueChecker as PVC

from numba import njit
import numpy as np
from datetime import datetime
import pvlib

class Irradiance(UserFieldBase):
    """
    Calculates irradiance just below the water surface.
    Assumes that model location is small enogh 
    such that the sun is approx. at the same angle for the full model domain.

    This field is designed to be used by the particle propery 'illumination'
    """
    def __init__(self):
        super().__init__()
        self.add_default_params(
                {
                    'latitude': PVC(0.0, float),
                    'longitude': PVC(0.0, float),
                    'timezone': PVC('UTC', str),
                    'albedo': PVC(0.1, float),
                    'is_time_varying': PVC(True,bool),
                    'num_components': PVC(1, int),
                    'is3D': PVC(False,bool)
                }
            )
    
    def check_requirements(self):
        pass

        # self.check_class_required_fields_prop_etc(
        #     required_fields_list=['turbidity'],
        #     requires3D=True
        #     )
    
        
    def initial_setup(self):
        super().initial_setup()

        # preparing irradiance data
        self.location = pvlib.location.Location(self.params['latitude'],
                                                          self.params['longitude'],
                                                          tz=self.params['timezone'])
        
    def initial_value_at_birth(self, active):
        # initial age is zero
        self.set_values(np.nan, active)





    def update(self, buffer_index):
        si = self.shared_info
        grid = si.classes['reader'].grid
        fields = si.classes['fields']
        self.calc_illumination(buffer_index,grid,self.location,self.params['albedo'],self.data)
        # self.calc_fiction_velocity(buffer_index, grid['zlevel'], grid['bottom_cell_index'], si.z0, fields['water_velocity'].data , self.data)


    @staticmethod
    # @njit()
    def calc_illumination(buffer_index,grid,location,albedo,out):
        for nt in buffer_index:
            # time in posix
            current_time = grid['time'][nt]
            # time in datetime
            current_time = datetime.utcfromtimestamp(current_time)
    
            # get Global horizontal irradiance (W/m^2) assuming clear sky
            ghi = location.get_clearsky(current_time, linke_turbidity=3)['ghi'][0]
            # account for water reflectivity
            subsurface_irradiance = (1-albedo)*ghi
            
            # set values for all grid nodes
            out[nt] = subsurface_irradiance