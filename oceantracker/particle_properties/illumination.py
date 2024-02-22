from oceantracker.particle_properties._base_particle_properties import ParticleProperty
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from oceantracker.particle_properties.util import particle_operations_util

from numba import njit
import numpy as np
import pvlib

class Illumination(ParticleProperty):
    def __init__(self):
        super().__init__()
        self.add_default_params(
                {
                    'name_of_turbidity_field': PVC('turbidity', str),
                    'name_of_irradiance_field': PVC('irradiance', str),
                    'c': PVC(0.15, float),
                    # 'is_time_varying': PVC(True,bool),
                    # 'num_components': PVC(1, int),
                    # 'is3D': PVC(False,bool)
                }
            )
        
    def check_requirements(self):

        self.check_class_required_fields_prop_etc(
            # required_fields_list=[
            #     'turbidity', 
            #     'irradiance', 
            #     'tide'
            # ],
            # required_props_list=[
            #     'x'
            # ],
            requires3D=True
            )
    
        
    def initial_setup(self):
        super().initial_setup()

        # preparing irradiance data
        # self.info['location'] = pvlib.location.Location(self.params['latitude'],
        #                                                   self.params['longitude'],
        #                                                   tz=self.params['timezone'])
        
    def initial_value_at_birth(self, active):
        # initial age is zero
        self.set_values(0, active)





    def update(self, active):
        si = self.shared_info
        # grid = si.classes['reader'].grid
        
        part_prop = self.shared_info.classes['particle_properties']
        fields = si.classes['fields']

        illumination = self.calc_illumination(active, 
                               part_prop['x'].data,
                               part_prop['turbidity'].data,
                               part_prop['irradiance'].data,
                               part_prop['tide'].data,
                               self.params['c'],
                               self.data)
        
        # particle_operations_util.set_value(self.data, illumination, active)
        self.data[active] = illumination


    @staticmethod
    # @njit()
    def calc_illumination(active,x,turbidity,irradiance,tide,c,out):
        z = x[active,2]
        tide = tide[active]

        depth = -(z-tide)
        np.clip(depth, 0, None)

        illumination = irradiance[active] * np.exp(- c * turbidity[active] * depth)

        return illumination
        
        
class AverageIllumination(Illumination):
    def __init__(self):
        super().__init__()
        self.add_default_params(
                {
                    "time_to_average": PVC(24*3600, float),
                }
            )

    def initial_setup(self):
        """ 
        To calculate the running average we need to get 
        the number of time steps taken during the averaging period.
        """
        super().initial_setup()

        self.info['time_steps_in_averaging_period'] = int(self.params['time_to_average'] / self.shared_info.settings['time_step'])


    def update(self, active):
        si = self.shared_info
        # grid = si.classes['reader'].grid
        
        particle = self.shared_info.classes['particle_properties']
        
        z = particle['x'].data[active,2]
        tide = particle['tide'].data[active]

        depth = -(z-tide)
        np.clip(depth, 0, None)

        current_illumination = particle['irradiance'].data[active] * np.exp(- self.params['c'] * particle['turbidity'].data[active] * depth)

        # calculate the running average
        recent_illumination = self.data[active]

        averaged_illumination = (recent_illumination * (self.info['time_steps_in_averaging_period']-1) + current_illumination) / self.info['time_steps_in_averaging_period']
        
        # particle_operations_util.set_value(self.data, illumination, active)
        self.data[active] = averaged_illumination