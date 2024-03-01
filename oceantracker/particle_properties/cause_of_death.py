from oceantracker.particle_properties._base_particle_properties import ParticleProperty
import numpy as np
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from oceantracker.particle_properties.util import particle_operations_util


class CauseOfDeath(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(0., float,doc_str='Particle property at the time of release')})
        self.class_doc(description='Particle cause of death, manually updated by culling trajectory modifiers.')


    def initial_setup(self):
        super().initial_setup()

        # self.flags = {
        #     'salinity': 1,
        #     'illumination': 2,
        #     'dryout': 3
        #     }

        # self.culled_by = {
        #     "salinity": set(),
        #     "illumination": set(),
        #     "dryout": set()
        # }


    def initial_value_at_birth(self, new_part_IDs):

        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values


    def update(self,active):
        # manually updated by culling classes
        pass
        
        # for item in self.culled_by:
        #     particle_operations_util.set_value(self.data,self.flags[item],
        #     np.array(list(self.culled_by[item]),dtype=int))
        

