
from oceantracker.util.parameter_checking import  ParamValueChecker as PVC
from oceantracker.particle_properties._base_properties import ParticleProperty

import numpy as np

class light_limitation(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({'name': PVC('light_limitation', str) ,
                                 'description': PVC( 'culling light starved plankton', str),
                                 'initial_value': PVC(None, float),
                                 'max_time_wo_light': PVC( 3600*24*10., float)})
        # get time step size (needed as a multiplier)
        #self.dt = self.shared_info.model_substep_timestep

    def initial_value_at_birth(self, new_part_IDs):
        if self.params['initial_value'] is None:
            s= new_part_IDs.shape[0]
            self.set_values(self.params['max_time_wo_light']*np.random.rand(s), new_part_IDs)
        else:
            self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # update decay prop each time step
        particle_properties = self.shared_info.classes['particle_properties']
        particle_status = self.shared_info.classes['particle_properties']['status'].data[active]
        dead_flag =self.shared_info.particle_status_flags['dead']
        
        particle_z = particle_properties['x'].get_values(active)[:,2]
        water_level = particle_properties['tide'].get_values(active)
        particle_depth = water_level - particle_z

        suffocating = np.where(particle_depth > 1)
        suffocating = active[suffocating]
        # yields the elements in first list that are NOT in second
        not_suffocating = np.setdiff1d(active,suffocating)
        
        self.add_values_to(self.shared_info.solver_info['model_timestep'], suffocating)
        self.add_values_to(-self.shared_info.solver_info['model_timestep'], not_suffocating)
        # only non-zero positive elements should be subtracted by dt
        self.set_values(0, active[np.where(self.data[active] < 0)])

        # set particles that are above the threshold to dead
        light_starved = np.where(self.data > self.params['max_time_wo_light'])
        
        particle_properties['status'].set_values(dead_flag, light_starved)
