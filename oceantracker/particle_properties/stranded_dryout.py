
# from oceantracker.util.parameter_checking import  ParamDictValueChecker as PVC
from oceantracker.util.parameter_checking import  ParamValueChecker as PVC
from oceantracker.particle_properties.util import particle_operations_util
from oceantracker.particle_properties._base_properties import ParticleProperty
import numpy as np


class StrandedDryout(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({'name': PVC('stranded_dryout', str) ,
                                 'description': PVC( 'culling stranded dryed particles', str),
                                 'initial_value': PVC(0., float),
                                 'max_time_stranded': PVC( 3600*24*10., float)})
        # get time step size (needed as a multiplier)
        #self.dt = self.shared_info.model_substep_timestep


    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # update decay prop each time step
        particle_properties = self.shared_info.classes['particle_properties']
        particle_status = particle_properties['status'].data[active]
        stranded_flag = self.shared_info.particle_status_flags['stranded_by_tide']
        dead_flag =self.shared_info.particle_status_flags['dead']
       
        # compare all particles (at once) and add dt or reset them
        stranded_particles = np.where(particle_status == stranded_flag)[0]
        # stranded_particles = active[stranded_particles]

        dt = self.shared_info.settings['time_step']
        
        # add values to stranded particles
        particle_operations_util.add_value_to(self.data,dt,stranded_particles)
        
        # reset values of not stranded particles
        not_stranded_particles = np.where(particle_status != stranded_flag)[0]
        not_stranded_particles = active[not_stranded_particles]
        
        particle_operations_util.set_value(self.data,0,not_stranded_particles)

        # set particles that are above the threshold to dead
        dried_out = active[np.where(self.data[active] > self.params['max_time_stranded'])]
        particle_properties['status'].set_values(dead_flag, dried_out)
