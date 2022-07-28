from  oceantracker.velocity_modifiers._base_velocity_modifer import VelocityModiferBase
from oceantracker.particle_properties.util import particle_operations_util
from oceantracker.particle_properties.particle_parameter_from_normal_distribution import  ParticleParameterFromNormalDistribution
from oceantracker.util.parameter_checking import ParamDictValueChecker as PVC

import numpy as np
from datetime import datetime,timezone,timedelta
import astral.sun

class AddTerminalVelocity(VelocityModiferBase):
    # add terminal velocity to particle velocity  < 0 is downwards ie sinking

    def __init__(self,):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'name': PVC('terminal_velocity',str),
                                 'mean': PVC(0.,float),
                                 'variance': PVC(0.,float, min=0.)})

        # only possible in in 3D so tweak flag

    def check_requirements(self):
        msg_list = self.check_class_required_fields_properties_grid_vars_and_3D(requires3D=True)
        return msg_list

    def initialize(self):
        super().initialize()
        particle= self.shared_info.classes['particle_group_manager']

        if self.params['variance'] > 0.:
           # set up individual particle terminal velocties
            p = ParticleParameterFromNormalDistribution()
            particle.create_particle_property('user',dict(name='terminal_velocity', instance=p,
                                             mean=self.params['mean'], variance=self.params['variance']))

    def modify_velocity(self,v, t, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        if self.params['variance'] == 0.:
            # constant fall vel
            particle_operations_util.add_value_to(v[:, 2], self.params['mean'] * si.model_direction, active)
        else:
            particle_operations_util.add_to(v[:, 2], si.classes['particle_properties']['terminal_velocity'].data, active, scale = si.model_direction)

class AddDielVelocity(AddTerminalVelocity):
    """
    Adds verticle velocities to the particles following a rectangular 
    wave pattern.
    Phase of the diel migration is calculated relative to the sun rise
    """
    def __init__(self):
        super().__init__()
        self.add_default_params({'period': PVC(24*60*60.,float,min=1),
                                 'phase': PVC(0.,float),
                                 'location': PVC([53.551086,9.993682],list)})

    def initialize(self):
        super().initialize()
        
        location = astral.LocationInfo(timezone="UTC",latitude=self.params['location'][0],
                                                      longitude=self.params['location'][1])
        model_start_time = datetime.fromtimestamp(self.shared_info.model_start_time)
        self.noon = astral.sun.sun(location.observer, date=model_start_time)['noon'].timestamp()


    def calculate_sun_phase(self,t):
        x = np.abs(t - self.noon + self.params['phase'])%self.params['period']
        x = x/self.params['period']
        v_velocity_scaling = np.cos(2*np.pi*x)
        v_velocity_scaling = np.sign(v_velocity_scaling)
        return v_velocity_scaling

    def modify_velocity(self,v, t, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        v_scaling = self.calculate_sun_phase(t)

        if self.params['variance'] == 0.:
            # constant fall vel
            particle_operations_util.add_value_to(v[:, 2],
                self.params['mean'] * si.model_direction * v_scaling, active)
        else:
            particle_operations_util.add_to(v[:, 2],
                si.classes['particle_properties']['terminal_velocity'].data, active,
                scale = si.model_direction)