from  oceantracker.velocity_modifiers._base_velocity_modifer import VelocityModiferBase
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from numba import njit
from oceantracker.util.numba_util import njitOT

import numpy as np
from datetime import datetime,timezone,timedelta
import astral.sun

class TerminalVelocity(VelocityModiferBase):
    # add terminal velocity to particle velocity  < 0 is downwards ie sinking

    def __init__(self,):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'value': PVC(0.,float, doc_str='Terminal velocity positive upwards, ie fall velocities ate negative'),
                                 'mean': PVC(0., float, obsolete='use "value" parameter'),
                                 'variance': PVC(None, float, min=0., doc_str='variance of normal distribution of terminal velocity, used to give each particles its own terminal velocity from random normal distribution'),
                                 })


    def check_requirements(self):
        self.check_class_required_fields_prop_etc(requires3D=True, required_props_list=['velocity_modifier'])


    def initial_setup(self):
        super().initial_setup()
        si = self.shared_info
        particle= si.classes['particle_group_manager']

        si.msg_logger.msg('When using a terminal velocity, ensure time step is small enough that vertical displacement is a small fraction of the water depth, ie vertical Courant number < 1',note=True)

        if self.params['variance'] is not None:
           # set up individual particle terminal velocties
           particle.add_particle_property('terminal_velocity','user',dict(
                                            class_name='oceantracker.particle_properties.particle_parameter_from_normal_distribution.ParticleParameterFromNormalDistribution',
                                                    value=self.params['value'], variance=self.params['variance']))

    def update(self, time_sec, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        part_prop = si.classes['particle_properties']
        velocity_modifier = part_prop['velocity_modifier']

        if self.params['variance'] is None:
            # constant fall vel
            self._add_constant_vertical_vel(velocity_modifier.data, self.params['value'] * si.model_direction, active)
        else:
            self._add_individual_vertical_vel(velocity_modifier.data, part_prop['terminal_velocity'].data,  si.model_direction, active)

    @staticmethod
    @njit
    def _add_constant_vertical_vel(v, w, sel):
        for n in sel:
            v[n, 2] += w

    @staticmethod
    @njit
    def _add_individual_vertical_vel(v, w, model_dir, sel):
        for n in sel:
            v[n, 2] += w[n]*model_dir


class TurbidityInducedSinking(TerminalVelocity):
    """
    We assume that a particle drifts through the water neutrally buoyent.
    Based on the local turbidity a sinking velocity is induced that changes 
    over time. This is a parameterisation of a collision based model."""

    def __init__(self):
        super().__init__()
        self.add_default_params({'turbidity_field': PVC('turbidity',str),
                                 'min_sinking_velocity': PVC(-0.01,float),
                                 'max_sinking_velocity': PVC(-0.1,float),
                                 'min_turbidity': PVC(0.1,float),
                                 'max_turbidity': PVC(1.,float)
                                 })
        
    def initial_setup(self):
        pass    

class TurbidityInducedPrecipitationLikeSinking(TerminalVelocity):
    """
    We assume that a particle drifts through the water neutrally buoyent.
    Based on turbidity it has a chance to collide with a dust particle.
    After the collisions it sinks with a constant velocity.
    Collision chances are calculated based on assumed particle sizes and 
    velocities.
    """
    def __init__(self):
        super().__init__()
        self.add_default_params({'turbidity_field': PVC('turbidity',str),
                                 'diameter_of_particles': PVC(1e-6,float),
                                 'diameter_of_dust': PVC(1e-5,float),
                                 'density_dust': PVC(200.0,float),
                                 'max_sinking_velocity': PVC(-0.01,float)
                                 })
        
    def initial_setup(self):
        super().initial_setup()
        si = self.shared_info
        particle= si.classes['particle_group_manager']

        si.msg_logger.msg('When using a terminal velocity, ensure time step is small enough that vertical displacement is a small fraction of the water depth, ie vertical Courant number < 1',warning=True)

        # set up individual particle terminal velocties
        particle.create_particle_property('terminal_velocity',
                                          'user',
                                           dict(class_name='oceantracker.particle_properties.particle_parameter_from_uniform_distribution.ParticleParameterFromUniformDistribution',
                                                value=0)
        )

        # collision cross section
        self.info['sigma'] = np.pi * (self.params['diameter_of_particles']/2 + self.params['diameter_of_dust']/2)**2

        # particles per kg/l of dust
        volume_of_ball = 4/3 * np.pi * (self.params['diameter_of_particles']/2)**3
        self.info['particles_per_concentration'] = (volume_of_ball * self.params['density_dust'])**-1

        # fix?
        self.info['average_relative_velocity'] = self.shared_info.classes['dispersion'].params['A_H']/2

    def update(self, time_sec, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        part_prop = si.classes['particle_properties']
        velocity_modifier = part_prop['velocity_modifier']

        particle_per_m3 = part_prop[self.params['turbidity_field']].data[active] * self.info['particles_per_concentration']

        cross_section = self.info['sigma']
        avg_velocity = self.info['average_relative_velocity']
        collision_frequency = particle_per_m3 * cross_section * avg_velocity
        
        collision_chance = 1 - np.exp(-collision_frequency*self.shared_info.solver_info['model_time_step'])


        rolling_for_collision = np.random.rand(len(active)) < collision_chance
        part_prop['terminal_velocity'].data[active[rolling_for_collision]] = self.params['max_sinking_velocity']
        

        # # randomly select if particle collides
        # self._roll_die_for_collision(velocity_modifier.data, 
        #                              active, 
        #                              part_prop[self.params['turbidity_field']].data,
        #                              cross_section,
        #                              avg_velocity,
        #                              self.info['particles_per_concentration'],
        #                              self.shared_info.solver_info['model_time_step'],
        #                              self.params['max_sinking_velocity'] * si.model_direction)

        self._add_individual_vertical_vel(velocity_modifier.data, part_prop['terminal_velocity'].data,  si.model_direction, active)




    @staticmethod
    # @njit
    def _roll_die_for_collision(v, active, turbidity, cross_section, avg_velocity ,particles_per_concentration, model_time_step, max_sinking_velocity):
        for n in active:

            # calculate collision chances
            particle_per_m3 = turbidity[n] * particles_per_concentration
            collision_frequency = particle_per_m3 * cross_section * avg_velocity
            # collision_chance = 1 - np.exp(-collision_frequency*model_time_step)
            collision_chance = collision_frequency**model_time_step

            if np.random.rand() < collision_chance:
                #set terminal velocity to max sinking velocity
                
                v[n, 2] = max_sinking_velocity   

    @staticmethod
    @njit
    def _add_constant_vertical_vel(v, w, sel):
        for n in sel:
            v[n, 2] += w

    @staticmethod
    @njit
    def _add_individual_vertical_vel(v, w, model_dir, sel):
        for n in sel:
            v[n, 2] += w[n]*model_dir






class DielVelocity(TerminalVelocity):
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
        model_start_time = datetime.fromtimestamp(self.shared_info.solver_info['model_start_time'])
        self.noon = astral.sun.sun(location.observer, date=model_start_time)['noon'].timestamp()


    def calculate_sun_phase(self,t):
        x = np.abs(t - self.noon + self.params['phase'])%self.params['period']
        x = x/self.params['period']
        v_velocity_scaling = np.cos(2*np.pi*x)
        v_velocity_scaling = np.sign(v_velocity_scaling)
        return v_velocity_scaling


    def update(self, nb, t, active):
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        v_scaling = self.calculate_sun_phase(t)
        velocity_modifier = part_prop['velocity_modifier']       

        if self.params['variance'] is None:
            self._add_constant_vertical_vel(velocity_modifier.data, self.params['mean'] * si.model_direction * v_scaling, active)
        else:
            self._add_individual_vertical_vel(velocity_modifier.data, part_prop['terminal_velocity'].data,  si.model_direction * v_scaling, active)
