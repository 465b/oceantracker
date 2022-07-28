import oceantracker.particle_velocity._velocity_modifer_base as baseModifier
import oceantracker.util.particle_operations as part_op
from  oceantracker.user_particle_properties.particleParameterFromNormalDistribution import  ParticleParameterFromNormalDistribution


class add_time_dependent_verticle_migration(baseModifier.VelocityModiferBase):
    # terminal velocity >0 is downwards ie sinking

    def __init__(self,):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'name': 'terminal_velocity','mean': None,'variance': None})

        # only possible in in 3D so tweak flag
        self.requires3D=True
        self.requirements_txt = 'Must be 3D hindcast'

    def initialize(self):
        super().initialize()
        particle= self.shared_info.pointers['particle']

        if self.params['variance'] is not None:
           # set up individual particle terminal velocties
            p = ParticleParameterFromNormalDistribution()
            particle.add_particle_property(name='terminal_velocity', type='parameter',
                                            initial_value=0., instance= p,
                                            mean=self.params['mean'],variance=  self.params['variance'])

            self.terminal_velocity = particle.get_partProp_dataPtr('terminal_velocity')

    def can_be_added(self, is3Dhindcast):
        return is3Dhindcast # only if 3D

    def time_dependent_vertical_migration(time,period=86400,phase_shift=3600):
        amplitude = np.sin((time+phase_shift)*2*np.pi/period)
        return amplitude

    def modify_velocity(self,v, t, active):
        # modify vertical velocity, if backwards, make negative
        if self.params['variance'] is None:
            # constant fall vel
            vertical_vel = self.params['mean']*self.shared_info.model_dir
            vertical_vel *= self.time_dependent_vertical_migration(t,self.params['period'],self.params['phase_shift'])
            part_op.add_value_to(v[:, 2], vertical_vel, active)
        else:
            part_op.add_to(v[:,2],self.terminal_velocity, active, scale = self.shared_info.model_dir)


