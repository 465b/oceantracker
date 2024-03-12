from oceantracker.particle_properties._base_particle_properties import ParticleProperty
import numpy as np
from oceantracker.util.parameter_checking import ParamValueChecker as PVC


class RadiusSpherical(ParticleProperty):
    """
    Assumes spherical particles with a homogenious density distribution.
    """
    
    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1e-6, float,doc_str='Particle property at the time of release')})
        self.class_doc(description='Particle fractal radius.')

    
    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        pass


class RadiusFractal(ParticleProperty):
    """
    Assumes spherical particles with a homogenious density distribution.
    """
    
    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1e-6, float,doc_str='Particle property at the time of release'),
                                  'alpha': PVC(22.25, float,doc_str='Fractal alpha'),
                                  'beta': PVC(0.429, float,doc_str='Fractal beta')})
        self.class_doc(description='Particle radius, manually updated.')

    

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        # calculate fractal properties based on Stemmann et al (2004)
        # presented in Adrian Burds (2013) paper

        # alpha_fractal = (4/3 * np.pi) ** (-1 / self.particle_fractal_dimension) * self.radius_unit_particle_m ** (1 - 3 / self.particle_fractal_dimension)
        # self.alpha_fractal = alpha_fractal * np.sqrt(0.6)
        # self.beta_fractal = 1. / self.particle_fractal_dimension

        alpha = self.params['alpha']
        beta = self.params['beta']

        volume_spherical = (4/3) * np.pi * part_prop['radius_spherical'].data[active]**3
        radius_fractal =  alpha * volume_spherical** beta

        self.set_values(radius_fractal, active)


class DensitySpherical(ParticleProperty):
    """
    Assumes spherical particles with a homogenious density distribution.
    """

    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1000., float,doc_str='Particle property at the time of release')})
        self.class_doc(description='Particle density, manually updated.')

    

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        pass


class DensityFractal(ParticleProperty):
    """
    Assumes fractal particles.
    """

    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1000., float,doc_str='Particle property at the time of release'),
                                  'alpha': PVC(22.25, float,doc_str='Fractal alpha'),
                                  'beta': PVC(0.429, float,doc_str='Fractal beta')})
        self.class_doc(description='Particle fractal density')

    

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        # calculate fractal properties based on Stemmann et al (2004)
        # presented in Adrian Burds (2013) paper

        # alpha_fractal = (4/3 * np.pi) ** (-1 / self.particle_fractal_dimension) * self.radius_unit_particle_m ** (1 - 3 / self.particle_fractal_dimension)
        # self.alpha_fractal = alpha_fractal * np.sqrt(0.6)
        # self.beta_fractal = 1. / self.particle_fractal_dimension

        alpha = self.params['alpha']
        beta = self.params['beta']

        volume_spherical = (4/3) * np.pi * part_prop['radius_spherical'].data[active]**3
        volume_fractal =  alpha * volume_spherical** beta

        density_fractal = part_prop['density_spherical'].data[active] * (volume_spherical / volume_fractal)

        self.set_values(density_fractal, active)


class StokesBasedBuoyancy(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({
            'initial_value': PVC(0, float,doc_str='Particle property at the time of release'),
            'mu': PVC(1e-3, float, doc_str='Dynamic viscosity of water'),
            'gravity': PVC(9.81, float, doc_str='Gravity acceleration'),            
        })
        self.class_doc(description='Particle buoyancy in m/s')

    def check_requirements(self):
        self.check_class_required_fields_prop_etc(required_props_list=['density_spherical', 'radius_spherical'])
  

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        si = self.shared_info
        part_prop = si.classes['particle_properties']
        
        # based on dynamic viscosty; assuming slow sinking
        # v = (2/9) * ((rho_p - rho_f) * g * radius**2) / mu
        # buoyancy = - (2/9) * ((part_prop['density_spherical'].data[active] - 1000) * self.params['gravity'] * part_prop['radius'].data[active]**2) / self.params['mu'] 
        
        # based on kinematic viscosty; assuming slow sinking
        # v = \frac{g}{\nu} \frac{(\rho_s - \rho_w)}{\rho_w} \frac{d^2}{18}
        # buoyancy = (1/18) * self.params['gravity'] * ((part_prop['density_spherical'].data[active] - 1000) / 1000) * (2*part_prop['radius'].data[active])**2
        # buoyancy = - (2/9) * self.params['gravity']/self.params['mu'] * ((part_prop['density_spherical'].data[active] - 1000) / 1000) * (part_prop['radius'].data[active])**2

        # alpha * (1/18) * (9.81 * (radius*2)**2 * (density - 1000)) / (kinematic_viscosity * density)
        buoyancy = - (0.25) * (9.81 * (part_prop['radius_spherical'].data[active]*2)**2 * (part_prop['density_spherical'].data[active] - 1000)) / (self.params['mu'] * part_prop['density_spherical'].data[active])

        self.set_values(buoyancy, active)


class PowerLawBasedBuoyancy(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({
            'radius': PVC('radius_spherical', str, doc_str='Name of the radius part prop to use'),
            'initial_value': PVC(0, float,doc_str='Particle property at the time of release'),
            'a': PVC(1, float, doc_str='scaling constant'),        
            'k': PVC(3, float, doc_str='exponential constant'),
        })
        self.class_doc(description='Particle buoyancy in m/s')

    # def check_requirements(self):
    #     self.check_class_required_fields_prop_etc(required_props_list=['density_spherical', 'radius'])
  
    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        radius = part_prop[self.params['radius']].data[active]
        buoyancy = - self.params['a'] * (radius**-self.params['k']) # m/d to m/s

        self.set_values(buoyancy, active)

class ParticleCollision(ParticleProperty):
    
    def __init__(self):
        super().__init__()
        self.add_default_params({ 'radius': PVC('radius_spherical', str, doc_str='Name of the radius part prop to use in collision'),
                                  'stickyness': PVC(0.1, float,doc_str='Chance of two colliding particles to stick together'),
                                  'spm_field': PVC('spm', str,doc_str='Name of the SPM field to use for collision detection'),
                                  'spm_radius': PVC(1e-6, float,doc_str='Radius of the SPM particles'),
                                  'spm_density': PVC(2650., float,doc_str='Density of the SPM particles'),
                                  'combine_density_method': PVC('fractal', str,possible_values=['spherical', 'fractal'],
                                                                doc_str='Method to use for combining densities. Options: spherical, fractal'),
                                  })
    
    # def check_requirements(self):redd_props_list=['density_spherical', 'radius', 'buoyancy'])

    def initial_setup(self):
        super().initial_setup()
        si = self.shared_info

        volume_of_ball = 4/3 * np.pi * (self.params['spm_radius'])**3
        mass_of_ball = volume_of_ball * self.params['spm_density']
        self.info['particles_per_kg'] = 1/mass_of_ball

        self._combine_density_method = self._combine_density_spherical


    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(0, new_part_IDs) # sets this properties values

        
    def update(self, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        local_spm_concentration = part_prop[self.params['spm_field']].data[active]
        # there is currently a "bug" in the input that makes the concentration of the SPM field
        # to be negative sometimes.
        # temporary fix: set negative values to 0
        local_spm_concentration[local_spm_concentration < 0] = 0

        particle_per_liter = local_spm_concentration * self.info['particles_per_kg']
        # make particles_per_liter
        particle_per_m3 = particle_per_liter * 1000
        
        # collision_kernel = self._collision_kernel_based_on_Delichatsios1975(self.params['spm_radius'], part_prop['radius'].data[active], self.info['average_relative_velocity'])
        # print(f'collision kernel by Delichatsios: {collision_kernel}')
        
        collision_radius = part_prop[self.params['radius']].data[active]
        collision_kernel = self._collision_kernel_based_on_burd2013(self.params['spm_radius'], collision_radius, 0.1)
        # print(f'collision kernel by Burd: {collision_kernel}')

        collision_frequency = (1/2) * particle_per_m3 * collision_kernel
        sticking_frequency = collision_frequency * self.params['stickyness']
        avg_coagulations = sticking_frequency*self.shared_info.settings['time_step']

        # we do not de-coagulate currently
        # to avoid large massive particles we stop particles above 1mm from coagulating
        avg_coagulations[part_prop['radius_spherical'].data[active] > 1e-3] = 0

        # roll for collision. 
        number_of_sticky_collisions = np.random.poisson(avg_coagulations)

        # add collision count to part_prop
        previous_collision_count = self.data[active]
        self.set_values(previous_collision_count + number_of_sticky_collisions, active)

        # adjust radius and density
        new_density, new_radius = self._combine_density_method(
            self.params['spm_radius'],
            part_prop['radius_spherical'].data[active],
            self.params['spm_density'],
            part_prop['density_spherical'].data[active],
            number_of_sticky_collisions)

        # update properties
        part_prop['radius_spherical'].set_values(new_radius, active)
        part_prop['density_spherical'].set_values(new_density, active)


    def _combine_density_spherical(self, radius_spm, radius_particle, density_spm, density_b, collisions_count):
        """
        Calculates the density of the combined sphere formed by merging
        two spheres with radii a and b and densities density_a and density_b.
        """
        # Volumes of the original spheres
        volume_spm = (4/3) * np.pi * radius_spm**3
        volume_particle = (4/3) * np.pi * radius_particle**3

        # Masses of the original spheres
        mass_spm = density_spm * volume_spm
        mass_particle = density_b * volume_particle

        # Combined volume
        combined_radius = (collisions_count*radius_spm**3 + radius_particle**3)**(1/3)
        combined_volume = (4/3) * np.pi * combined_radius**3

        # Combined density
        combined_density = (collisions_count*mass_spm + mass_particle) / combined_volume
        
        return combined_density, combined_radius


    # This is wrong/buggy as it updates the presumably spherical radius to a fractal radius
    # at each time step even tho is is a "fractal radius" after the first iteration.
    # def _combine_density_fractal(self, radius_spm, radius_particle, density_spm, density_b, collisions_count):
    #     """
    #     Calculates the density of the combined sphere formed by merging
    #     two spheres with radii a and b and densities density_a and density_b.
    #     """
    #     # Volumes of the original spheres
    #     volume_spm = (4/3) * np.pi * radius_spm**3
    #     volume_particle = (4/3) * np.pi * radius_particle**3

    #     # Masses of the original spheres
    #     mass_spm = density_spm * volume_spm
    #     mass_particle = density_b * volume_particle

    #     """
    #     Calculate the radius of a fractal particle.
    #     Partly based on Stemmann et al (2004).
    #     """
    #     alpha = self.alpha_fractal # 22.25
    #     beta = self.beta_fractal # 0.429

    #     combined_volume_spherical = (collisions_count*volume_spm + volume_particle)
    #     combined_volume_fractal =  alpha * combined_volume_spherical** beta

    #     combined_radius_fractal = (3/4 * combined_volume_fractal / np.pi)**(1/3)

    #     # Combined density
    #     combined_density = (collisions_count*mass_spm + mass_particle) / combined_volume_fractal
    #     return combined_density, combined_volume_fractal
    


    @staticmethod
    def _collision_kernel_based_on_burd2013(spm_radius, test_particle_radius, shear_gradient):
        """
        Calculate the collision kernel based on Adrians Burd model.

        Parameters:
        radius_i (float): The radius of particle i.
        radius_j (float): The radius of particle j.
        shear_gradient (float): The shear gradient of the fluid.

        Returns:
        float: The collision kernel.
        """

        return (4/3) * shear_gradient * (spm_radius + test_particle_radius)**3
        

    @staticmethod
    def _collision_kernel_based_on_Delichatsios1975(spm_radius, test_particle_radius, avg_velocity):
        
        cross_section =  np.pi * (spm_radius + test_particle_radius)**2
        kernel = cross_section * avg_velocity
        
        return kernel

