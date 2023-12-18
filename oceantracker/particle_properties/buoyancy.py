from typing import Any
from oceantracker.particle_properties._base_properties import ParticleProperty
import numpy as np
from oceantracker.util.parameter_checking import ParamValueChecker as PVC


class Density(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1000., float,doc_str='Particle property at the time of release')})
        self.class_doc(description='Particle density, manually updated.')

    

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        pass


class Radius(ParticleProperty):
    
    def __init__(self):
        super().__init__()
        self.add_default_params({ 'initial_value': PVC(1e-6, float,doc_str='Particle property at the time of release')})
        self.class_doc(description='Particle radius, manually updated.')

    

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        pass


class Buoyancy(ParticleProperty):

    def __init__(self):
        super().__init__()
        self.add_default_params({
            'initial_value': PVC(0, float,doc_str='Particle property at the time of release'),
            'mu': PVC(1e-3, float, doc_str='Dynamic viscosity of water'),
            'gravity': PVC(9.81, float, doc_str='Gravity acceleration'),            
        })
        self.class_doc(description='Particle buoyancy in m/s')

    def check_requirements(self):
        self.check_class_required_fields_prop_etc(required_props_list=['density', 'radius'])
  

    def initial_value_at_birth(self, new_part_IDs):
        self.set_values(self.params['initial_value'], new_part_IDs) # sets this properties values

    def update(self,active):
        # manually updated
        si = self.shared_info
        part_prop = si.classes['particle_properties']
        
        # v = (2/9) * ((rho_p - rho_f) * g * radius**2) / mu
        buoyancy = - (2/9) * ((part_prop['density'].data[active] - 1000) * self.params['gravity'] * part_prop['radius'].data[active]**2) / self.params['mu'] 

        self.set_values(buoyancy, active)


class ParticleCollision(ParticleProperty):
    
    def __init__(self):
        super().__init__()
        self.add_default_params({ 'stickyness': PVC(0.1, float,doc_str='Chance of two colliding particles to stick together'),
                                  'spm_field': PVC('spm', str,doc_str='Name of the SPM field to use for collision detection'),
                                  'spm_radius': PVC(1e-6, float,doc_str='Radius of the SPM particles'),
                                  'spm_density': PVC(2650., float,doc_str='Density of the SPM particles')})
    
    # def check_requirements(self):redd_props_list=['density', 'radius', 'buoyancy'])

    def initial_setup(self):
        super().initial_setup()
        si = self.shared_info
        particle= si.classes['particle_group_manager']

        volume_of_ball = 4/3 * np.pi * (self.params['spm_radius'])**3
        mass_of_ball = volume_of_ball * self.params['spm_density']
        self.info['particles_per_kg'] = 1/mass_of_ball

        horizontal_dispersion = self.shared_info.classes['dispersion'].params['A_H']
        vertical_dispersion = 0

        self.info['average_relative_velocity'] = np.sqrt(horizontal_dispersion**2 + vertical_dispersion**2)*0.798

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
        cross_section =  np.pi * (self.params['spm_radius'] + part_prop['radius'].data[active])**2
        avg_velocity = self.info['average_relative_velocity']

        collision_frequency = particle_per_m3 * cross_section * avg_velocity
        sicking_frequency = collision_frequency * self.params['stickyness']
        avg_coagulations = collision_frequency*self.params['stickyness']*self.shared_info.solver_info['model_time_step']

        # we do not de-coagulate currently
        # to avoid large massive particles we stop particles above 1mm from coagulating
        avg_coagulations[part_prop['radius'].data[active] > 1e-3] = 0

        # roll for collision. 
        number_of_sticky_collisions = np.random.poisson(avg_coagulations)

        # any collision?
        #any_collisions = (number_of_sticky_collisions > 0).any()

        # add collision count to part_prop
        previous_collision_count = part_prop['collision_very_fine_silt'].data[active]
        self.set_values(previous_collision_count + number_of_sticky_collisions, active)

        # adjust radius and density
        new_density, new_radius = self._combine_density(
            self.params['spm_radius'],
            part_prop['radius'].data[active],
            self.params['spm_density'],
            part_prop['density'].data[active],
            number_of_sticky_collisions)

        # update properties
        part_prop['radius'].set_values(new_radius, active)
        part_prop['density'].set_values(new_density, active)


    @staticmethod
    def _combine_density(radius_spm, radius_particle, density_spm, density_b, collisions_count):
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

