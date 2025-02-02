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

        buoyancy = settling_velocity_kriest_dPAM(part_prop[self.params['radius']].data[active],self.params['a'],self.params['k'])

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
                                  'coagulation_kernel': PVC('curviliniar_shear', str,possible_values=['rectilinear_shear', 'curviliniar_shear', 'curviliniar_shear & curvilinear_diff_settling'],
                                                             doc_str='Method to use for calculating the coagulation kernel aka as the probability of two particles colliding and sticking together. Options: rectilinear_shear, curviliniar_shear'),
                                  'initial_radius': PVC(1e-5, float,doc_str='Organic aggregate size at the time of release'),})
    
    # def check_requirements(self):redd_props_list=['density_spherical', 'radius', 'buoyancy'])

    def initial_setup(self):
        super().initial_setup()
        si = self.shared_info

        volume_of_ball = 4/3 * np.pi * (self.params['spm_radius'])**3
        mass_of_ball = volume_of_ball * self.params['spm_density']
        self.info['particles_per_kg'] = 1/mass_of_ball

        self._combine_density_method = self._combine_density_spherical

        # set coagulation kernel based on the self.params['coagulation_kernel']
        if self.params['coagulation_kernel'] == 'rectilinear_shear':
            self.coagulation_kernel = self._rectilinear_shear
        elif self.params['coagulation_kernel'] == 'curviliniar_shear':
            self.coagulation_kernel = self._curviliniar_shear
        elif self.params['coagulation_kernel'] == 'curviliniar_shear & curvilinear_diff_settling':
            self.coagulation_kernel = self._curviliniar_shear_and_curvilinear_diff_settling
        else:
            # raise error
            raise ValueError(f"Invalid coagulation kernel method: {self.params['coagulation_kernel']}")





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
        
        collision_kernel = self.coagulation_kernel(self.params['spm_radius'],active, initial_particle_size=self.params['initial_radius'])

        collision_frequency = particle_per_m3 * collision_kernel
        sticking_frequency = collision_frequency * self.params['stickyness']
        avg_coagulations = sticking_frequency*self.shared_info.settings['time_step']

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
        combined_volume = (collisions_count*volume_spm + volume_particle)
        combined_radius = (3/4 * combined_volume / np.pi)**(1/3)

        # Combined density
        combined_density = (collisions_count*mass_spm + mass_particle) / combined_volume
        
        return combined_density, combined_radius
   

    def _rectilinear_shear(self, spm_radius, test_particle_radius, active, shear_gradient=0.1):
        """
        Calculate the collision kernel based on Adrians Burd model.

        Parameters:
        radius_i (float): The radius of particle i.
        radius_j (float): The radius of particle j.
        shear_gradient (float): The shear gradient of the fluid.

        Returns:
        float: The collision kernel.
        """
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        test_particle_radius = part_prop[self.params['radius']].data[active]

        beta = (4/3) * shear_gradient * (spm_radius + test_particle_radius)**3

        return beta


    def _curviliniar_shear(self, spm_radius, active, epsilon = 0.0026964394, nu = 1e-6, initial_particle_size = 1e-5):
        """
        Calculate coagulation kernel with particles "avoiding" other particles
        due to streamline curving around the particles.

        Parameters:
        - radius_i: Radius of particle i
        - radius_j: Radius of particle j
        - epsilon: turbulence kinetic energy dissipation rate
        - nu: kinematic viscosity
        
        Returns:
        - beta: coagulation kernel value
        """

        si = self.shared_info
        part_prop = si.classes['particle_properties']

        test_particle_radius = part_prop[self.params['radius']].data[active]
        tke_diss = part_prop['TKE_dissipation_rate'].data[active]
        # set all elements in tke_diss to zero if negative (hotfix)
        tke_diss[tke_diss < 0] = 0
                

        ratio_organic_inorganic = initial_particle_size**3 / test_particle_radius**3 # volume ratio
        radius_gyration = (spm_radius + test_particle_radius) * 1.3 #self.radius_of_sphere_to_radius_of_gyration
        particle_ratio = np.minimum(spm_radius,test_particle_radius) / np.maximum(spm_radius,test_particle_radius)
        coag_efficiency = 1 - (1 + 5*particle_ratio + 2.5*particle_ratio**2) / (1 + particle_ratio)**5
        beta = np.sqrt(8*np.pi*tke_diss/15/nu) * coag_efficiency * ratio_organic_inorganic * radius_gyration**3

        return beta
        

    def _curvilinear_diff_settling(self, spm_radius, active, initial_particle_size = 1e-5):
        """
            \beta^{C}_{ds} &= \frac{1}{2}\pi r^{2}_{i} | v_{i} - v_{j} |
        """

        si = self.shared_info
        part_prop = si.classes['particle_properties']

        test_particle_radius = part_prop[self.params['radius']].data[active]

        settling_velocity_spm = settling_velocity_kriest_dPAM(spm_radius)
        settling_velocity_test_particle = settling_velocity_kriest_dPAM(test_particle_radius)

        ratio_organic_inorganic = initial_particle_size**3 / test_particle_radius**3 # volume ratio
        radius_gyration = (spm_radius + test_particle_radius) * 1.3 #self.radius_of_sphere_to_radius_of_gyration

        beta = ratio_organic_inorganic * np.pi * radius_gyration**2 * np.abs(settling_velocity_spm - settling_velocity_test_particle)

        return beta

    def _curviliniar_shear_and_curvilinear_diff_settling(self, spm_radius, active, epsilon = 0.0026964394, nu = 1e-6, initial_particle_size = 5e-5):
        return self._curviliniar_shear(spm_radius, active, epsilon, nu, initial_particle_size) \
               + self._curvilinear_diff_settling(spm_radius, active, initial_particle_size)


def settling_velocity_kriest_dPAM(radius, a = 942, k = 1.17):

        # transform it to be consistent with kriest 2002 (dense particle)
        radius = radius * 1e2 # [cm]
        
        buoyancy = - a * (radius**k) # [m/d]
        buoyancy = buoyancy / 86400 # [m/s]

        return buoyancy


def settling_velocity_sedimorph(radius):

    # | Very fine silt    | > 8       | 8 – 4        | 0.31·10⁻³       |
    if radius == 6e-6/2:
        return 3.1e-6
    # | Fine silt         | > 7       | 16 – 8       | 0.11·10⁻²        |
    elif radius == 12e-6/2:
        return 1.1e-5
    # | Medium silt       | > 6       | 31 – 16      | 0.51·10⁻²        |
    elif radius == 24e-6/2:
        return 5.1e-5
    # | Coarse silt       | > 5       | 62 – 31      | 0.19·10⁻¹        |
    elif radius == 47e-6/2:
        return 1.9e-4
    # | Very fine sand    | > 4       | 125 – 62     | 0.07              |
    elif radius == 94e-6/2:
        return 7e-4
    else:
        raise ValueError(f"Invalid radius: {radius}")        

