import numpy as np
import oceantracker.util.particle_operations as part_op
from oceantracker.user_trajectory_modifiers._trajectory_modifers_base import \
    _BaseTrajectoryModifier
from oceantracker.util.polygonUtil import InsidePolygon

from numba import njit


class gathering(_BaseTrajectoryModifier):
    # fallows particles to freeze if inside a polygon
    def __init__(self, param_dict={}):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'name': self.__class__.__name__,
                                 'requires_3D': False,
                                 'particle mass': 1
                                 })
        self.update_params(param_dict)  # this is required in children
        self.requiements_txt = 'requirements not given'

    def initialize(self, **kwargs):
        super().initialize(**kwargs)

        particle = self.shared_info.pointers['particle']
        interp = self.shared_info.pointers['interp']
        
        # add arbitrary (mass) to particles (for dev purposes)
        try:
            particle.part_prop['particle mass']
        except KeyError:
            particle.add_time_varying_prop_manualUpdate(
                'particle mass', {'data_type': np.int8})
        
        # calculates the list of surrounding triangles for each triangles
        try:
            interp.surrounding
        except AttributeError:
            self.code_timer.start('build_surrounding_cells')
            interp.surrounding = build_surrounding_cells(interp.grid_data['triangles'])
            self.code_timer.stop('build_surrounding_cells')

    def can_be_added(self, hindcast_is3D):
        # required method to decide if the operation can be add to simulation
        if not hindcast_is3D:
            return True

    def update(self, buffer_index, time, active):
        solver = self.shared_info.pointers['solver']
        x2 = solver.x_new

        # calculates velocity change due to proximity to other particles
        v_aggregation = self.aggregation_trajectory_modifier(active)

        # and adds it to net velocity by a simple euler forward
        dt = self.shared_info.model_time_step*solver.info['n_sub_steps']
        part_op.add(x2, x2, v_aggregation, dt, active)

    def find_surrounding_particles(self, active):
        """
        Looks thru all 'surrounding' triangles - aka touching triangles - to 
        find and return all particles within the surrounding of each particle
        """
        
        self.code_timer.start('find_surrounding_particles')

        # index of cell containing each particles
        particle = self.shared_info.pointers['particle']
        interp = self.shared_info.pointers['interp']

        # get the index of cell in which a particle currently is
        n_cell = particle.get_partProp_dataPtr('n_cell')[active]
        # surrounding cells of each cell
        surrounding_cells = interp.surrounding

        surrounding_particles = []
        # Carefull! I am not sure that it is guaranteed that the cells in
        # surrounding_cells have the the same ordering as the one in "grid"
        for ii, each_particle_cell in enumerate(n_cell):
            # maybe i should index the particles by cell to save time?
            cells_to_check = surrounding_cells[each_particle_cell]
            surrounding_current_particle = np.array([],dtype=np.uint8)
            for each_surrounding_cell in cells_to_check:
                    particles_in_current_cell = active[n_cell == each_surrounding_cell]
                    surrounding_current_particle = np.append(surrounding_current_particle,particles_in_current_cell)
            # removing the current particle from the surrounding particles
            # (technically it is within the distance limit of itself)
            surrounding_current_particle = surrounding_current_particle[surrounding_current_particle != ii]
            surrounding_particles.append(surrounding_current_particle)

        self.code_timer.stop('find_surrounding_particles')
        
        return surrounding_particles

    def aggregation_trajectory_modifier(self, active):
        """
        calculates the velocity change of particles induces by the proximity to 
        other particles
        """

        self.code_timer.start('aggregation_trajectory_modifier')

        self.surrounding_particles = self.find_surrounding_particles(active)
        x_old = self.shared_info.pointers['solver'].x_old

        v_aggregation = np.zeros(x_old.shape)

        for ii, current_particle in enumerate(x_old[active]):
            # calculate distance between particles
            # current surrounding particle
            cur_sur_part = x_old[self.surrounding_particles[ii]]
            if len(cur_sur_part) == 0:
                v_aggregation[active[ii]] = 0
            else:
                particle_dist_vec = cur_sur_part - current_particle
                # currently we keep all particles

                # calculate force based on distance and some attraction parameter
                forces = lennard_jones_potential(particle_dist_vec, 1e0, 1e15)

                # sum forces
                force = np.sum(forces, axis=0)

                # translate forces into a net velocity change
                #   this one depends on the concept of mass
                #   currently we assume that the mass is one, making each particle
                #   equally heavy.
                v_aggregation[active][ii] = force
        
        self.code_timer.stop('aggregation_trajectory_modifier')

        return v_aggregation


class merging(_BaseTrajectoryModifier):
    # fallows particles to freeze if inside a polygon
    def __init__(self, param_dict={}):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'name': self.__class__.__name__,
                                 'requires_3D': False,
                                 'merging_distance': 1e3,
                                 'merging_probability': 1
                                 })
        self.update_params(param_dict)  # this is required in children
        self.requiements_txt = 'requirements not given'

    def initialize(self, **kwargs):
        super().initialize(**kwargs)

        particle = self.shared_info.pointers['particle']
        interp = self.shared_info.pointers['interp']
        
        # add arbitrary (mass) to particles (for dev purposes)
        try:
            particle.part_prop['particle mass']
        except KeyError:
            particle.add_time_varying_prop_manualUpdate(
                'particle mass', {'data_type': np.int8})

        # calculates the list of surrounding triangles for each triangles
        try:
            interp.surrounding
        except AttributeError:
            self.code_timer.start('build_surrounding_cells')
            interp.surrounding = build_surrounding_cells(interp.grid_data['triangles'])
            self.code_timer.stop('build_surrounding_cells')

    def can_be_added(self, hindcast_is3D):
        # required method to decide if the operation can be add to simulation
        if not hindcast_is3D:
            return True

    def update(self, buffer_index, time, active):
        # joins close particles together
        self.merge_particles(active)

    def merge_particles(self,active):


        self.code_timer.start('merge_particles')

        particle    = self.shared_info.pointers['particle']
        interp      = self.shared_info.pointers['interp']
        status      = particle.part_prop['status'].data
        mass        = particle.part_prop['particle mass'].data
        x = self.shared_info.pointers['solver'].x_old
        
        
        #for ii,current_active in enumerate(active):
        ii = 0
        while ii < len(active):
            current_active = active[ii]
            
            self.surrounding_particles = gathering().find_surrounding_particles(active)
            # index in relation to particle list

            ## ii or current_active?
            considered_mergers = self.surrounding_particles[ii]
            #considered_mergers = active[index_potential_mergers]
            if len(considered_mergers) != 0:
                print(f'initial considered: {considered_mergers}') 
                print(f'current active: {current_active}')
                if True: # selection of mergers based on range and chance
                # distance between considered particles
                    particle_dist_vec = x[considered_mergers] - x[current_active]
                    # considered particles inside range
                    considered_mergers = \
                        considered_mergers[np.linalg.norm(particle_dist_vec,axis=1) 
                                            < self.params['merging_distance']]

                    # merge with a certain probability
                    considered_mergers = \
                        considered_mergers[np.random.rand(len(considered_mergers))
                                            < self.params['merging_probability']]

                    selected_mergers = considered_mergers
                    print(f'new: {selected_mergers}')

                if False: # pure random selection
                    considered_mergers = self.surrounding_particles[ii]
                    # select a random amount of particles to merge
                    n_to_merge = np.random.randint(len(considered_mergers))
                    # indices in relation to considered mergers
                    selected_indices = np.random.choice(np.arange(len(considered_mergers)),
                                                        n_to_merge,replace=False)
                    # indices in relation to particle list
                    selected_mergers = considered_mergers[selected_indices]
                    print(f'old: {selected_mergers}')

                #cur_sur_part = self.particle_location[self.surrounding_particles[ii]]
                #particle_dist_vec = cur_sur_part - current_particle

                merger_particle_masses = \
                    particle.part_prop['particle mass'].data[selected_mergers]

                joined_mass = np.sum(merger_particle_masses)
                part_op.add_scaler(mass, joined_mass, np.array([current_active]))
                part_op.scalar_set(mass, 0 , selected_mergers)
                part_op.scalar_set(status,particle.part_status['dead'],
                                selected_mergers)
                
                # i don't know if this works like this
                #for jj in range(len(self.surrounding_particles)):
                #    # drop those indices (rel to p) which have been merged
                #    for kk in selected_mergers:  
                #        self.surrounding_particles[jj] = \
                #            self.surrounding_particles[jj][self.surrounding_particles[jj] != kk]
                for jj in selected_mergers:
                    active = active[active != jj]

            ii += 1
        
        self.code_timer.stop('merge_particles')


@njit
def build_surrounding_cells(triangles):
    """ 
    Takes a list of triangles which contain the nodes of each triangle
    for which if finds all triangles that "touch" it i.e. have a common
    node

    Returns:
    --------
        surrounding : list
            list of triangles which lists all touching triangles
    """

    # take a look at the draft from ross
    surrounding = []
    for t1 in triangles:
        # tri has the nodes of the triangles
        surrounding_t1 = []

        for ii in np.arange(len(triangles)):
            t2 = triangles[ii]
            # if they have a common node - they touch
            touching_idx = _are_tri_touching(ii,t1,t2)
            if touching_idx != -1: surrounding_t1.append(ii)
        surrounding.append(surrounding_t1)

    return surrounding

@njit
def _are_tri_touching(ii,t1,t2):
    for t2_node in t2:
        if t2_node in t1:
            return ii
    return -1


def lennard_jones_potential(vector, sigma, eta):
    # catch warning if an item in r has 0 distance
    r = np.linalg.norm(vector, axis=1)
    potential = 4*eta*(sigma/r)**12 - (sigma/r)**6
    force = potential.reshape((len(r), 1)) * vector/r.reshape((len(r), 1))
    return force
