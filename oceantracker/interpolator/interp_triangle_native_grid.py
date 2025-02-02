# linear interploation for triangles in both space and time
#todo  are BC cords as np.float32, faster as lower memory transfer demand and good enough?
import numpy as np
from scipy.spatial import cKDTree

from oceantracker.interpolator._base_interp import _BaseInterp
from oceantracker.util import basic_util
from oceantracker.util.profiling_util import function_profiler
from time import perf_counter
from oceantracker.util import numpy_util
from oceantracker.common_info_default_param_dict_templates import  cell_search_status_flags, particle_info
from oceantracker.interpolator.util import triangle_interpolator_util as tri_interp_util ,  triangle_eval_interp
from oceantracker.particle_properties.util import  particle_operations_util

from oceantracker.util.parameter_checking import  ParamValueChecker as PVC


class  InterpTriangularNativeGrid_Slayer_and_LSCgrid(_BaseInterp):

    # uses tweaked sci py which allows using start triangle location

    def __init__(self):
        # set up info/attributes
        super().__init__()
        self.add_default_params({'bc_walk_tol': PVC(1.0e-4, float,min = 0.),
                                 'max_search_steps': PVC(500,int, min =1)})
        self.info['current_buffer_index'] = np.zeros((2,), dtype=np.int32)

    #@function_profiler(__name__)
    def initial_setup(self, grid):
        super().initial_setup()  # children must call this parent class to default shared_params etc
        params = self.params
        si= self.shared_info


        # make barcentric transform matrix for the grid,  typee to match numba signature
        t0 = perf_counter()
        grid['bc_transform'] = tri_interp_util.get_BC_transform_matrix(grid['x'].data, grid['triangles'].data).astype(np.float64)
        si.msg_logger.progress_marker('built barycentric-transform matrix', start_time=t0)
        # build kd tree for initial triangle find node nearest node locations
        #todo delete xy_centriod = np.mean(grid['x'][grid['triangles']], axis=1)
        #todo deleteself.KDtree = cKDTree(xy_centriod)
        self.KDtree = cKDTree(grid['x'])

    def add_part_prop_for_book_keeping(self):
        si = self.shared_info
        pgm = si.classes['particle_group_manager']

        pgm.add_particle_property('n_cell', 'manual_update', dict(write=False, dtype=np.int32, initial_value=0))  # start with cell number guess of zero
        pgm.add_particle_property('n_cell_last_good', 'manual_update', dict(write=False, dtype=np.int32, initial_value=0))  # start with cell number guess of zero
        pgm.add_particle_property('cell_search_status', 'manual_update', dict(write=False, initial_value=cell_search_status_flags['ok'], dtype=np.int8))

        pgm.add_particle_property('bc_cords', 'manual_update', dict(write=False, initial_value=0., vector_dim=3, dtype=np.float32))

        # BC walk info
        if si.is3D_run:
            # space to record vertical cell for each particles' triangle at two timer steps  for each node in cell containing particle
            # used to do 3D time dependent interpolation
            pgm.add_particle_property('nz_cell', 'manual_update', dict(write=False, dtype=np.int32, initial_value=0))  # todo  create  initial serach for vertical cell
            pgm.add_particle_property('z_fraction', 'manual_update', dict(write=False, dtype=np.float32, initial_value=0.))
            pgm.add_particle_property('z_fraction_water_velocity', 'manual_update', dict(write=False, dtype=np.float32, initial_value=0., description=' thickness of bottom layer in metres, used for log layer velocity interp in bottom layer'))

    def final_setup(self, grid):

        # set up a grid class,part_prop and vertical cell find functions to minimise numba function arguments
        si = self.shared_info
        info = self.info

        # create particle properties to  store history of current triangle  for reuse


        # set up place for walk info failures
        info['tri_walk_full_failures'] = []

        # set up walk count vector and map to info
        self.walk_counts=np.zeros((8,),dtype=np.int64)
        wc = self.walk_counts
        info.update(dict(   particles_located_by_walking = wc[0:1],
                            number_of_triangles_walked = wc[1:2],
                            longest_triangle_walk = wc[2:3],
                            nans_encountered_triangle_walk = wc[3:4],
                            triangle_walks_retried = wc[4:5],
                            particles_killed_after_triangle_walk_retry_failed = wc[5:6],
                            total_vertical_steps_walked = wc[6:7],
                            longest_vertical_walk = wc[7:8])
                            )

        # location of each vertex
        grid['x_vertex'] = np.stack((grid['x'][grid['triangles'], 0], grid['x'][grid['triangles'], 1]), axis=2)

        # build triangle walk array of structures
        grid['tri_walk_AOS'] =numpy_util.numpy_array_of_structures_from_dict(
                                     dict(bc_transform= grid['bc_transform'],
                                        adjacency=  grid['adjacency']),
                                        )



    def fix_bad_cell_search(self, current_buffer_steps,fractional_time_steps, open_boundary_type,  active):
        si = self.shared_info
        part_prop = si.classes['particle_properties']

        # deal with open boundary
        sel = part_prop['cell_search_status'].find_subset_where(active, 'eq', cell_search_status_flags['outside_domain'], out=self.get_partID_subset_buffer('cell_status'))
        if sel.size > 0:
            if open_boundary_type > 0:
                part_prop['status'].set_values(particle_info['status_flags']['outside_open_boundary'], sel)
                part_prop['n_cell'].copy('n_cell_last_good', sel)  # move back the cell, but not the location
            else:
                # outside so move back
                part_prop['x'].copy('x_last_good', sel)  # move back location
                part_prop['n_cell'].copy('n_cell_last_good', sel)  # move back the cell

        # do move backs for blocked and bad
        sel = part_prop['cell_search_status'].find_subset_where(active, 'lt', cell_search_status_flags['ok'], out=self.get_partID_subset_buffer('cell_status'))
        if sel.size > 0:
            part_prop['x'].copy('x_last_good', sel)  # move back location
            part_prop['n_cell'].copy('n_cell_last_good', sel)  # move back the cell


            # debug_util.plot_walk_step(xq, si.classes['reader'].grid, part_prop)



    #@function_profiler(__name__)
    def _interp_field2D(self,field_name, field_instance,grid, current_buffer_steps, fractional_time_steps,
                                n_cell, bc_cords,  output, active):
        # interp reader field_name inplace to particle locations to same time and memory
        # output can optionally be redirected to another particle property name different from  reader's field_name
        # particle_prop_name
        # in place evaluation of field interpolation
        si = self.shared_info

        triangles = grid['triangles']


        if field_instance.is_time_varying():
            self._time_dependent_2Dfield(current_buffer_steps, fractional_time_steps,output,
                                                   field_instance, triangles,
                                                   n_cell, bc_cords, active)
        else:
            self._time_independent_2Dfield(output, field_instance, triangles,
                                                 n_cell, bc_cords,
                                                 active)

    def _time_independent_2Dfield(self,output, field_instance, triangles, n_cell, bc_cords, active):
        data =  field_instance.data
        if data.shape[3] ==1:
            #scalar, eg water depth
            triangle_eval_interp.time_independent_2Dfield_scalar(output,data, triangles, n_cell, bc_cords, active)
        elif data.shape[3] == 2:
                # vector, is this used?
                triangle_eval_interp.time_independent_2Dfield_vector(output, data, triangles, n_cell, bc_cords, active)

    def _time_dependent_2Dfield(self,current_buffer_steps, fractional_time_steps,output,
                                                   field_instance, triangles,
                                                   n_cell, bc_cords, active):
        si = self.shared_info
        data =  field_instance.data

        if data.shape[3] ==1:
            #scalar, eg water depth
            triangle_eval_interp.time_dependent_2Dfield_scalar(current_buffer_steps, fractional_time_steps,output,
                                                   data, triangles,
                                                   n_cell, bc_cords, active)
        elif data.shape[3] == 2:
                # vector, eg 2D water velocity
                triangle_eval_interp.time_dependent_2Dfield_vector(current_buffer_steps, fractional_time_steps,output,
                                                   data, triangles,
                                                   n_cell, bc_cords, active)
        else:
            si.msg_logger.msg(f'2D interpolation currently only works for scalar and vector fields "{field_instance.info["name"]}" has {data.shape[3]} components',
                              fatal_error=True, exit_now=True)

    #@function_profiler(__name__)
    def _interp_field3D(self, field_name, field_instance,grid, current_buffer_steps,fractional_time_steps,
                         n_cell,nz_cell,z_fraction, bc_cords,output,active):
        # interp reader field_name inplace to particle locations to same time and memory
        # output can optionally be redirected to another particle property name different from  reader's field_name
        # particle_prop_name
        # in place evaluation of field interpolation
        si = self.shared_info

        triangles = grid['triangles']

        part_prop = si.classes['particle_properties']
        n_cell = part_prop['n_cell'].data
        bc_cords = part_prop['bc_cords'].data
        nz_cell = part_prop['nz_cell'].data

        if field_name == 'water_velocity':
            # fractions for water vel. are log layer in bottom cell
            z_fraction = part_prop['z_fraction_water_velocity'].data
        else:
            z_fraction = part_prop['z_fraction'].data


        if field_instance.is_time_varying():
            if 'sigma' in grid:
                self._time_dependent_3Dfield_sigma_grid(current_buffer_steps, fractional_time_steps,
                            field_instance, triangles,
                           n_cell, bc_cords, nz_cell, z_fraction,
                          output, active)
            else:
                triangle_eval_interp.time_dependent_3Dfield_LSC_grid(current_buffer_steps, fractional_time_steps,
                                                                     field_instance.data,
                                                                     triangles,  grid['bottom_cell_index'],
                                                                     n_cell, bc_cords, nz_cell, z_fraction,
                                                                     basic_util.atLeast_Nby1(output), active)
        else:
            # 3D non-time varying
            # todo eval_interp3D_timeIndependent not implented for 3D non-time varying fields
            raise Exception('eval_field_interpolation_at_particle_locations : spatial interp using eval_interp3D_timeIndependent not implemented yet ')

    def _time_dependent_3Dfield_sigma_grid(self,current_buffer_steps, fractional_time_steps,
                            field_instance,       triangles,
                           n_cell, bc_cords, nz_cell, z_fraction,
                           output, active):
        si = self.shared_info

        data = field_instance.data
        if data.shape[3] == 1:            # scalar
            triangle_eval_interp.time_dependent_3Dfield_scalar_sigma_grid(current_buffer_steps, fractional_time_steps,
                                                           data, triangles,
                                                           n_cell, bc_cords, nz_cell, z_fraction,
                                                           output, active)
        elif data.shape[3] == 3:
            # scalar
            triangle_eval_interp.time_dependent_3Dfield_vector_sigma_grid(current_buffer_steps, fractional_time_steps,
                                                                          data, triangles,
                                                                          n_cell, bc_cords, nz_cell, z_fraction,
                                                                          output, active)
        else:
            si.msg_logger.msg(f'3D interpolation currently only works for scalar and vector fields "{field_instance.info["name"]}" has {data.shape[3]} components',
                              fatal_error=True, exit_now=True)

    #@function_profiler(__name__)
    def eval_field_interpolation_at_given_locations(self,field_name, field_instance,grid, reader, x,
                                                    current_buffer_steps, fractional_time_steps,time=None, output=None, n_cell=None):
        # in  evaluation of field interpolation at specific locations, ie not particle locations
        # todo only time_dependent_2Dfield  working - eval_field_interpolation_at_given_locations
        # todo add time dependence/ time fractions
        # does this over write paricle props??
        si = self.shared_info
        info = self.info

        part_prop = si.classes['particle_properties']

        # is no output name given particle property for output is same as hindcast field_name
        if output is None:
            output = np.full((x.shape[0], field_instance.data.shape[3]), np.nan)  if field_instance.data.shape[3] > 1else np.full((x.shape[0],), np.nan)

        if n_cell is None:
            n_cell = self.initial_horizontal_cell(grid, x)

        # get bc cords for the cells
        bc_cords = self.get_bc_cords(grid, x, n_cell)
        active = np.arange(x.shape[0])
        if time is  None:
            # not time varying
            nb = None
            nt = None
        else:
            nt = reader.time_to_hydro_model_index(time)
            nb = reader.buffer_index_to_buffer_offset(nt)

        if field_instance.is3D():
            nz_cell = part_prop['nz_cell'].data
            z_fraction = part_prop['z_fraction'].data
            bottom_cell_index = grid['bottom_cell_index']
            #todo not working in 3D or equal sigma layers
            self._interp_field3D(field_instance, n_cell, bc_cords, nz_cell,
                                 z_fraction,bottom_cell_index, active, output)
        else:

            if field_instance.is_time_varying():
                triangle_eval_interp.time_dependent_2Dfield(nb,  fractional_time_steps, basic_util.atLeast_Nby1(output),
                                                            field_instance.data, grid['triangles'],
                                                            n_cell, bc_cords,
                                                            active)
            else:
                self._time_independent_2Dfield(output,field_instance, grid['triangles'],
                                                    n_cell, bc_cords, active)
        return output

    
    def find_vertical_cell(self, grid, fields, xq, current_buffer_steps,fractional_time_steps, active):
        # locate vertical cell in place

        si= self.shared_info

        part_prop = si.classes['particle_properties']
        n_cell= part_prop['n_cell_last_good'].data
        status = part_prop['status'].data
        bc_cords = part_prop['bc_cords'].data


        t0 = perf_counter()
        nz_cell = part_prop['nz_cell'].data
        z_fraction = part_prop['z_fraction'].data
        z_fraction_water_velocity = part_prop['z_fraction_water_velocity'].data

        if 'sigma' in grid:
                tri_interp_util.get_depth_cell_sigma_layers(xq,
                                        grid['triangles'],
                                        #grid['water_depth_triangles'],
                                        fields['water_depth'].data.ravel(),
                                        fields['tide'].data,
                                        si.minimum_total_water_depth,
                                        grid['sigma'], grid['sigma_map_nz_interval_with_sigma'],
                                        n_cell, status, bc_cords, nz_cell, z_fraction, z_fraction_water_velocity,
                                        current_buffer_steps, fractional_time_steps,
                                        active, si.z0)
        else:
            # natve slayer option
            tri_interp_util.get_depth_cell_time_varying_Slayer_or_LSCgrid(xq,
                                        grid['triangles'],grid['zlevel'],grid['bottom_cell_index'],
                                        #grid['triangles'], grid['zlevel_vertex'], grid['bottom_cell_index'],
                                        n_cell, status, bc_cords,nz_cell,z_fraction,z_fraction_water_velocity,
                                        current_buffer_steps,fractional_time_steps,
                                        self.walk_counts,
                                        active,  si.z0)
        si.block_timer('Find cell, vertical walk', t0)

    def find_hori_cell(self,grid, fields,  xq, current_buffer_steps,fractional_time_steps,open_boundary_type, active):
        # do cell walk
        si= self.shared_info
        t0= perf_counter()
        info = self.info
        part_prop = si.classes['particle_properties']

        n_cell = part_prop['n_cell'].data
        bc_cords = part_prop['bc_cords'].data
        cell_search_status = part_prop['cell_search_status'].data
        params = self.params

        # used 2D or 3D walk chosen above
        if False:
            x_last_good = part_prop['x_last_good'].data
            n_cell_last_good = part_prop['n_cell_last_good'].data
            status = part_prop['status'].data
            tri_interp_util.BCwalk_with_move_backs(
                           xq,
                           grid['tri_walk_AOS'], grid['dry_cell_index'],
                           x_last_good, n_cell, n_cell_last_good, status, bc_cords,
                           self.walk_counts,
                           params['max_search_steps'], params['bc_walk_tol'], open_boundary_type,si.settings['block_dry_cells'],
                           active)
            sel = part_prop['status'].find_subset_where(active, 'eq', si.particle_status_flags['cell_search_failed'], out=self.get_partID_subset_buffer('B1'))

        else:
            tri_interp_util.BCwalk(
                xq,
                grid['tri_walk_AOS'], grid['dry_cell_index'],
                n_cell, cell_search_status,bc_cords,
                self.walk_counts,
                params['max_search_steps'], params['bc_walk_tol'], open_boundary_type, si.settings['block_dry_cells'], active)

            sel = part_prop['cell_search_status'].find_subset_where(active, 'eq', cell_search_status_flags['failed'], out=self.get_partID_subset_buffer('B1'))

        if sel.size > 0:
            # si.msg_logger.msg(f'Search retried for {sel.size} cells')
            info['triangle_walks_retried'] += sel.size
            new_cell = self.initial_horizontal_cell(grid, xq[sel, :])
            part_prop['n_cell'].set_values(new_cell, sel)
            tri_interp_util.BCwalk(
                xq,
                grid['tri_walk_AOS'], grid['dry_cell_index'],
                n_cell, cell_search_status, bc_cords,
                self.walk_counts,
                params['max_search_steps'], params['bc_walk_tol'], open_boundary_type, si.settings['block_dry_cells'], sel)


            # recheck for additional failures
            sel = part_prop['cell_search_status'].find_subset_where(active, 'eq', cell_search_status_flags['failed'], out=self.get_partID_subset_buffer('B1'))

            if sel.size > 0:
                wf = {'x0': part_prop['x_last_good'].get_values(sel),
                      'xq': part_prop['x'].get_values(sel)}

                info['tri_walk_full_failures'].append(wf)
                info['particles_killed_after_triangle_walk_retry_failed'] += sel.size  # total failed walks
                si.msg_logger.msg('walks too long after kd retry- killed ' + str(sel.shape[0]) + ' particles', warning=True, tabs=0,
                                  hint='Try decreasing time step or increasing interpolator parameter "max_search_steps", current value =' + str(self.params['max_search_steps']))
                # make notes for log file enabling follow up
                si.msg_logger.msg('particle locations of failed walks, first 3 or less ', warning=True, tabs=2)
                si.msg_logger.msg(' location xq =' + str(xq[sel[:3], :].tolist()), warning=True, tabs=2)
                si.msg_logger.msg(' x_old =' + str(part_prop['x_last_good'].data[sel[:3], :].tolist()), warning=True, tabs=2)
                # kill particles
                part_prop['status'].set_values(si.particle_status_flags['dead'], sel)

        si.block_timer('Find cell, horizontal walk', t0)



    #@function_profiler(__name__)
    def initial_horizontal_cell(self, grid, xq):
        # find nearest cell
        si=self.shared_info
        t0 = perf_counter()

         # find nearest node
        dist, nodes = self.KDtree.query(xq[:, :2])
        nodes = nodes.astype(np.int32)  # KD tree gives int64,need for compatibility of types

        # look in triangles attached to each node for tri containing the point
        #t0= perf_counter()
        n_cell= tri_interp_util.check_if_point_inside_triangle_connected_to_node(xq[:, :2], nodes,
                                     grid['node_to_tri_map'], grid['tri_per_node'],  grid['bc_transform'], self.params['bc_walk_tol'])
        #print('node find ', perf_counter()-t0)

        # if x is nan dist is infinite
        n_cell[~np.isfinite(dist)] = -1
        si.block_timer('Initial cell guess', t0)
        return n_cell

    def are_points_inside_domain(self,grid, xq):
        n_cell  = self.initial_horizontal_cell(grid, xq)
        bc = self.get_bc_cords(grid, xq,n_cell)
        is_inside=  np.all(np.logical_and(bc >= -self.params['bc_walk_tol'], bc  <= 1.+self.params['bc_walk_tol']),axis=1)
        return is_inside, n_cell, bc  # is inside if  magnitude of all BC < 1



    def get_bc_cords(self,grid, x,n_cells):
        # get BC cords for given x's
        si= self.shared_info
        bc_cords = np.full((x.shape[0],3), 0.)
        tri_interp_util.calc_BC_cords_numba(x, n_cells, grid['bc_transform'], bc_cords)
        return bc_cords

    def update_dry_cell_index(self,grid, current_buffer_steps, fractional_time_steps):
        # update 0-255 dry cell index
        triangle_eval_interp.update_dry_cell_index(grid['is_dry_cell'], grid['dry_cell_index'],
                                                   current_buffer_steps, fractional_time_steps)
    def close(self):
        si=self.shared_info
        info = self.info
        # transfer walk / step info from to class attributes to write in case_info file

        # convert mapped array counts to int
        for key in info.keys():
            if isinstance(info[key], np.ndarray) and info[key].size==1:
                info[key] = int(info[key])

        info['average_number_of_triangles_walked']= info['number_of_triangles_walked']/max(1,info['particles_located_by_walking'])
        if si.is3D_run:
            info['average_vertical_walk_steps'] = info['total_vertical_steps_walked'] / max(1, info['particles_located_by_walking'])

        f = f" Triangle walk summary: Of  {info['particles_located_by_walking']:6,d} particles located, "
        f += f" {info['triangle_walks_retried']:1d} walks were too long and were retried, "
        f += f" of these  {info['particles_killed_after_triangle_walk_retry_failed']:1d} failed after retrying and were discarded"
        si.msg_logger.progress_marker(f)

        if info['triangle_walks_retried'] > 100:
            si.msg_logger.msg(f" Of {info['particles_located_by_walking']:3d} particles located {info['triangle_walks_retried']:3d}, were to long and restated from initial guess",
                              crumbs='Interpolator calls, InterpTriangularNativeGrid_Slayer_and_LSCgrid',
                              hint= f"This number of retries can be reduced by decreasing the time step.")

        if info['particles_killed_after_triangle_walk_retry_failed'] > 0:
            si.msg_logger.msg(f" Of {info['particles_located_by_walking']:3d} particles located {info['particles_killed_after_triangle_walk_retry_failed']:3d}, failed to find cell",
                              crumbs='Interpolator calls, InterpTriangularNativeGrid_Slayer_and_LSCgrid',
                              hint= f"Try increasing interpolator parameter 'max_search_steps', current value ={self.params['max_search_steps']:3d}")



