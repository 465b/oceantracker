from oceantracker.particle_statistics._base_location_stats import _BaseParticleLocationStats
from oceantracker.util.parameter_checking import  ParamDictValueChecker as PVC, ParameterListChecker as PLC
from oceantracker.util.message_and_error_logging import GracefulExitError,FatalError

from oceantracker.particle_release_groups.polygon_release import PolygonRelease
from oceantracker.util.polygon_util import InsidePolygon, inside_ray_tracing_single_point
import numpy as np
from numba import njit

class ResidentInPolygon(_BaseParticleLocationStats):
    def __init__(self):
        # set up info/attributes
        super().__init__()

        self.add_default_params({'count_release_group':  PVC(None, int, min=1,is_required=True,
                                     doc_str='Numer of polygon,release group to count particles inside, 1-N'),
                                 'role_output_file_tag': PVC('residence', str),
                                 'z_range': PLC([], [float, int], min_length=2, doc_str='z range = [zmin, zmax] count particles in this z range in 3D'),
                                 })

    def initialize(self, **kwargs):
        si= self.shared_info
        params = self.params
        # do standard stats initialize
        super().initialize()  # set up using regular grid for  stats

        # get associated release group
        if params['count_release_group'] > len(si.classes['particle_release_groups']):
            si.case_log.write_msg('Parameter count_release_group be  range 1- number of release groups, given  value=' + str(params['count_release_group'], exception=GracefulExitError))

        self.release_group = si.classes_as_lists['particle_release_groups'][params['count_release_group'] - 1]

        if not isinstance( self.release_group, PolygonRelease):
            #todo fixed
            si.case_log.write_warning(
                'Parameter count_release_group be  range 1- number of release groups, given  value='
                    + str(params['count_release_group']), exception=GracefulExitError)

        self.polygon = InsidePolygon(self.release_group.params['points'])

        # tag file with release group number
        #params['role_output_file_tag'] += '_RG%3.0f ' % params['count_release_group']
        self.open_output_file()

        self.set_up_time_bins(self.nc)
        self.set_up_binned_variables(self.nc)
        self.set_up_part_prop_lists()

        #todo move to base location stats, so all can use depth range
        if len(self.params['z_range'])==0:
            self.params['z_range'] = [-1.0e30, 1.0e30]

        self.params['z_range']= np.asarray(self.params['z_range'])

    def check_requirements(self):
        si= self.shared_info
        params = self.params
        msg_list = self.check_class_required_fields_properties_grid_vars_and_3D()
        return msg_list

    def set_up_binned_variables(self, nc):
        si = self.shared_info
        if not self.params['write']: return
        rg_info = self.release_group.info['release_info']

        dim_names = ('time', 'pulse')
        num_pulses= len(rg_info['release_schedule_times'])
        nc.add_a_Dimension('pulse', dim_size=num_pulses)
        nc.create_a_variable('count', dim_names,
                             {'notes': 'counts of particles in each pulse of release group inside release polygon at given times'},
                             np.int32)
        nc.create_a_variable('count_all_particles', ['time', 'pulse'],
                             {'notes': 'counts of particles in each, whether inside polygon or not at given times'}, np.int32)
        # set up space for requested particle properties
        # working count space
        self.count_time_slice = np.full((num_pulses,), 0, np.int32)
        self.count_all_particles_time_slice = np.full_like(self.count_time_slice, 0, np.int32)

        for p_name in self.params['particle_property_list']:
            if p_name in si.classes['particle_properties']:
                self.sum_binned_part_prop[p_name] = np.full(num_pulses, 0.)  # zero for  summing
                nc.create_a_variable('sum_' + p_name, dim_names,
                                     {'notes': 'sum of particle property inside polygon  ' + p_name}, np.float64)
            else:
                si.case_log.write_msg(
                    'Part Prop "' + p_name + '" not a particle property, ignored and no stats calculated')

    def set_up_spatial_bins(self,nc ): pass

    def update(self,**kwargs):
        si= self.shared_info
        part_prop = si.classes['particle_properties']
        rg  = self.release_group
        poly= self.polygon

        # update time stats  recorded
        time = kwargs['time']
        self.record_time_stats_last_recorded(time)

        sel = self.select_particles_to_count(self.get_particle_index_buffer())

        # do counts
        self.do_counts_and_summing_numba(poly.line_bounds, poly.slope_inv, poly.polygon_bounds,
                                    part_prop['IDrelease_group'].data,
                                    part_prop['IDpulse'].data,
                                    rg.info['instanceID'],
                                    self.params['z_range'],
                                    part_prop['x'].data,
                                    self.count_time_slice, self.count_all_particles_time_slice,
                                    self.prop_list, self.sum_prop_list, sel)

        self.write_time_varying_stats(self.nWrites,time)
        self.nWrites += 1

    def info_to_write_at_end(self):pass  # nothing extra to write

    @staticmethod
    @njit
    def do_counts_and_summing_numba(lb, slope_inv, bounds,
                                    release_group_ID, pulse_ID, required_release_group,zrange, x, count,
                                    count_all_particles, prop_list, sum_prop_list, active):
        # count those of each pulse inside release polygon

        # zero out counts in the count time slices
        count[:] = 0
        count_all_particles[:] = 0
        for m in range(len(prop_list)):
            sum_prop_list[m][:] = 0.

        for n in active:
            if  release_group_ID[n] == required_release_group:
                pulse= pulse_ID[n]

                # count all particles not whether in volume or not
                count_all_particles[pulse] += 1  # all particles count whether in a polygon or not , whether in required depth range or not

                if x.shape[1] == 3 and not (zrange[0] <= x[n, 2] <= zrange[1]): continue

                if inside_ray_tracing_single_point(x[n,:], lb, slope_inv, bounds):

                    count[pulse] += 1
                    # sum particle properties
                    for m in range(len(prop_list)):
                        sum_prop_list[m][pulse] += prop_list[m][n]