from oceantracker.util.parameter_base_class import ParameterBaseClass
from oceantracker.util.parameter_checking import ParamValueChecker as PVC, ParameterListChecker as PLC
from oceantracker.util.basic_util import nopass
from oceantracker.common_info_default_param_dict_templates import particle_info
from oceantracker.util.ncdf_util import NetCDFhandler
from oceantracker.util import output_util
from os import path
from datetime import datetime
import numpy as np

class _BaseTriangleProperties(ParameterBaseClass):

    def __init__(self):
        super().__init__()
        # set up info/attributes
        self.add_default_params({'class_name' : PVC(None, str,is_required=True),
                                 'particle_properties_to_track': PLC([],[str],  make_list_unique=True),
                                 'write': PVC(True, bool),
                                 'role_output_file_tag': PVC('_concentrations_', str),
                                 'count_status_equal_to': PVC(None, str, possible_values=particle_info['status_flags'].keys()),
                                 'only_update_concentrations_on_write': PVC(True, bool),
                                 'update_interval': PVC(3600., int, min=1, units='sec', doc_str='the time in model seconds between writes (will be rounded to model time step)'),
                                 'update_values_every_time_step': PVC(False, bool, min=1, units='sec', doc_str='update values in memory every time step, needed if using concentrations within modelling to change particle behaviour or properties. Output interval still sep by update_interval')
                                 })
        self.info['time_last_stats_recorded'] = -np.inf

    def initial_setup(self):
        # set up data buffer and output variables
        self.set_up_data_buffers()
        self.set_up_output_file()

    def set_up_data_buffers(self): nopass('')

    def set_up_output_file(self):
        # set up output file
        si = self.shared_info
        grid = si.classes['field_group_manager'].grid

        self.info['output_file'] = si.output_file_base + '_' + self.params['role_output_file_tag'] + '_' + self.info['name'] + '.nc'
        si.msg_logger.progress_marker('opening concentrations output to : ' + self.info['output_file'])

        self.nc = NetCDFhandler(path.join(si.run_output_dir, self.info['output_file']), 'w')
        nc = self.nc
        nc.write_global_attribute('created', str(datetime.now().isoformat()))
        nc.add_dimension('time_dim', None)
        nc.add_dimension('triangle_dim', grid['triangles'].shape[0])

        nc.create_a_variable('time', ['time_dim'], np.float64,description='time in seconds since 1970-01-01')
        self.time_steps_written = 0
        # need to add other variables in children
        self.info['time_last_stats_recorded'] = si.time_of_nominal_first_occurrence

    def select_particles_to_count(self):
        si= self.shared_info
        part_prop =  si.classes['particle_properties']
        return part_prop['status'].compare_all_to_a_value('gteq', si.particle_status_flags['frozen'], out=self.get_partID_buffer('B1'))

    def write(self, time_sec): nopass()


    def update(self, time_sec):nopass()

    def check_requirements(self): pass

    def record_time_stats_last_recorded(self, t): self.info['time_last_stats_recorded'] = t

    def close(self):

        if hasattr(self,'nc'):
            si = self.shared_info
            nc = self.nc

            # add attributes mapping release index to release group name
            output_util.add_release_group_ID_info_to_netCDF(nc, si.classes['release_groups'])

            nc.close()






