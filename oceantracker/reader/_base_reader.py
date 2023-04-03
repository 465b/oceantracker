import numpy as np
from oceantracker.util.parameter_base_class import ParameterBaseClass, make_class_instance_from_params
from oceantracker.util.parameter_checking import ParamDictValueChecker as PVC, ParameterListChecker as PLC
from oceantracker.util import time_util
from os import path, walk
from glob import glob
from oceantracker.util.ncdf_util import NetCDFhandler
from time import perf_counter
from oceantracker.fields.util import fields_util
from oceantracker.util.basic_util import nopass
from oceantracker.reader.util.reader_util import append_split_cell_data
from oceantracker.util import  cord_transforms
from oceantracker.reader.util import shared_reader_memory_util


from oceantracker.reader.util import reader_util

class _BaseReader(ParameterBaseClass):

    def __init__(self):
        super().__init__()  # required in children to get parent defaults and merge with give params
        self.add_default_params({'input_dir': PVC(None, str),
                                 'file_mask': PVC(None, str, is_required=True, doc_str='Mask for file names, eg "scout*.nc", is joined with "input_dir" to give full file names'),
                                 'grid_file': PVC(None, str, doc_str='File name with hydrodynamic grid data, as path relative to input_dir, default is get grid from first hindasct file'),
                                 'coordinate_projection' : PVC(None, str, doc_str='string map project for meters grid for use by pyproj module, eg  "proj=utm +zone=16 +datum=NAD83" '),
                                 'minimum_total_water_depth': PVC(0.25, float, min=0.0,doc_str= 'Min. water depth used to decide if stranded by tide and which are dry cells to block particles from entering'),
                                 'time_zone': PVC(None, int, min=-12, max=23),
                                 'cords_in_lat_long': PVC(False, bool),
                                 'time_buffer_size': PVC(48, int, min=2),
                                 'required_file_variables' :PLC([], [str]),
                                 'required_file_dimensions': PLC([], [str]),
                                 'water_density': PVC(48, int, min=2),
                                 'depth_average': PVC(False, bool),  # turns 3D hindcast into a 2D one
                                 'field_variables_to_depth_average': PLC([], [str]),  # list of field_variables that are depth averaged on the fly
                                 'one_based_indices' :  PVC(False, bool,doc_str='indcies in hindcast start at 1, not zero, eg. triangulation nodes start at 1 not zero as in python'),
                                 'grid_variables': {'time': PVC('time', str, is_required=True),
                                                    'x': PLC(['x', 'y'], [str], fixed_len=2),
                                                    'zlevel': PVC(None, str),
                                                    'bottom_cell_index': PVC(None, str),
                                                    'is_dry_cell': PVC(None, np.int8, doc_str='Time variable flag of when cell is dry, 1= is dry cell')},
                                 'field_variables': {'water_velocity': PLC(['u', 'v', None], [str, None], fixed_len=3,is_required=True),
                                                     'water_depth': PVC(None, str),
                                                     'tide': PVC(None, str),
                                                     'water_temperature': PVC(None, str),
                                                     'salinity': PVC(None, str),
                                                     'wind_stress': PVC(None, str),
                                                     'bottom_stress': PVC(None, str),
                                                     },

                                 'dimension_map': {'time': PVC('time', str, is_required=True), 'node': PVC('node', str), 'z': PVC(None, str),
                                                   'vector2Ddim': PVC(None, str), 'vector3Ddim': PVC(None, str)},
                                 'isodate_of_hindcast_time_zero': PVC('1970-01-01', 'iso8601date'),
                                 'search_sub_dirs': PVC(False, bool),
                                 'max_numb_files_to_load': PVC(10 ** 7, int, min=1)
                                 })  # list of normal required dimensions
        self.grid = {}
        self.grid_time_buffers = {} # for time varying grid variables

        # store instances of shared memory classes for variables shared between processes
        self.shared_memory= {'grid' :{}, 'fields':{},'control':{}}

    #required read methods non time dependent variables
    def read_nodal_x_as_float64(self, nc): nopass('reader method: read_x is required')
    def read_bottom_cell_index_as_int32(self, nc):nopass('reader method: read_bottom_cell_index_as_int32 is required for 3D hindcasts')

    # required methods time dependent variables, also require a set up method
    def read_zlevel_as_float32(self, nc, file_index, zlevel_buffer, buffer_index): nopass('reader method: read_zlevel_as_float32 is required for 3D hindcasts')

    # checks on first hindcast file
    def additional_setup_and_hindcast_file_checks(self, nc,msg_logger): pass

    def make_non_time_varying_grid(self,nc, grid): nopass('setup_grid required')

    # required variable  structure query methods
    def is_var_in_file_3D(self, nc, var_name_in_file):
        return self.params['dimension_map']['z'] is not None and nc.is_var_dim(var_name_in_file, self.params['dimension_map']['z'])

    def get_num_vector_components_in_file_variable(self,nc,file_var_name):
        dm = self.params['dimension_map']
        if dm[ 'vector2Ddim'] is not None and  nc.is_var_dim(file_var_name, dm['vector2Ddim']):
            n_comp = 2
        elif dm[ 'vector3Ddim'] is not None and  nc.is_var_dim(file_var_name, dm['vector3Ddim']):
            n_comp = 3
        else:
            n_comp = 1
        return  n_comp

    def is_file_variable_time_varying(self,nc, var_name_in_file): return nc.is_var_dim(var_name_in_file, self.params['dimension_map']['time'])

    def get_number_of_z_levels(self,nc): return nc.get_dim_size(self.params['dimension_map']['z'])

    def is_hindcast3D(self, nc): nopass('must define method to test if hindcast is 3D')

    # working methods

    def preprocess_field_variable(self, nc,name, data): return data # allows tweaks to named fields, eg if name=='depth:

    def _add_grid_attributes(self, grid): pass

    def initialize(self):
        # map variable internal names to names in NETCDF file
        # set update default value and vector variables map  based on given list
        si = self.shared_info

    def get_list_of_files_and_hindcast_times(self, input_dir):
        # get list of files matching mask
        if self.params['search_sub_dirs']:
            # search hindcast sub dirs
            file_names = []
            for root, dirs, files in walk(input_dir):
                # add matching files in root folder to list
                new_files = glob(path.join(root, self.params['file_mask']))
                if len(new_files) > 0: file_names += new_files
        else:
            file_names = glob(path.normpath(path.join(input_dir, self.params['file_mask'])))

        file_info = {'names': file_names, 'n_time_steps': [], 'date_start': [], 'date_end': []}
        for n, fn in enumerate(file_names):
            # get first/second/last time from each file,
            nc = NetCDFhandler(fn, 'r')
            time = self.read_datetime(nc)
            nc.close()
            file_info['date_start'].append(time[0])
            file_info['date_end'].append(time[-1]) # -1 guards against there being only one time step in the file
            file_info['n_time_steps'].append(time.shape[0])
            if n + 1 >= self.params['max_numb_files_to_load']: break

        # check some files found
        if len(file_info['names']) == 0:
            self.msg_logger.msg('reader: cannot find any files matching mask "' + self.params['file_mask']
                           + '"  in input_dir : "' + self.params['input_dir'] + '"', fatal_error=True)
            # convert file info to numpy arrays for sorting
        keys = ['names', 'n_time_steps', 'date_start', 'date_end']
        for key in keys:
            file_info[key] = np.asarray(file_info[key])

        # sort files in time order
        file_order = np.argsort(file_info['date_start'], axis=0)
        for key in keys:
            file_info[key] = file_info[key][file_order]
        file_info['names'] = file_info['names'].tolist()

        return file_info

    def get_hindcast_files_info(self):
        # read through files to get start and finish times of each file
        # create a time sorted list of files given by file mask in file_info dictionary
        # sorts based on time from read time,  assumes a global time across all files
        # note this is only called once by OceantrackRunner to form file info list,
        # which is then passed to  OceanTrackerCaseRunner
        msg_logger = self.msg_logger
        # build a dummy non-initialise reader to get some methods and full params
        # add defaults from template, ie get reader class_name default, no warnings, but get these below
        # check cals name
        file_info = self.get_list_of_files_and_hindcast_times(self.params['input_dir'])
        info = {}
        # checks on hindcast using first hindcast file
        nc = NetCDFhandler(file_info['names'][0], 'r')
        self._basic_file_checks(nc, msg_logger)

        self.additional_setup_and_hindcast_file_checks(nc, msg_logger)

        nc.close()

        # get index at start and end on files
        cs = np.cumsum(file_info['n_time_steps'])
        file_info['nt_starts'] = cs - file_info['n_time_steps']
        file_info['n_time_steps_in_hindcast'] = np.sum(file_info['n_time_steps'], axis=0)

        # checks on hindcast
        if file_info['n_time_steps_in_hindcast'] < 2:
            msg_logger.msg('Hindcast must have at least two time steps, found ' + str(file_info['n_time_steps_in_hindcast']), fatal_error=True)

        # check for large time gaps between files
        info['first_date'] = file_info['date_start'][0]
        info['last_date'] = file_info['date_end'][-1]
        info['duration'] = info['last_date'] - info['first_date']
        info['hydro_model_timedelta'] = info['duration'] / (file_info['n_time_steps_in_hindcast'] - 1)
        info['first_time']  = np.float64(info['first_date'])
        info['last_time'] = np.float64(info['last_date'])
        info['hydro_model_time_step'] = np.float64(info['hydro_model_timedelta'])

        # check if time diff between starts of file and end of last are larger than average time step
        if len(file_info['date_start']) > 1:
            dt_gaps = file_info['date_start'][1:] - file_info['date_end'][:-1]
            sel = np.abs(dt_gaps.astype(np.float64)) > 1.8 * info['hydro_model_time_step']
            if np.any(sel):
                msg_logger.msg('Some time gaps between hindcast files is are > 1.8 times average time step, check hindcast files are all present??', hint='check hindcast files are all present and times in files consistent', warning=True)
                for n in np.flatnonzero(sel):
                    msg_logger.msg('file gaps between ' + file_info['names'][n] + ' and ' + file_info['names'][n + 1], tabs=1)

        msg_logger.exit_if_prior_errors('exiting from _get_hindcast_files_info, in setting up readers')
        info['file_info'] = file_info
        return info

    def make_grid_builder(self, grid,grid_time_buffers,reader_build_info):
        # make share memory builder for the non time varying grid variables
        #todo put reader_build_info['use_shared_memory'] test around all below
        reader_build_info['grid_constant_arrays_builder'] = {}
        for key, item in grid.items():
            if item is not None and type(item) == np.ndarray:
                if reader_build_info['use_shared_memory']:
                    sm = shared_reader_memory_util.create_shared_arrayy(values=item)
                    self.shared_memory['grid'][key] = sm  # retains a reference to keep sm alive in windows, othewise will quickly be deleted
                    reader_build_info['grid_constant_arrays_builder'][key] = sm.get_shared_mem_map()
                else:
                    reader_build_info['grid_constant_arrays_builder'][key] = {'shape': item.shape, 'dtype': item.dtype}

        # now make info to build time buffers, eg time, zlevel
        reader_build_info['grid_time_buffers_builder'] = {}
        for key, item in grid_time_buffers.items():
            if reader_build_info['use_shared_memory']:  # make shared moemory for shared_reader
                sm = shared_reader_memory_util.create_shared_arrayy(values=item)
                self.shared_memory['grid'][key] = sm  # retains a reference to keep sm alive in windows, othewise will quickly be deleted
                reader_build_info['grid_time_buffers_builder'][key] = sm.get_shared_mem_map()
            else:
                reader_build_info['grid_time_buffers_builder'][key] = {'shape': item.shape, 'dtype': item.dtype}

        return reader_build_info

    def maker_field_builder(self,nc, reader_build_info):
        # loop over reader field params
        reader_build_info['field_builder'] ={}
        for name, field_variable_comps in self.params['field_variables'].items():
            if field_variable_comps is not None:
                field_params, comp_info = self.get_field_variable_info(nc, name, field_variable_comps)
                reader_build_info['field_builder'][name] = {'field_params': field_params,'variable_info': comp_info}
                pass

        return reader_build_info


    def setup_reader_fields(self, reader_build_info):
        si = self.shared_info
        fm = si.classes['field_group_manager']
        self.code_timer.start('build_hindcast_reader')
        nc = NetCDFhandler(reader_build_info['info']['file_info']['names'][0], 'r')

        # setup reader fields from their named components in field_variables param ( water depth ad tide done earlier from grid variables)
        for name, item in reader_build_info['field_builder'].items():
            i = make_class_instance_from_params(item['field_params'], si.msg_logger)
            i.info['variable_info']= item['variable_info']
            si.add_class_instance_to_interator_lists('fields', 'from_reader_field', i, crumbs='Adding Reader Field "' + name + '"')
            i.initialize()  # require variable_info to initialise

            if reader_build_info['use_shared_memory']:
                #todo make this part of reader field intialize???
                self.shared_memory['fields'][name] = shared_reader_memory_util.create_shared_arrayy(values=i.data)

            if not i.params['is_time_varying']:
                # if not time dependent field read in now,
                data = self.assemble_field_components(nc, i)
                data = self.preprocess_field_variable(nc,name, data)
                i.data[:] = data

            # set up depth averaged version if requested
            if name in self.params['field_variables_to_depth_average']:
                # tweak shape to fit depth average of scalar or 3D vector
                p = {'class_name':'oceantracker.fields.reader_field.DepthAveragedReaderField',
                     'name': name + '_depth_average','num_components': min(2, i.params['num_components']),
                     'is_time_varying': i.params['is_time_varying'],
                    'is3D': False}
                i2 = make_class_instance_from_params(p,si.msg_logger)
                si.add_class_instance_to_interator_lists('fields', 'depth_averaged_from_reader_field', i2,
                                                         crumbs='Adding Reader Depth Averaged Field "' + name + '"')
                i2.initialize()
        nc.close()

        # rinf buffer ono, needed to force read at first time step read to make
        bi = self.info['buffer_info']
        bi['n_filled'] = 0
        bi['buffer_size'] = self.params['time_buffer_size']
        bi['buffer_available'] = bi['buffer_size']
        bi['nt_buffer0'] = 0

        self.code_timer.stop('build_hindcast_reader')

    def _basic_file_checks(self, nc, msg_logger):
        # check named variables are in first file
        # check dim
        for name, d in self.params['dimension_map'].items():
            if d is not None and not nc.is_dim(d):
                msg_logger.msg('Cannot find dimension_map dimension "' + name + ' ", file dimension given is "' + d + '"',
                               hint='Dimensions in hydro-model file = ' + str(nc.get_dims()),
                               fatal_error=True)

        # check variables are there
        for vm in ['grid_variables', 'field_variables']:
            for name, d in self.params[vm].items():
                if type(d)== list:
                    for vf in d:
                        if vf is not None and not nc.is_var(vf):
                            msg_logger.msg(' For  "' + vm + '" for param   "' + name + ' ",  cannot find variable in file  "' + vf + '"', fatal_error=True)

                elif d is not None and not nc.is_var(d) :
                    msg_logger.msg('For "' + vm + '" for param,  "' + name + ' ", cannot find variable in file "' + str(d) + '"', fatal_error=True)


        # check if all required dims and non-feilds variables present
        for v in self.params[ 'required_file_variables']:
            if not nc.is_var(v):
                msg_logger.msg('Cannot find required variable in hydro model output file "' + v + '"', fatal_error=True)

        for v in self.params[ 'required_file_dimensions']:
            if not nc.is_dim(v):
                msg_logger.msg( 'Cannot find required dimension in hydro model outptut file "' + v + '"', fatal_error=True)




    def get_field_variable_info(self, nc, name,var_list):
        # get info from list of component eg ['temp'], ['u','v']
        si= self.shared_info

        if type(var_list) is not list: var_list=[var_list] # if a string make a list of 1
        var_list = [v for v in var_list if v != None]
        var_file_name0=var_list[0]

        is3D_in_file = self.is_var_in_file_3D(nc, var_file_name0)
        is_time_varying = self.is_file_variable_time_varying(nc, var_file_name0)

        var_info= { 'component_list':[],
                    'is3D_in_file': is3D_in_file,
                    'is_time_varying': is_time_varying,
                    'requires_depth_averaging': self.params['depth_average'] and is3D_in_file}

        # work out number of components in list of variables
        n_total_comp=0
        for file_var_name in var_list:
            n_comp =self.get_num_vector_components_in_file_variable(nc,file_var_name)
            n_total_comp += n_comp
            var_info['component_list'].append({'name_in_file':file_var_name, 'num_components': n_comp })

        # if a 3D var and vector then it must have  3D components
        # eg this allows for missing vertical velocity, to have zeros in water_velocity
        if is3D_in_file and n_total_comp> 1: n_total_comp =3

        params = {'name': name,
                  'class_name':'oceantracker.fields.reader_field.ReaderField',
             'is_time_varying':  is_time_varying,
             'num_components' : n_total_comp,
             'is3D' :  False if var_info['requires_depth_averaging'] else is3D_in_file
            }

        return params, var_info
    


    def fill_time_buffer(self, nt0_hindcast):
        # fill as much of  hindcast buffer as possible starting at global hindcast time step nt0_buffer
        # fill buffer starting at hindcast time step nt0_buffer
        # todo change so does not read current step again after first fill of buffer

        self.code_timer.start('reading_to_fill_time_buffer')
        t0 = perf_counter()

        si = self.shared_info
        grid = self.grid
        grid_time_buffers = self.grid_time_buffers

        fi = self.info['file_info']
        bi = self.info['buffer_info']
        buffer_size = bi['buffer_size']


        # get hindcast global time indices of first block, loads in model order
        # ie if backtracking are still moving forward in buffer
        total_read = 0
        bi['buffer_available'] = buffer_size
        bi['nt_buffer0'] = nt0_hindcast

        n_file = np.argmax(nt0_hindcast < fi['nt_starts'] + fi['n_time_steps']) # find first file with time step

        # always move forward in buffer, by loading hindcast backwards if backtracking
          # from current step to end of hindcast
         # buffer is always moving forward
        if si.backtracking:
            # from start to current step available,by loading time steps in reverse order
            nt_total_available = min(nt0_hindcast, bi['buffer_available'])
            nt_hindcast_required = nt0_hindcast - np.arange(nt_total_available)
            nt_file_required = fi['n_time_steps_in_hindcast']- 1 - nt_hindcast_required
        else:
            # forward tracking, time steps same as hindcast time steps, test
            nt_total_available = min(fi['n_time_steps_in_hindcast'] - nt0_hindcast, bi['buffer_available'])
            nt_file_required = nt_hindcast_required.copy()

        while len(nt_file_required) > 0 and 0 <= n_file < len(fi['names']):

            t0_file = perf_counter()
            nc = NetCDFhandler(fi['names'][n_file], 'r')

            if si.backtracking:
                    # from first remaining requested step to beginning of file is avaiable
                    nt_available = nt_file_required[0] -  fi['nt_starts'][n_file]
            else:
                # from first requested file time step to end of file is available
                nt_available = fi['nt_starts'][n_file]+ fi['n_time_steps'][n_file] - nt_file_required[0]

            num_read = min(nt_available, nt_file_required.size, bi['buffer_available'])
            file_offsets = nt_file_required[:num_read] - fi['nt_starts'][n_file] # steps in current file
            buffer_index = nt_hindcast_required[:num_read]  % bi['buffer_size']     # make ring buffer index in by modula maths

            s =  f'Reading-file-{(n_file+1):02d}  {path.basename(fi["names"][n_file])}, steps in file {fi["n_time_steps"][n_file]:3d},'
            s += f'reading  {num_read:2d} of {nt_total_available:2d} steps, '
            s += f' for hindcast time steps {nt_hindcast_required[0]:02d}:{nt_hindcast_required[-1]:02d}, '
            s += f' available {fi["nt_starts"][n_file]:03d}:{fi["nt_starts"][n_file]+fi["n_time_steps"][n_file]-1:03d}, '
            s += f' read hindcast indices {nt_hindcast_required[0]:03d}:{nt_hindcast_required[num_read-1]:03d}'
            s += f' file  offsets  {file_offsets[0]:03d}:{file_offsets[num_read - 1]:03d} to ring buffer offsets {buffer_index[0]:03}:{buffer_index[num_read - 1]:03d} '
            si.msg_logger.write_progress_marker(s)

            grid_time_buffers['nt_hindcast'][buffer_index] = nt_file_required[:num_read]  # add a grid variable with global hindcast time steps

            # read time varying vector and scalar reader fields
            for name, field in si.class_interators_using_name['fields']['from_reader_field'].items():
                if field.is_time_varying():
                    data_added_to_buffer = self.assemble_field_components(nc, field, buffer_index=buffer_index, file_index=file_offsets)
                    data_added_to_buffer = self.preprocess_field_variable(nc, name, data_added_to_buffer) # do any customised tweaks

                    field.data[buffer_index, ...] = data_added_to_buffer

                    if name in self.params['field_variables_to_depth_average']:
                       si.classes['fields'][name + '_depth_average'].data[buffer_index, ...] = fields_util.depth_aver_SlayerLSC_in4D(data_added_to_buffer, grid_time_buffers['zlevel'], grid['bottom_cell_index'])

            # read grid time, zlevel
            # do this after reading fields as some hindcasts required tide field to get zlevel, eg FVCOM
            self.read_time_variable_grid_variables(nc, buffer_index,file_offsets)

            nc.close()

            # now all  data has been read from file, now
            # update user fields from newly read fields and data
            for field_types in ['derived_from_reader_field','user']:
                for field in si.class_interators_using_name['fields'][field_types].values():
                    if field.is_time_varying():
                        field.update(buffer_index)

            total_read += num_read

            # set up for next step
            bi['buffer_available'] -= num_read
            n_file += int(si.model_direction)
            nt_file_required = nt_file_required[num_read:]
            nt_hindcast_required = nt_hindcast_required[num_read:]

        si.msg_logger.write_progress_marker( f' read {total_read:3d} time steps in  {perf_counter()-t0:3.1f} sec',tabs=2)


        # record useful info/diagnostics


        bi['n_filled'] = total_read


        self.code_timer.stop('reading_to_fill_time_buffer')

    def assemble_field_components(self,nc, field, buffer_index=None, file_index=None):
        # read scalar fields / join together the components which make vector from component list

        grid = self.grid
        grid_time_buffers = self.grid_time_buffers

        m= 0 # num of vector components read so far
        var_info = field.info['variable_info']
        for component_info in var_info['component_list']:
            data = self.read_file_field_variable_as4D(nc, component_info,var_info['is_time_varying'], file_index)

            if var_info['requires_depth_averaging']:
                data = fields_util.depth_aver_SlayerLSC_in4D(data, grid_time_buffers['zlevel'], grid['bottom_cell_index'])

            m1 = m + component_info['num_components']

            # get view of where in buffer data is to be placed
            if var_info['is_time_varying']:
                field.data[buffer_index, :, :, m:m1] = data
            else:
                field.data[0, :, :, m:m1] = data

            m += component_info['num_components']

        # return a view of data added to buffer to allow pre-processing
        data_added_to_buffer= field.data[buffer_index, ...] if field.params['is_time_varying'] else field.data[0,...]
        return data_added_to_buffer

    def read_file_field_variable_as4D(self, nc, file_var_info,is_time_varying, file_index=None):
        # reformat file variable into 4D time,node,depth, components  form
        var_name= file_var_info['name_in_file']

        data = nc.read_a_variable(var_name, sel= file_index if is_time_varying else None).astype(np.float32) # allow for time independent data

        # reorder dim to time,node,depth, components order, if present
        # default is in correct order

        # now reshape in 4D
        if not self.is_file_variable_time_varying(nc,var_name): data = data[np.newaxis,...]
        if not self.is_var_in_file_3D(nc, var_name):    data = data[:, :, np.newaxis,...]
        if file_var_info['num_components'] == 1:             data = data[:, :, :, np.newaxis]

        return data


    def read_dry_cell_data(self,nc,file_index,is_dry_cell_buffer, buffer_index):
        # calculate dry cell flags, if any cell node is dry
        si = self.shared_info
        grid = self.grid
        grid_time_buffers = self.grid_time_buffers
        fields = si.classes['fields']

        if self.params['grid_variables']['is_dry_cell'] is None:
            if grid_time_buffers['zlevel'] is None:
                reader_util.set_dry_cell_flag_from_tide( grid['triangles'],
                                                        fields['tide'].data, fields['water_depth'].data,
                                                        si.minimum_total_water_depth, is_dry_cell_buffer,buffer_index)
            else:
                reader_util.set_dry_cell_flag_from_zlevel( grid['triangles'],
                                                          grid_time_buffers['zlevel'], grid['bottom_cell_index'],
                                                          si.minimum_total_water_depth, is_dry_cell_buffer,buffer_index)
        else:
            # get dry cells for each triangle allowing for splitting quad cells
            data_added_to_buffer = nc.read_a_variable(self.params['grid_variables']['is_dry_cell'], file_index)
            is_dry_cell_buffer[buffer_index, :] = append_split_cell_data(grid, data_added_to_buffer, axis=1)
            #grid['is_dry_cell'][buffer_index, :] =  np.concatenate((data_added_to_buffer, data_added_to_buffer[:, grid['quad_cell_to_split_index']]), axis=1)

    def read_open_boundary_data_as_boolean(self, grid):
        is_open_boundary_node = np.full((grid['x'].shape[0],), False)
        return is_open_boundary_node


    def time_to_hydro_model_time_step(self, time_sec):
        #convert date time to global time step in hindcast just before/after when forward/backtracking
        # always move forward through buffer, but file info is always forward in time
        si = self.shared_info
        info = self.info

        hindcast_fraction= (time_sec - info['first_time']) / (info['last_time'] - self.info['first_time'])
        nt = (info['file_info']['n_time_steps_in_hindcast'] - 1) *  hindcast_fraction

        return np.int32(np.floor(nt))


    def get_buffer_index_from_hindcast_global_time_step(self, nt):
        # ring buffer mapping
        return nt % self.info['buffer_info']['buffer_size']

    def time_to_buffer_index(self, time_sec):
        nt = self.time_to_hydro_model_time_step(time_sec)
        return self.get_buffer_index_from_hindcast_global_time_step(nt)

    def are_time_steps_in_buffer(self, time_sec):
        # check if next two steps of remaining  hindcast time steps required to run  are in the buffer
        si = self.shared_info
        bi = self.info['buffer_info']


        # find first and last times in ring buffer
        nt_buffer = self.time_to_hydro_model_time_step(time_sec)
        out = bi['nt_buffer0'] <= nt_buffer <  (bi['nt_buffer0'] + bi['n_filled'] - 1) # ensure this and the next are in the buffer
        print('xx', bi['nt_buffer0'] , nt_buffer ,  (bi['nt_buffer0'] + bi['buffer_size'] - 1),out)

        return out

    def _open_grid_file(self,reader_build_info):
        if self.params['grid_file']:
            file_name = path.join(self.params['input_dir'],self.params['grid_file'])
        else:
            file_name= reader_build_info['info']['file_info']['names'][0]

        nc =NetCDFhandler(file_name, 'r')
        return nc

    def set_up_shared_grid_memory(self, reader_build_info):
        if 'shared_memory' not in reader_build_info:
            # build shared memory and add to reader_build_info
            reader_build_info['shared_memory'] = {'grid': {},'fields':{}}
            for key, item in self.grid.items():
                if item is not None:
                    sm = shared_memory.SharedMemArray(values=item)
                    self.shared_memory['grid'][key] = sm
                    reader_build_info['shared_memory']['grid'][key] = sm.get_shared_mem_map()
        else:
            # build grid variables from reader_build_info shared_memory info
            for key, item in reader_build_info['shared_memory']['grid'].items():
                self.shared_memory['grid'][key] = shared_memory.SharedMemArray(sm_map=item)
                self.grid[key] = self.shared_memory['grid'][key].data  # grid variables is shared version
        return reader_build_info

    def convert_lon_lat_to_meters_grid(self, x):

        if self.params['coordinate_projection'] is None:
            x_out, self.cord_transformer= cord_transforms.WGS84_to_UTM( x, out=None)
        else:
            #todo make it work with users transform?
            x_out = cord_transforms.WGS84_to_UTM(x, out=None)
        return x_out

    def close(self):
        # release any shared memory
        sm_info=self.shared_memory
        for sm in list(sm_info['grid'].values())+ list(sm_info['grid'].values()):
            sm.delete()