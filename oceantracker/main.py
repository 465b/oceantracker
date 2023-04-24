# method to run ocean tracker from parameters
# eg run(params)
import sys

code_version = '0.4.00.001 2023-04-11'

# todo kernal/numba based RK4 step
# todo short name map requires unique class names in package, this is checked on startup,add checks of uniqueness of user classes added from outside package

# Dev notes
# line debug?? python3.6 -m pyinstrument --show-all plasticsTrackOnLine_Main.py
# python -m cProfile
# python -m vmprof  <program.py> <program parameters>
# python -m cProfile -s cumtime

# do first to ensure its right
import multiprocessing

from copy import deepcopy
from datetime import datetime

from os import path, makedirs
from sys import version, version_info
import shutil
from time import perf_counter
from copy import  copy
import numpy as np

from oceantracker.util.ncdf_util import NetCDFhandler
from oceantracker.util import basic_util
from oceantracker.util import json_util ,yaml_util

from oceantracker.util.parameter_checking import merge_params_with_defaults
from oceantracker.oceantracker_case_runner import OceanTrackerCaseRunner

from oceantracker import common_info_default_param_dict_templates as common_info

from oceantracker.util.parameter_util import make_class_instance_from_params
from oceantracker.util.messgage_logger import GracefulError, MessageLogger

import subprocess
import traceback


def run(raw_params):
    if type(raw_params) is not dict:
        raise GracefulError('Parameter must be a dictionary or json/yaml file readable as a dictionary, given parameters are type=' + str(type(raw_params)),
                            hint='check parameter file or parameter variable is in dictionary form')

    OT = _RunOceanTrackerClass()
    if 'shared_params' in raw_params:
        raw_params = convert_ver03_to_04_params(raw_params)

    full_runInfoJSON_file_name, has_errors = OT.run(raw_params)

    return  full_runInfoJSON_file_name, has_errors


def decompose_params(params, msg_logger, crumbs=''):
    if type(params) != dict:
        msg_logger.msg(f'Unpacking given parameters, params must be a dictionary , got type ={str(type(params))}', fatal_error=True,
                       crumbs=crumbs, exit_now=True)
    if 'case_list' in params and type(params['case_list']) != list:
        msg_logger.msg(f'Unpacking given parameters, "case_list"" must be a list , got type ={str(type(params["case_list"]) )}', fatal_error=True,
                       crumbs=crumbs, exit_now=True)

    param_parts = {'settings': {}, 'core_classes': {}, 'class_lists': {}}
    for key in params.keys():
        if key in common_info.shared_params:
            param_parts['settings'][key] = params[key]

        elif key in common_info.run_params:
                param_parts['settings'][key] = params[key]

        elif key.removesuffix('_class') in common_info.core_classes:
            param_parts['core_classes'][key.removesuffix('_class')] = params[key]

        elif key.removesuffix('_list') in common_info.class_lists:
            param_parts['class_lists'][key.removesuffix('_list')] = params[key]

        else:
            msg_logger.msg('Unexpected parameter key "' + key, warning=True, crumbs=crumbs)
    return param_parts

class _RunOceanTrackerClass(object):

    def run(self, raw_params):
        self.msg_logger = MessageLogger('M:')
        ml = self.msg_logger
        ml.insert_screen_line()
        ml.msg(common_info.package_fancy_name + ' preliminary setup')

        # keep clean parameter copy to write out
        params = deepcopy(raw_params)

        # split into run_prams and case list
        if 'case_list' in params:
            cl = params['case_list']
            params.pop('case_list')
            # unpack case list params
            case_list_params =[]
            for c in cl:
                case_list_params.append(decompose_params(c,msg_logger=ml))
        else:
            case_list_params = []

        # split  top level params into useful blocks
        run_builder = decompose_params(params, ml, crumbs='top level params check-')


        # merge top level run_params params with defaults
        run_builder['settings'] = merge_params_with_defaults(run_builder['settings'],
                                            common_info.run_params, {},
                                            ml, crumbs='merging top level parameters with defaults')
        ml.exit_if_prior_errors('Errors in top level settings')

        ml.max_warnings = run_builder['settings']['advanced_settings']['max_warnings']

        # prelim set up to get output dir, to make message logger file
        run_builder =  self.setup_output_folders_and_msg_loger_file(run_builder, raw_params)

        # complete the setup with message logger
        ml.insert_screen_line()
        ml.msg('Running ' + common_info.package_fancy_name + ' started ' + str(datetime.now()))
        ml.progress_marker('Starting: ' + run_builder['settings']['output_file_base'])
        run_builder = self.setup_check_python_version(run_builder, ml)

        # merge core classes

        # merge class names only for core classses
        #todo check for unknown top level keys?
        run_builder['core_classes']= self.check_top_level_class_keys(run_builder['core_classes'],
                                                common_info.core_classes, ml, required_keys=['reader'], crumbs='checking top level core classes')

        # get info to build a reader, and create a dummy reader instance
        run_builder, reader = self._setup_reader_builder(run_builder, ml)
        ml.progress_marker('Input directory: ' + reader.params['input_dir'])

        run_builder['is_3D_run'] = run_builder['reader_build_info']['hindcast_is3D'] and not run_builder['settings']['run_as_depth_averaged']

    # merge all param with defaults and make full case_builders
        ml.progress_marker('Setting up params for each case ')
        # make full list of cases to run with fully merged params
        case_params_list = self.setup_build_full_case_params(run_builder,case_list_params, ml)

        # now run cases
        full_runInfoJSON_file_name, has_errors = self._run(run_builder, case_params_list,reader)

        ml.close()
        return full_runInfoJSON_file_name, has_errors

    def setup_output_folders_and_msg_loger_file(self, run_builder, raw_params):
        #setus up params, opens log files/ error handling, required befor mesage loger can be used

        run_params = run_builder['settings']
        msg_logger= self.msg_logger

        # get output files location
        root_output_dir = path.abspath(path.normpath(run_params['root_output_dir']))
        run_output_dir = path.join(root_output_dir, run_params['output_file_base'])

        if run_params['add_date_to_run_output_dir']:
            run_output_dir += datetime.now().strftime("_%Y-%m-%d_%H-%M")

        # kill existing folder
        if path.isdir(run_output_dir):  shutil.rmtree(run_output_dir)

        try:
            makedirs(run_output_dir)  # make  and clear out dir for output
        except OSError as e:
            # path may already exist, but if not through other error, exit
            msg_logger.msg(f'Failed to make run output dir:{run_output_dir}',
                           exception = e, traceback_str=traceback.print_exc())

        # make message/error logger for main code  with log and error output files
        runLog_file, run_error_file=  msg_logger.set_up_files(run_output_dir, run_params['output_file_base'])

        # write a copy of user given parameters, to help with debugging and code support
        fb = 'users_params_' +  run_params['output_file_base']
        json_util.write_JSON(path.join(run_output_dir,fb), raw_params)

        run_builder['output_files'] = {'root_output_dir': root_output_dir,
                        'run_output_dir' :  run_output_dir,
                        'output_file_base': run_params['output_file_base'],
                        'runInfo_file'  : run_params['output_file_base'] + '_runInfo.json',
                        'runLog_file':  run_params['output_file_base']+ '_runScreen.log',
                        'run_error_file':  run_params['output_file_base']+ '_run.err',
                        'users_params_json' : fb + '.json',
                         'runLog_file': runLog_file,
                         'run_error_file' : run_error_file
                }

        return run_builder



    def setup_check_python_version(self,run_builder,msg_logger):
        # set up log files for run

        msg_logger.insert_screen_line()
        msg_logger.msg('Starting ' + common_info.package_fancy_name + '  Version ' + code_version)
        msg_logger.msg('Python version: ' + version, tabs=1)

        vi = version_info
        install_hint = 'Install Python 3.10 or used environment.yml to build a Conda virtual environment named oceantracker'

        # todo check share mem version
        if False and run_builder['working_params']['shared_params']['shared_reader_memory']:
            # for shared reader python version must be >=3.9
            if not (vi.major == 3 and vi.major >= 9):
                msg_logger.msg('To use shared reader memory ' +
                             common_info.package_fancy_name + ' requires Python 3 , version >= 3.9',
                             hint=install_hint, warning=True, tabs=1)

        if (vi.major == 3 and vi.major >= 11):
            msg_logger.msg(common_info.package_fancy_name + ' is not yet compatible with Python 3.11, as not al imported packages have been updated, eg Numba ',
                         hint=install_hint, fatal_error=True, exit_now= True,tabs=1)
        if vi.major < 3:
            msg_logger.msg(common_info.package_fancy_name + ' requires Python version 3 ', hint=install_hint,  fatal_error = True,exit_now=True, tabs = 1)

        msg_logger.exit_if_prior_errors()

        return run_builder

    def _run(self,run_builder, case_builders_list, reader):

        d0 = datetime.now()
        t0 = perf_counter()
        run_info = {'user_note': {}, 'screen_log': [],
                    'run_started': str(d0),
                    'run_ended': None,
                    'elasped_time': None,
                    'performance': {},
                    'output_files': {},
                    }


        output_files= run_builder['output_files']

        ml = self.msg_logger

        ml.insert_screen_line()

        # run the cases, return list of case info json files which make up the run of all cases
        #----------------------------------------------------------------------------------------------
        output_files['case_info'], case_error = self._run_case_list(case_builders_list, reader)
        # ----------------------------------------------------------------------------------------------

        # tidy up run
        tnow= datetime.now()
        run_info.update({
            'user_note': run_builder['settings']['user_note'],
            'run_ended':tnow ,
            'elasped_time': str(tnow-d0),
            'code_version_info': self._U2_code_version_info(),
            'computer': basic_util.get_computer_info(),
            'performance' : self._U1_get_all_cases_performance_info(output_files['case_info'], case_error, run_builder['settings'],
                                                                    output_files['run_output_dir'], t0),
            'output_files': output_files
                })
        # write runInfo JSON
        full_runInfoJSON_file_name = path.join(output_files['run_output_dir'], output_files['runInfo_file'])
        json_util.write_JSON(full_runInfoJSON_file_name, run_info)

        # the end
        ml.show_all_warnings_and_errors()
        ml.insert_screen_line()
        ml.progress_marker('Finished ' + '-  started: ' + str(d0) + '- ended: ' + str(datetime.now()))
        ml.msg('Elapsed time =' + str(datetime.now() - d0), tabs=3)
        ml.msg('Output in ' + output_files['run_output_dir'], tabs=4)
        ml.insert_screen_line()
        ml.close()

        has_errors= any(case_error)
        return full_runInfoJSON_file_name, has_errors

    def setup_build_full_case_params(self, run_builder,case_list_params, msg_logger):
        # make set of case params merged with defaults and checked

        # fill in any missing class lists
        run_builder['class_lists'] = self.check_top_level_class_keys(run_builder['class_lists'],
                                                  common_info.class_lists, msg_logger,
                                                        crumbs='checking class_lists')

        # merge params of base case class lists
        for key, item in  run_builder['class_lists'].items():
            for n in range(len(item)):
                i = make_class_instance_from_params(run_builder['class_lists'][key][n], msg_logger, class_type_name=key,
                                      crumbs=' merging run_params of top level class_list "' + key + '" ')
                run_builder['class_lists'][key][n] = i.params

        if len(case_list_params) ==0 :
            case_list_params=[decompose_params({}, msg_logger)] # ensure there is at least one empty case

        # build each case params
        processor_number = 0
        full_list_case_params=[]

        for n_case, case in enumerate(case_list_params):
            tag = 'case #' +str(n_case+1)
            # parameter list to run same particle on multiple threads
            c  = deepcopy(run_builder)  # a set of parameters for this case

            # fill in missing top level values
            c['core_classes'] = self.check_top_level_class_keys(c['core_classes'],
                                                               common_info.core_classes, msg_logger,
                                                               crumbs='checking class_lists  ' + tag)
            c['class_lists'] = self.check_top_level_class_keys(c['class_lists'],
                                                    common_info.class_lists, msg_logger,
                                                     crumbs='checking class_lists  ' + tag)

            # check if trying to change param which are shared by all runs
            for key in common_info.shared_params:
                if key in case['settings']:
                    msg_logger.msg(f'Cannot set parameter {key} within a case, it must be the same for all cases, ignoring this value' ,
                                   warning=True, hint ='move this param to top level or delete  from ' + tag)
                    case['settings'].pop(key)

            # merge params
            c['settings']= merge_params_with_defaults(case['settings'], common_info.run_params,
                                       c['settings'],  msg_logger,
                                       crumbs='merging run_params into case ' + tag)
            # merge core class params
            for key, params in case['core_classes'].items():
                i = make_class_instance_from_params(params, msg_logger, class_type_name=key,
                                                    crumbs=' merging run_params of core_class "' + key + '" '  + tag,
                                                    base_case_params=run_builder['core_classes'][key])
                c['core_classes'][key] = i.params

            # merge case's  class_list params and append to top level class_list
            for key, item in case['class_lists'].items():
                for n in range(len(item)):
                    i = make_class_instance_from_params(case['class_lists'][key][n], msg_logger,
                                                        class_type_name=key,
                                                        crumbs=' merging run_params of top level class_list "' + key + '" ')
                    c['class_lists'][key].append(i.params)

            # form case output file base name
            if len(case_list_params) > 1 :
                c['output_files']['output_file_base'] += '_C%03.0f' % (n_case+1)
            if c['settings']['case_output_file_tag'] is not None :
                c['output_files']['output_file_base'] += '_' + c['settings']['case_output_file_tag']

            # now build full params  for a single run
            processor_number += 1  # is one base
            
            c.update({ 'processorID' : n_case,
                    'processor_number': n_case + 1,
                    'code_version_info' : self._U2_code_version_info(),                    
                                })
            full_list_case_params.append(c)
            
        msg_logger.exit_if_prior_errors('Errors in merging case parameters')
        return full_list_case_params


    def _setup_reader_builder(self,run_builder , msg_logger):
        # created a dict which can be used to build a reader

        reader_params =  run_builder['core_classes']['reader']
        if 'class_name' not in  reader_params:
            msg_logger.msg(' Cannot find  "class_name" in reader_class parameter, cannot continue',fatal_error=True, exit_now=True)

        reader = make_class_instance_from_params(reader_params, msg_logger,  class_type_name='reader')
        msg_logger.exit_if_prior_errors() # class name missing or missimg requied variables

        run_builder['core_classes']['reader'] = reader.params # copy full params back to reader
        reader.shared_info.settings = run_builder['settings'] # reader needs access to share parms for set up

        # construct reader_build_info to be used by case_runners to build their reader
        run_builder['reader_build_info'] =  {'reader_params': reader.params}
        reader_build_info = run_builder['reader_build_info']
        reader_build_info['file_info'] , reader_build_info['hindcast_is3D'] = reader.get_hindcast_files_info() # get file lists

        msg_logger.progress_marker('Finished sorting hyrdo-model  files ', tabs=3)

        # read and set up reader grid now as  required for writing grid file
        # also if requested shared grid memory is set up
        # and shared memory info will be also added to reader_build_info for case runner to build reader
        grid = reader.grid
        grid_time_buffers = reader.grid_time_buffers
        nc = reader._open_grid_file(reader_build_info)

        grid = reader.make_non_time_varying_grid(nc, grid)
        grid_time_buffers = reader.make_grid_time_buffers(nc,grid,grid_time_buffers)

        # add information to build grid and grid_time_buffers to reader_build_info
        run_builder = reader.make_grid_builder(grid, grid_time_buffers, run_builder)

        # get information required to build reader fields
        run_builder = reader.maker_field_builder(nc, run_builder)

        # special case read water depth field, so it can be wrtten to grid file
        #add_a_reader_field(self, nc, name, field_variable_comps)

        nc.close()

        if run_builder['settings']['shared_reader_memory']:

            reader_build_info = reader.set_up_shared_grid_memory(reader_build_info)

            # add to class to shared info, alows fileds to be built
            si = reader.shared_info
            si.reset() # needed to clear out old run on windows
            si.classes['reader'] = reader
            si.reader_build_info = reader_build_info

            # set up reader fields
            reader.setup_reader_fields()
            sm_fields = {}
            for name, sm in reader.shared_memory['fields'].items():
                sm_fields[name] = sm.get_shared_mem_map()

        # write grid and save file names
        output_files = run_builder['output_files']
        if run_builder['settings']['write_output_files']:
            output_files['grid'], output_files['grid_outline'] = self._write_run_grid_netCDF(output_files, reader_build_info, reader)

        run_builder['reader_build_info'] = reader_build_info
        return run_builder, reader



    def _write_run_grid_netCDF(self, output_files, reader_build_info, reader):
        # write a netcdf of the grid from first hindcast file
   
        grid= reader.grid

        # get depth from first hincast file
        hindcast = NetCDFhandler(reader_build_info['file_info']['names'][0], 'r')
        depth_var = reader.params['field_variables']['water_depth']
        if depth_var is not None and hindcast.is_var(depth_var):
            # world around to ensure depth read in right format
            field_params,var_info = reader.get_field_variable_info(hindcast,'water_depth',reader.params['field_variables']['water_depth'])
            water_depth = reader.read_file_field_variable_as4D(hindcast,var_info['component_list'][0],var_info['is_time_varying'], file_index=None)
            water_depth = reader.preprocess_field_variable(hindcast,'water_depth',water_depth)
            water_depth= water_depth[0, :, 0, 0] # guard against  water depth being time dependent, so only write first time step of the 4D depth field
            grid['water_depth'] = np.squeeze(water_depth)

        hindcast.close()

        # write grid file
        grid_file = output_files['output_file_base'] + '_grid.nc'
        nc = NetCDFhandler(path.join(output_files['run_output_dir'], grid_file), 'w')
        nc.write_global_attribute('Notes', ' all indices are zero based')
        nc.write_global_attribute('created', str(datetime.now().isoformat()))

        nc.write_a_new_variable('x', grid['x'], ('node_dim', 'vector2D'))
        nc.write_a_new_variable('triangles', grid['triangles'], ('triangle_dim', 'vertex'))
        nc.write_a_new_variable('triangle_area', grid['triangle_area'], ('triangle_dim',))
        nc.write_a_new_variable('adjacency', grid['adjacency'], ('triangle_dim', 'vertex'))
        nc.write_a_new_variable('node_type', grid['node_type'], ('node_dim',), attributesDict={'node_types': ' 0 = interior, 1 = island, 2=domain, 3=open boundary'})
        nc.write_a_new_variable('is_boundary_triangle', grid['is_boundary_triangle'].astype(np.int8), ('triangle_dim',))
        nc.write_a_new_variable('water_depth', grid['water_depth'], ('node_dim',))
        nc.close()

        grid_outline_file = output_files['output_file_base'] + '_grid_outline.json'
        json_util.write_JSON(path.join(output_files['run_output_dir'], grid_outline_file), grid['grid_outline'])

        return grid_file, grid_outline_file

    def _run_case_list(self, case_param_list, reader):
        # do run of all cases
        ml = self.msg_logger
        shared_params = case_param_list[0]['settings']

        if shared_params['shared_reader_memory']:
            case_results = self.run_with_shared_reader(case_param_list, reader)

        elif shared_params['processors'] == 1:
            # serial or non-parallel mode,   run a single case non parallel to makes debuging easier and allows reruns of single case
            case_results = self.run_serial(case_param_list)
        else:
            # run parallel
            case_results = self.run_parallel(case_param_list)

        # unpack case output
        case_info_file_list, case_errors_list = [],[]
        n=0
        for case_file,case_error  in case_results:
            case_info_file_list.append(case_file)
            case_errors_list.append(case_error)
            if case_file is None and case_error is not None:
                ml.msg('CaseInfo files missing  for case ' + str(n)
                                   + ', or other error, case may have not completed, check for .err file!!!!!!, error type= ' + case_error.__class__.__name__, warning=True)
            n += 1

        return case_info_file_list,  case_errors_list

    def run_serial(self,case_param_list):
        # serial or non-parallel mode,   run a single case non parallel to makes debuging easier and allows reruns of single case
        case_results = []
        for c in case_param_list:
            a = self._run1_case(c)
            case_results.append(a)
        return case_results

    def run_parallel(self,case_param_list):
        ml = self.msg_logger
        settings = case_param_list[0]['settings']
        settings['processors'] = min(settings['processors'], len(case_param_list))
        ml.progress_marker('oceantracker:multiProcessing: processors:' + str(settings['processors']))

        with multiprocessing.Pool(processes=settings['processors']) as pool:
            case_results = pool.map(self._run1_case, case_param_list)

        ml.progress_marker('parallel pool complete')

        return case_results

    def run_with_shared_reader(self, case_param_list, reader):
        # complete reader build then run cases asynchronously with shared reader

        # shared grid is already to allow grid to be written  and case

        # shared reader control arrays
        time_steps_in_buffer = shared_memory.SharedMemArray(shape=(2,0),dtype=np.int2)
        c = {'time_steps_in_buffer':time_steps_in_buffer.get_shared_mem_map() }

        # add field shared memory build info to each case
        for case in case_param_list:
            case['reader_build_info']['shared_memory']['fields'] = sm_fields
            case['reader_build_info']['shared_memory']['control'] = c
        #todo make shared fields and grid arrays read only in children sm.data.setflags(write=False)
        pass

    @staticmethod
    def _run1_case(run_params):
        # run one process on a particle based on given family class parameters
        # by creating an independent instances of  model classes, with given parameters
        ot = OceanTrackerCaseRunner()
        caseInfo_file, case_error = ot.run(deepcopy(run_params))
        return caseInfo_file, case_error

    def _U1_get_all_cases_performance_info(self, case_info_files, case_error_list, sparams, run_output_dir, t0):
        # finally get run totals of steps and particles

        n_time_steps = 0.
        total_alive_particles = 0
        # load log files to get info on run from solver info
        for n, case_file, case_error in zip(range(len(case_info_files)), case_info_files,case_error_list) :
            if case_file is not None :
                c= json_util.read_JSON(path.join(run_output_dir, case_file))
                sinfo = c['class_info']['solver']
                n_time_steps += sinfo['time_steps_completed']
                total_alive_particles += sinfo['total_alive_particles']

        num_cases = len(case_info_files)

        # JSON parallel run info data
        d = {'processors': sparams['processors'],
             'num_cases': num_cases,
                'elapsed_time' :perf_counter() - t0,
            'average_active_particles': total_alive_particles / num_cases if num_cases > 0 else None,
             'average_number_of_time_steps': n_time_steps/num_cases  if num_cases > 0 else None,
             'particles_processed_per_second': total_alive_particles /(perf_counter() - t0)
             }

        # put parallel info in first file base
        return d

    def _U2_code_version_info(self):
        try:
            git_revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=path.dirname(path.realpath(__file__))).decode().replace('\n','')
        except:
            git_revision = 'unknown'
        return { 'version': code_version, 'git_revision': git_revision, 'python_version':version}

    def check_top_level_class_keys(self, params, template, msg_logger, required_keys=[], crumbs=None):
        # ensure top level parameter dict has all keys, and any required ones
        required_type = type(list(template.values())[0])
        tag = '_class' if required_type == dict else '_list'

        # check for required keys
        for key in required_keys:
            if key not in params:
                msg_logger.msg('Required param key  "' + key + tag + '" is missing', crumbs=crumbs, fatal_error=True)

        # make sure all template keys are present
        for key, item in template.items():
            if key not in params:
                params[key] = item
            elif  type(params[key]) != type(item):
                msg_logger.msg('Param key = "' + key + '" must be type ' + str(type(item)) + ' not type' + str(type(params[key])),crumbs=crumbs, fatal_error=True)

        # key for unexpected keys
        for key in params.keys():
            if key not in template.keys():
                msg_logger.msg( 'Unexpected key = "' + key + tag + '" in parameters', warning=True, crumbs=crumbs,
                      hint= 'must be one of keys = ' + str(list(template.keys())), tabs=1)

        return params

def convert_ver03_to_04_params(params):
    p= {}
    p.update(params['shared_params'])
    p['reader_class'] = params['reader']

    if 'base_case_params' in params:
        p.update(convert03_case_params_to04(params['base_case_params']))
        p.update(params['shared_params'])
    if 'case_list_params' in params:
        p['case_list_params'] = []
        for c in params['case_list_params']:
            p['case_list_params'].append(convert03_case_params_to04(c))
    return p

def convert03_case_params_to04(case_params):
    p ={}
    for key in case_params.keys():
        if key=='settings':
            p.update(case_params[key])
        else:
            name = copy(key)
            if key in common_info.class_lists: name += '_list'
            if key in common_info.core_classes: name += '_class'
            p[name] = case_params[key]
    return p
