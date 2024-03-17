# method to run ocean tracker from parameters
# eg run(params)
import sys

#


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
import shutil
from time import perf_counter
from copy import  copy
import numpy as np
import difflib

from oceantracker.util import setup_util
from oceantracker import common_info_default_param_dict_templates as common_info
from oceantracker.util.class_importing_util import ClassImporter

from oceantracker.util.parameter_checking import merge_params_with_defaults

from oceantracker.util import json_util ,yaml_util, get_versions_computer_info
from oceantracker.util.messgage_logger import GracefulError, MessageLogger
from oceantracker.reader.util import get_hydro_model_info

import traceback
import os
OTname = common_info.package_fancy_name
help_url_base = 'https://oceantracker.github.io/oceantracker/_build/html/info/'

def run(params):
    ot = _OceanTrackerRunner()
    case_info_files = ot.run(params)
    return case_info_files

def run_parallel(base_case_params, case_list_params=[{}]):
    ot= _OceanTrackerRunner()
    case_info_files  = ot.run(base_case_params, case_list_params)
    return case_info_files

class OceanTracker():
    def __init__(self,params=None):
        self.params= param_template() if params is None else params
        self.case_list_params=[]

        self.msg_logger = MessageLogger('helper')
        self.msg_logger.print_line()
        self.msg_logger.msg('Starting OceanTracker helper class')

    # helper methods
    def settings(self,case=None, **kwargs):

        # work out if to add to base params or case list params
        p = self._check_case(case, f'adding settings "{str(kwargs.keys())}"')
        for key in kwargs:
            p[key]= kwargs[key]

    def add_class(self, class_role =None,case=None, **kwargs):
        ml = self.msg_logger
        known_class_roles = common_info.class_dicts_list + common_info.core_class_list
        if class_role is None:
            ml.msg('oceantracker.add_class, must give first parameter as class role, eg. "release_group"', fatal_error=True, caller =self)
            return

        if type(class_role) != str:
            ml.msg(f'oceantracker.add_class, class_role must be a string', fatal_error=True, caller=self,
                   hint='Given type =' + str(type(class_role)))
            return

        if class_role not in known_class_roles:
            ml.msg(f'oceantracker.add_class, class_role parameter is not recognised, value ="{class_role}"', fatal_error=True,
                   hint= 'Must be one of-' +str(known_class_roles))
            return

        if class_role in common_info.core_class_list:
            # single class
            self.params[class_role] = kwargs

        elif class_role in common_info.class_dicts_list:
            # can have more than one class
            if 'name' not in kwargs:
                ml.msg(f'add_class_dict() : for class_role"{class_role}" must have a name parameter, eg. name=my_{class_role}1, ignoring this class' ,
                       fatal_error=True,
                       hint= f'the can be more than one {"{class_role}"}, so each must have a unique name')
                return
            name =kwargs["name"]
            if name in self.params[class_role]:
                ml.msg(f'class type "{class_role}" already has a class named "{name}", ignoring later versions', warning=True)
            else:
                # add params to  class type of given name and case
                p = self._check_case( case, f'adding class "{name}"') # find if added to base case or case list # case

                p[class_role][name]={} # add blank param dict
                for key in kwargs:
                    if key != 'name':
                        p[class_role][name][key] = kwargs[key]
        else:
            ml.spell_check('ignoring class helper function add_class(),',class_role, known_class_roles,
                        crumbs=f' in add_class() class type "{class_role}"', caller = self)
        pass

    def _check_case(self,case, crumbs):
        # check and make space for case
        ml = self.msg_logger

        # no case number given add to ordinary params
        if case is None: return self.params

        if type(case) != int or case < 0:
            ml.msg(f'Parallel case number nust be a integer >0, got case # "{str(case)}"', fatal_error=True, crumbs=crumbs)
            return

        # expand case list to fit case if needed
        if len(self.case_list_params) < case+1 :
            ml.progress_marker(f'Adding parallel case number # "{str(case)}"')
            for n in range(len(self.case_list_params), case + 1):
                self.case_list_params.append(deepcopy(param_template()))

        return self.case_list_params[case]

    def run(self):
        self.msg_logger.progress_marker('Starting run using helper class')
        ot= _OceanTrackerRunner()
        # todo print helper message here at end??
        ot.helper_msg_logger = self.msg_logger  # used to print helper messages at end and write to file

        case_info_file = ot.run(self, self.params, self.case_list_params)

        self.msg_logger.close()

        return case_info_file

class _OceanTrackerRunner(object):
    def __init__(self):
        self.msg_logger = MessageLogger('main')

    def run(self,params, case_list_params = None):
        ml = self.msg_logger
        run_builder= dict()
        self.start_t0 = perf_counter()
        self.start_date = datetime.now()

        ml.print_line()
        ml.msg(f'{OTname} starting main:')
        ml.progress_marker('Output dir set up.')
        # setup output dir and msg files
        o = setup_util.setup_output_dir(params, self.msg_logger, crumbs='setting up output dir')
        o['run_log'], o['run_error_file'] = ml.set_up_files(o['run_output_dir'], o['output_file_base'])
        run_builder['output_files'] = o

        setup_util.write_raw_user_params(run_builder['output_files'],params,ml, case_list=case_list_params)

        # set numba config environment variables, before any import of numba, eg by readers,
        setup_util.config_numba_environment(params, ml, crumbs='main setup', caller=self)  # must be done before any numba imports

        params = setup_util.merge_settings_with_defaults(params, ['debug'],ml,crumbs='Debug setting ') # make sure needed values are set or use defaults

        ml.msg(f'Output is in dir "{o["run_output_dir"]}"', hint='see for copies of screen output and user supplied parameters, plus all other output',tabs=1)
        ml.print_line()
        ml.msg(f' {OTname} version {common_info.code_version} - preliminary setup')

        self._prelimary_checks(params)
        ml.exit_if_prior_errors('parameters have errors')

        # get list of files etc in time order
        run_builder = self._get_hindcast_file_info(params, run_builder)

        if case_list_params is  None or len(case_list_params) == 0:
            # no case list
            case_info_file = self._run_single(params, run_builder)
        else:
            # run // case list with params as base case defaults for each run
            case_info_file = self._run_parallel(params, case_list_params, run_builder)

        ml.close()
        return case_info_file

    def _run_single(self, user_given_params, run_builder):

        ml = self.msg_logger


        run_builder['caseID'] = 0 # tag case
        run_builder['working_params'] = setup_util.decompose_params(user_given_params, ml, crumbs='_run_single>', caller=self)

        # check defaults of all settings
        run_builder['working_params']['settings'] = merge_params_with_defaults(run_builder['working_params']['settings'],
                                                            common_info.all_setting_defaults, ml, crumbs='_runs_single checking defaults>',
                                                            caller=self, check_for_unknown_keys=True)


        # try catch is needed in notebooks to ensure mesage loger file is close,
        # which allows rerunning in notebook without  permission file errors
        try:
            # keep oceantracker_case_runner out of main namespace
            from oceantracker.oceantracker_case_runner import OceanTrackerCaseRunner
            # make instance of case runer and run it with decomposed working params
            ot_case_runner = OceanTrackerCaseRunner()
            case_info_file, case_msg = ot_case_runner.run_case(run_builder)

        except Exception as e:
            # ensure message loggers are closed

            print(str(e))
            traceback.print_exc()
            return None

        # check is case ran
        if case_info_file is None:
            ml.msg('case_info_file is None, run may not have completed', fatal_error=True)

        self._main_run_end(case_info_file, len(case_msg['errors']), len(case_msg['warnings']), len(case_msg['notes']))

        return case_info_file

    def _prelimary_checks(self,params):
        ml = self.msg_logger

        setup_util.check_python_version(ml)
        # check for compulsory classes
        # check reader params
        if 'reader' not in params or len(params['reader']) < 2:
            ml.msg('Parameter "reader" is required, or missing required parameters',
                           hint='Add a "reader" top level key to parameters with a dictionary containing  at least "input_dir" and "file_mask" keys and values',
                           fatal_error=True, crumbs='case_run_set_up', caller=self)


    def _main_run_end(self,case_info_files, num_case_errors,num_case_warnings,num_case_notes):
        # final info output
        ml = self.msg_logger
        self._write_run_info_json(case_info_files, self.start_t0)

        ml.show_all_warnings_and_errors()

        # rewite any help class error/warnings
        if hasattr(self,'helper_msg_logger'):
            ml_helper = self.helper_msg_logger
            for l in ml_helper.errors_list:
                ml.msg(l)

            for l in ml_helper.warnings_list:
                ml.msg(l)

        ml.print_line()
        ml.msg(f'OceanTracker summary:  elapsed time =' + str(datetime.now() - self.start_date),)

        ml.msg(f'Cases - {num_case_errors:3d} errors, {num_case_warnings:3d} warnings, {num_case_notes:3d} notes, check above', tabs=3)
        if hasattr(self, 'helper_msg_logger'):
            ml.msg(f'Helper- {len(ml_helper.errors_list):3d} errors, {len(ml_helper.warnings_list):3d} warnings, {len(ml_helper.notes_list):3d} notes, check above', tabs=3)
        ml.msg(f'Main  - {len(ml.errors_list):3d} errors, {len(ml.warnings_list):3d} warnings, {len(ml.notes_list):3d} notes, check above', tabs=3)

        ml.print_line()
        ml.close()

    def _run_parallel(self,base_case_params, case_list_params, run_builder):
        # run list of case params
        ml = self.msg_logger
        self.helper_msg_logger =  ml # keep references to write message at end as runs has main message logger

        # split base_case_params into shared and case params
        base_working_params = setup_util.decompose_params(base_case_params, ml, caller=self, crumbs='_run_parallel decompose base case params')

        base_working_params['settings'] = setup_util.merge_settings_with_defaults(base_working_params['settings'], 'processors', ml,
                                                                                                crumbs='setting processor number', caller=self)

        case_run_builder_list=[]

        for n_case, case_params in enumerate(case_list_params):
            case_working_params = setup_util.decompose_params(case_params, msg_logger=ml, caller=self)

            case_working_params =  setup_util.merge_base_and_case_working_params(base_working_params, case_working_params, ml,
                                                                                 crumbs=f'_run_parallel case #[{n_case}]', caller=None)

            # get any missing settings from defaults after merging with base case settings
            case_working_params['settings'] = merge_params_with_defaults(case_working_params['settings'],
                                                                        common_info.all_setting_defaults, ml, crumbs=f'_run_parallel case #[{n_case}]',
                                                                        caller=self, check_for_unknown_keys=True)


            ml.exit_if_prior_errors(f'Errors in setting up case #{n_case}')
            case_run_builder = deepcopy(run_builder)
            case_run_builder['caseID'] = n_case
            case_run_builder['working_params'] = case_working_params

            # add and tweak output file info
            case_run_builder['output_files']['output_file_base'] += '_C%03.0f' % (n_case)

            # now add builder to list to run
            case_run_builder_list.append(case_run_builder)



        # do runs
        num_proc = base_working_params['settings']['processors']
        ml.progress_marker(' oceantracker:multiProcessing: processors:' + str(num_proc))

        # run // cases
        with multiprocessing.Pool(processes=num_proc) as pool:
            case_results = pool.map(self._run1_case, case_run_builder_list)

        ml.progress_marker('parallel pool complete')

        # get case files and  error/warning counts
        case_msg=[]
        num_warnings= 0
        num_errors = 0
        num_notes = 0
        case_info_files =[]
        for n, c in  enumerate(case_results):
            if c[0] is None:
                ml.msg(f'Case #{n:03}, failed to finish', fatal_error=True, hint='hint see above, any case errors detected are listed below')
                for m in c[1]['errors']:
                    ml.msg(m, tabs =2)
            case_info_files.append(c[0])
            case_msg.append(c[1])
            num_errors += len(c[1]['errors'])
            num_warnings += len(c[1]['warnings'])
            num_notes += len(c[1]['notes'])

        self._main_run_end(case_info_files,num_errors,num_warnings,num_notes)
        return case_info_files

    @staticmethod
    def _run1_case(working_params):
        # run one process on a particle based on given family class parameters
        # by creating an independent instances of  model classes, with given parameters

        # keep oceantracker_case_runner out of main namespace
        from oceantracker.oceantracker_case_runner import OceanTrackerCaseRunner

        ot = OceanTrackerCaseRunner()
        caseInfo_file, return_msgs= ot.run_case(deepcopy(working_params))
        return caseInfo_file, return_msgs


    def _get_hindcast_file_info(self, params, run_builder ):
        # created a dict which can be used to build a reader
        t0= perf_counter()
        ml = self.msg_logger
        class_importer= ClassImporter(path.dirname(__file__), msg_logger=ml)
        reader_params =  params['reader']

        if 'input_dir' not in reader_params or 'file_mask' not in reader_params:
            ml.msg('Reader class requires settings, "input_dir" and "file_mask" to read the hindcast',fatal_error=True, exit_now=True )
        # check input dir exists

        if path.isdir(reader_params['input_dir']):
            ml.progress_marker(f'Found input dir "{reader_params["input_dir"]}"')
        else:
            ml.msg(f' Could not find input dir "{reader_params["input_dir"]}"',
                   hint ='Check reader parameter "input_dir"', fatal_error=True, exit_now=True)

        reader_params, file_list = get_hydro_model_info.find_file_format_and_file_list(reader_params,class_importer, ml)


        reader = class_importer.new_make_class_instance_from_params(reader_params,'reader',  default_classID='reader', crumbs='primary reader>')

        ml.exit_if_prior_errors() # class name missing or missing required variables
        run_builder['reader_builder']=dict(params=reader_params,
                                              file_info= reader.get_hindcast_files_info(file_list, ml)
                                              )

        ml.progress_marker('sorted hyrdo-model files in time order', start_time=t0)

        # get file info for nested readers
        run_builder['nested_reader_builders']= {}
        if 'nested_readers' not in params: params['nested_readers'] ={}

        for name, params in params['nested_readers'].items():
            t0 = perf_counter()
            nested_params, nested_file_list = get_hydro_model_info.find_file_format_and_file_list(params,class_importer, ml)
            nested_reader = class_importer.new_make_class_instance_from_params( nested_params,'reader', default_classID='reader', crumbs=f'nested reader{name}>')

            d= dict(params=nested_params,
                    file_info= nested_reader.get_hindcast_files_info(nested_file_list, ml)
                    )
            run_builder['nested_reader_builders'][name]=d
            ml.progress_marker(f'sorted nested hyrdo-model files in time order{name}', start_time=t0)

        return run_builder





    def _write_run_info_json(self, case_info_files, t0):
        # read first case info for shared info
        ml = self.msg_logger
        ci = deepcopy(case_info_files) # dont alter input
        if type(ci) is not list: ci= [ci]

        # finally get run totals of steps and particles across al cases and write
        n_time_steps = 0.
        total_alive_particles = 0
        case_info_list=[]
        # load log files to get info on run from solver info
        for n, case_file  in enumerate(ci) :

            if case_file is not None :
                c= json_util.read_JSON(case_file)
                n_time_steps +=  c['run_info']['time_steps_completed']
                total_alive_particles += c['run_info']['total_alive_particles']
                case_info_list.append(path.basename(case_file))
            else:
                case_info_list.append((None))
                ml.msg(f'Case #{n:d} has no case info file, likely has crashed',warning=True)

        num_cases = len(ci)

        # JSON parallel run info data
        d = {'output_files' :{},
            'version_info': get_versions_computer_info.get_code_version(),
            'computer_info': get_versions_computer_info.get_computer_info(),
            'num_cases': num_cases,
            'elapsed_time' :perf_counter() - t0,
            'average_active_particles': total_alive_particles / num_cases if num_cases > 0 else None,
            'average_number_of_time_steps': n_time_steps/num_cases  if num_cases > 0 else None,
            'particles_processed_per_second': total_alive_particles /(perf_counter() - t0)
             }

        # get output file names
        c0= json_util.read_JSON(ci[0])
        o = c0['output_files']
        d['output_files'] = {'root_output_dir': o['root_output_dir'],
                            'run_output_dir': o['run_output_dir'],
                            'output_file_base': o['output_file_base'],
                            'runInfo_file': o['runInfo_file'],
                            'runLog_file': o['runLog_file'],
                            'run_error_file': o['run_error_file'],
                            'users_params_json': o['raw_user_params'],
                             'caseInfo_files':case_info_list
                             }
        json_util.write_JSON(path.join(o['run_output_dir'],o['runInfo_file']),  d)
        ml.msg('run summary with case file names   "' + o['runInfo_file'] + '"',  tabs=2, note=True)

    def close(self):
        pass

def param_template():
    # return an empty parameter dictionary

    d = {}
    all_settings=list(common_info.shared_settings_defaults.keys())+ list(common_info.case_settings_defaults.keys())
    for key in all_settings:
        d[key] = None

    for key in sorted(common_info.class_dicts_list):
        d[key] = {}
    for key in sorted(common_info.core_class_list):
        d[key] = {}
    return deepcopy(d)