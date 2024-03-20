import os
from copy import deepcopy
from os import path, environ, remove
from oceantracker.util.parameter_base_class import ParameterBaseClass

from time import  perf_counter
from oceantracker.util.messgage_logger import MessageLogger, GracefulError
from oceantracker.util import profiling_util, get_versions_computer_info
import numpy as np
from oceantracker.util import time_util, numba_util
from oceantracker.util import json_util
from datetime import datetime
from time import sleep
import traceback
from oceantracker.util.class_importing_util import ClassImporter
from oceantracker.util.setup_util import config_numba_environment
from oceantracker import common_info_default_param_dict_templates as common_info

from oceantracker.shared_info import SharedInfo as si


# note do not import numba here as its environment  setting must ve done first, import done below

from oceantracker.util.package_util import get_all_parameter_classes

class OceanTrackerCaseRunner(ParameterBaseClass):
    # this class runs a single case
    def __init__(self):
        # set up info/attributes
        super().__init__()  # required

    def run_case(self, run_builder):
        si._setup() # clear out classes from class instance of SharedInfo if running series of mains
        si.msg_logger.set_screen_tag(f'C{run_builder["caseID"]:03d}')
        si.msg_logger.settings(max_warnings=run_builder['working_params']['settings']['max_warnings'])

        d0 = datetime.now()
        t_start = perf_counter()
        # used to write a unaltered version to case_info.json
        self.raw_working_params = deepcopy(run_builder['working_params'])

        # set up message logging
        output_files = run_builder['output_files']

        # set numba config environment variables, before any import of numba, eg by readers,
        # also done in main but also needed here for parallel runs
        config_numba_environment(run_builder['working_params']['settings'],  si.msg_logger, crumbs='Starting case_runner', caller=self)  # must be done before any numba imports

        # basic  shortcuts
        si.run_builder = run_builder
        si.caseID = run_builder['caseID']
        si.working_params = run_builder['working_params']
        si.settings = si.working_params['settings']

        si.output_files = run_builder['output_files']
        si.run_output_dir = si.output_files['run_output_dir']
        si.output_file_base = si.output_files['output_file_base']
        output_files['case_log_file'], output_files['case_error_file'] = \
                    si.msg_logger.set_up_files(output_files['run_output_dir'], output_files['output_file_base'] + '_caseLog')

        si.msg_logger.print_line()
        si.msg_logger.msg('Starting case number %3.0f, ' % si.caseID + ' '
                                      + si.output_files['output_file_base']
                                      + ' at ' + time_util.iso8601_str(datetime.now()))
        si.msg_logger.print_line()

        # setup class importer
        t0=perf_counter()
        si.class_importer = ClassImporter(path.dirname(__file__), msg_logger=si.msg_logger)
        si.msg_logger.progress_marker('Scanned OceanTracker to build short name map to the full class_names', start_time=t0)

        # other useful shared values
        si.backtracking = si.settings['backtracking']
        si.model_direction = -1 if si.backtracking else 1
        si.time_of_nominal_first_occurrence = -si.model_direction * 1.0E36

        si.write_output_files = si.settings['write_output_files']
        si.z0 = si.settings['z0']
        si.minimum_total_water_depth = si.settings['minimum_total_water_depth']
        si.computer_info = get_versions_computer_info.get_computer_info()

        # set up profiling
        profiling_util.set_profile_mode(si.settings['profiler'])

        case_info_file = None
        return_msgs = {'errors': si.msg_logger.errors_list, 'warnings': si.msg_logger.warnings_list, 'notes': si.msg_logger.notes_list}

        # run info
        si.run_info = dict(time_step = si.settings['time_step'],
                        backtracking=si.settings['backtracking'],
                        model_direction= -1 if si.settings['backtracking'] else 1)


        ri = si.run_info
        #useful short cuts



        # case set up
        try:
            self._set_up_run()
            self._do_pre_processing()
            self._make_core_classes()

            # add any intergrated model, to get params changes, but don't initialise
            # there can only be one model added
            params = si.working_params['core_classes']['integrated_model']
            if len(params) > 0:
                im = si.add_core_role('integrated_model', params, crumbs='adding  "integrated_model" class',initialise=False)
                # any changes to params by model,
                # e.g. release groups from model which affect start and end times
                im.have_called_add_settings_and_class_params = False
                im.add_settings_and_class_params()
                im.have_called_add_settings_and_class_params = True # this is used to block use of im.add_class and im.settings

            #todo remerge/check settings changed by model here?

            # now  params commlet cabn set up

            # set up start time and duration based on particle releases
            self._setup_particle_release_groups_and_start_end_times(si.working_params['class_dicts']['release_groups'])
            self._make_and_initialize_user_classes()
            self._finalize_classes()  # any setup actions that mus be done after other actions, eg shedulers

            # model may depend on ther classes, so intilialise after all other claases are setup
            if si.core_roles.integrated_model is not None:
                im.initial_setup()
                im.final_setup()

        except GracefulError as e:
            si.msg_logger.show_all_warnings_and_errors()
            si.msg_logger.write_error_log_file(e)
            si.msg_logger.msg(f' Case Funner graceful exit from case number [{si.caseID:2}]', hint ='Parameters/setup has errors, see above', fatal_error= True)

            if si.settings['debug']:
                si.msg_logger.write_error_log_file(e)
                si.msg_logger.write_error_log_file(traceback.print_exc())

            si.msg_logger.close()
            return None, return_msgs

        except Exception as e:
            si.msg_logger.show_all_warnings_and_errors()
            si.msg_logger.write_error_log_file(e)
            si.msg_logger.write_error_log_file(traceback.print_exc())
            si.msg_logger.msg(f' Unexpected error in case number [{si.caseID:2}] ', fatal_error=True,hint='check above or .err file')
            si.msg_logger.close()
            return  None, return_msgs


        # check particle properties have other particle properties, fields and other compatibles they require
        self._do_run_integrity_checks()


        try:
            self._do_a_run()
            case_info = self._get_case_info(d0,t_start)

            if si.settings['write_output_files']:
                # write grid if first case
                if si.caseID == 0:
                    si.core_roles.field_group_manager.write_hydro_model_grid()

                case_info_file = si.output_file_base + '_caseInfo.json'
                json_util.write_JSON(path.join(si.run_output_dir, case_info_file), case_info)

        except GracefulError as e:
            si.msg_logger.show_all_warnings_and_errors()
            si.msg_logger.write_error_log_file(e)

            si.msg_logger.msg(f' Case Funner graceful exit from case number [{si.caseID:2}]', hint ='Parameters/setup has errors, see above', fatal_error= True)

            if si.settings['debug']:
                si.msg_logger.write_error_log_file(e)
                si.msg_logger.write_error_log_file(traceback.print_exc())
            si.msg_logger.close()
            return None, return_msgs

        except Exception as e:
            si.msg_logger.show_all_warnings_and_errors()
            si.msg_logger.write_error_log_file(e)
            si.msg_logger.write_error_log_file(traceback.print_exc())
            si.msg_logger.msg(f' Unexpected error in case number [{si.caseID:2}] ', fatal_error=True,hint='check above or .err file')
            si.msg_logger.close()
            return  None, return_msgs

        # reshow warnings
        si.msg_logger.show_all_warnings_and_errors()

        si.msg_logger.print_line()
        si.msg_logger.progress_marker('Finished case number %3.0f, ' % si.caseID + ' '
                                      + si.output_files['output_file_base']
                                      + ' started: ' + str(d0)
                                      + ', ended: ' + str(datetime.now()))
        si.msg_logger.msg('Elapsed time =' + str(datetime.now() - d0), tabs=3)
        si.msg_logger.print_line()

        self.close()  # close al classes and msg logger
        si.msg_logger.close()
        case_info_file = path.join(si.run_output_dir,case_info_file)
        return case_info_file,  return_msgs

    def _set_up_run(self):
        # builds shared_info class variable with data and classes initialized  ready for run
        # from single run case_runner_params

        # get short class names map
        # delay  start, which may avoid occasional lockup at start if many cases try to read same hindcast file at same time
        if si.settings['multiprocessing_case_start_delay'] is not None:
            delay = si.settings['multiprocessing_case_start_delay'] * (si.caseID % si.settings['processors'])
            si.msg_logger.progress_marker('Delaying start by  ' + str(delay) + ' sec')
            sleep(delay)
            si.msg_logger.progress_marker('Starting after delay  of ' + str(delay) + ' sec')


        if si.settings['use_random_seed']:
            np.random.seed(0) # set numpy
            from oceantracker.util.numba_util import seed_numba_random
            seed_numba_random(0)
            si.msg_logger.msg('Using numpy.random.seed(0), makes results reproducible (only use for testing developments give the same results!)',warning=True)

    #@function_profiler(__name__)
    def _do_a_run(self):
        # build and run solver from parameter dictionary
        # run from a given dictionary to enable particle tracking on demand from JSON type parameter set
        # also used for parallel  version
        info= self.info
        info['model_run_started'] = datetime.now()

        solver = si.core_roles.solver

        # fill and process buffer until there is less than 2 steps
        si.msg_logger.print_line()
        si.msg_logger.progress_marker('Starting ' + si.output_file_base + ',  duration: ' + time_util.seconds_to_pretty_duration_string(si.run_info['duration']))

        t0 = perf_counter()
        si.msg_logger.progress_marker('Initialized Solver Class', start_time=t0)

        # ------------------------------------------
        solver.solve()
        # ------------------------------------------
        pass

    def _do_run_integrity_checks(self):
         # check all have required, fields, part props and grid data
        for i in si.all_class_instance_pointers_iterator():
            i.check_requirements()

        si.msg_logger.exit_if_prior_errors()



    def _do_pre_processing(self):
        # do pre-processing, eg read polygons from files
       pass

    def _setup_particle_release_groups_and_start_end_times(self, particle_release_groups_params_dict):
        # particle_release groups setup and instances,
        # find extremes of  particle existence to calculate model start time and duration
        t0 = perf_counter()
        pgm = si.core_roles.particle_group_manager
        ri = si.run_info
        # set up to start end times based on release_groups

        # set up release groups and find first release time to start model
        first_time = []
        last_time = []
        md = si.model_direction

        for name, rg_params in particle_release_groups_params_dict.items():
            # make instance and initialise
            rg = pgm.add_release_group(name, rg_params)
            start, life_span = rg.get_start_time_and_life_span()

            first_time.append(start)
            last_time.append( start + md * life_span)

        # set model run start/end time allowing for back tracking
        ri['start_time'] = np.min(md * np.asarray(first_time)) * md
        ri['end_time']   = np.max(md * np.asarray(last_time)) * md
        # clip end time to be within hincast
        if si.backtracking:
            ri['end_time'] = max(ri['end_time'],si.hindcast_info['start_time'])
        else:
            ri['end_time'] = min(ri['end_time'], si.hindcast_info['end_time'])

        # make release times array
        ri['times']    = np.arange(ri['start_time'], ri['end_time'], ri['time_step'] * md)
        ri['end_time'] = ri['times'][-1] # adjust end to nearest time step

        # setup scheduler for each release group, how start and model times known
        #   this will round start times and release interval to be integer number of model time steps after the start
        for name, rg in si.roles.release_groups.items():
              si.add_sheduler_to_class('release_scheduler', rg,
                                        start= rg.params['start'],  # set above from start date
                                        interval= rg.params['release_interval'],
                                        end     = rg.params['end'],
                                        duration= rg.params['duration'], caller=rg)

        if len(si.roles.release_groups) == 0:
            # guard against there being no release groups
            si.msg_logger.msg('No valid release groups' , fatal_error=True, exit_now=True, caller=self)

        # useful information
        ri['duration'] = abs(ri['times'][-1] - ri['times'][0])
        ri['start_date'] = time_util.seconds_to_isostr(ri['times'][0])
        ri['end_date'] = time_util.seconds_to_isostr(ri['times'][-1])

        si.msg_logger.progress_marker('Set up run start and end times, plus release groups and their schedulers', start_time=t0)
    def _make_core_classes(self):
        # initialise all classes, order is important!
        # shortcuts
        t0 = perf_counter()

        # start with setting up field gropus, which set up readers
        # as it has info on whether 2D or 3D which  changes class options'
        # reader prams should be full and complete from oceanTrackerRunner, so dont initialize
        # chose fiel manager for normal or nested readers
        if len(si.working_params['class_dicts']['nested_readers']) > 0:
            # use devopment nested readers class
            si.working_params['core_classes']['field_group_manager'].update(dict(class_name='oceantracker.field_group_manager.dev_nested_fields.DevNestedFields'))

        # set up feilds
        fgm = si.add_core_role('field_group_manager', si.working_params['core_classes']['field_group_manager'], crumbs=f'adding core class "field_group_manager" ')
        fgm.initial_setup()  # needed here to add reader fields inside reader build
        si.is3D_run = fgm.info['is3D']
        si.run_info['is3D_run']  = fgm.info['is3D']

        fgm.setup_dispersion_and_resuspension()  # setup files required for didpersion etc, depends on what variables are in the hydro-files
        fgm.final_setup()

        si.hindcast_info = fgm.get_hydro_model_info()


        if si.settings['write_tracks']:
            si.add_core_role('tracks_writer',si.working_params['core_classes']['tracks_writer'], initialise=True)
        else:
            si.core_roles.tracks_writer = None

        pgm = si.add_core_role('particle_group_manager', si.working_params['core_classes']['particle_group_manager'], crumbs=f'adding core class "particle_group_manager" ')
        pgm.initial_setup()  # needed here to add reader fields inside reader build


        # set up particle properties associated with fields etc
        fgm.add_part_prop_from_fields_plus_book_keeping()
        core_role_params = si.working_params['core_classes']

        # make other core classes
        si.add_core_role('solver', core_role_params['solver'], crumbs='core class solver ')


        if  si.settings['time_step'] >  si.hindcast_info['time_step']:
            si.msg_logger.msg(f'Results may not be accurate as, time step param={si.settings["time_step"]:2.0f} sec,  > hydo model time step = {si.hindcast_info["time_step"]:2.0f}',
                              warning=True)

    def _finalize_classes(self):
        # finalise the classeds
        #todo , more needed here to finalise othe classes?
        t0 = perf_counter()
        ri = si.run_info
        # do final setp which may depend on settings from intitial st up

        # order matters, must do interpolator after particle_group_manager, to get stucted arrays and solver last
        si.core_roles.particle_group_manager.final_setup()
        si.core_roles.solver.final_setup()

        si.msg_logger.progress_marker('final set up of core classes', start_time=t0)

    def _make_and_initialize_user_classes(self):
        # complete build of particle by adding reade, custom properties and modifiers
        pgm = si.core_roles.particle_group_manager

        # any custom particle properties added by user
        for name, p in si.working_params['class_dicts']['particle_properties'].items():
            pgm.add_particle_property(name, 'user',p)


        # build and initialise other user classes, which may depend on custom particle props above or reader field, not sure if order matters
        for user_type in ['velocity_modifiers','trajectory_modifiers',
                             'particle_statistics', 'particle_concentrations', 'event_loggers']:
            for name, params in si.working_params['class_dicts'][user_type].items():
                i = si.add_user_class(user_type,name, params,class_type= 'user',   crumbs=' making class type ' + user_type + ' ')
                i.initial_setup()  # some require instanceID from above add class to initialise


        pass

    # ____________________________
    # internal methods below
    # ____________________________
    def _get_case_info(self, d0, t0):
        pgm= si.core_roles.particle_group_manager
        info = self.info
        info['date_of_time_zero'] = time_util.seconds_to_datetime64(np.asarray([0.]))

        info['backtracking'] = si.backtracking
        elapsed_time_sec = perf_counter() -t0
        info.update(dict(started=str(d0), ended=str(datetime.now()),
                         duration=str(datetime.now() - d0),
                         elapsed_time_sec=elapsed_time_sec,
                         number_particles_released= pgm.info['particles_released']))
        info.update(si.run_info)
        # base class variable warnings is common with all descendents of parameter_base_class
        d = {'user_note': si.settings['user_note'],
             'output_files': si.output_files,
             'caseID' : si.caseID,
             'block_timings': [],
             'version_info': get_versions_computer_info.get_code_version(),
            'computer_info': get_versions_computer_info.get_computer_info(),
            'file_written': datetime.now().isoformat(),
             'run_info' : info,
             'hindcast_info': si.core_roles.field_group_manager.info,
             'working_params': self.raw_working_params,
             'full_case_params': si.working_params,
             'particle_status_flags': si.particle_status_flags.as_dict(),
             'release_groups' : {},
             'particle_release_group_user_maps': si.core_roles.particle_group_manager.get_release_group_userIDmaps(),
             'errors': si.msg_logger.errors_list,
             'warnings': si.msg_logger.warnings_list,
             'notes': si.msg_logger.notes_list,
             'function_timers': {},
             'class_roles_info': {}, }

        for key, i in si.roles.as_dict().items():
            if i is None : continue

            d['class_roles_info'][key] = {}
            d['output_files'][key] = {}
            # interate over dict
            for key2, i2 in i.items():
                d['class_roles_info'][key][key2]= i2.info
                d['output_files'][key][key2]= i2.info['output_file'] if 'output_file' in i2.info else None
                if hasattr(i2,'scheduler_info'):
                    d['class_roles_info'][key][key2]['scheduler_info'] = i2.scheduler_info

        # core roles
        for key, i in si.core_roles.as_dict().items():
            if i is None: continue
            d['class_roles_info'][key] = i.info
            d['output_files'][key] = i.info['output_file'] if 'output_file' in i.info else None
            if hasattr(i, 'scheduler_info'):
                d['class_roles_info'][key]['scheduler_info'] = i.scheduler_info

        # add basic release group info
        #do do
        for key, item in si.roles.release_groups.items():
            # grid does not have points pararm

            rginfo={'points': item.points if hasattr(item,'points') else item.params['points'],
                    'is_polygon': hasattr(item,'polygon')}
            d['release_groups'][key]= rginfo

        # sort function times into order
        keys= []
        times=[]
        for key, f in profiling_util.func_timings.items():
            times.append(f['time'])
            keys.append(key)
            f['msec_per_call'] = 1000*f['time']/f['calls']
            f['% of total time'] = 100*f['time']/elapsed_time_sec

        # add in reverse timing order
        for n in np.argsort(np.asarray(times))[::-1]:
            d['function_timers'][keys[n]]= profiling_util.func_timings[keys[n]]

        # block timings in time order
        b = si.block_timers
        times = np.asarray( [item['time'] for key, item in b.items()])
        order = np.argsort(times)[::-1]
        for key in [list(b.keys())[i] for i in order]:
            l = f' {100*b[key]["time"]/elapsed_time_sec:5.1f}% {key} : calls {b[key]["calls"]:4d}, {time_util.seconds_to_pretty_duration_string(b[key]["time"])}'
            d['block_timings'].append(l)
        d['block_timings'].append(f'--- Total time {time_util.seconds_to_pretty_duration_string(elapsed_time_sec)}')


        # check numba code for SIMD
        if True:
            from numba.core import config
            d['numba_code_info'] = dict(signatures={},SMID_code = {},
                                        config={key:val for key, val in config.__dict__.items()
                                                    if not key.startswith('_') and type(val) in [None,int,str,float]
                                                }
                                        )
            for name, func in  numba_util.numba_func_info.items():
                if hasattr(func,'signatures') : # only code that has been compiled has a sig
                    sig = func.signatures
                    d['numba_code_info']['signatures'][name] = str(sig)
                    d['numba_code_info']['SMID_code'][name] = []
                    for nsig in range(len(sig)):
                        d['numba_code_info']['SMID_code'][name].append(numba_util.count_simd_intructions(func, sig=nsig))
                    pass
        return d

    def close(self):
        # close all instances, eg their files if not close etc
        for i in si.all_class_instance_pointers_iterator():
            try:
                i.close()

            except Exception as e:
                si.msg_logger.msg(f'Unexpected error closing class ="{ i.info["name"]}"', fatal_error= True, exception=e)

