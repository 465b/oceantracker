from typing import TypedDict

from oceantracker import definitions
from oceantracker.util.parameter_checking import ParameterCoordsChecker as PCC, ParameterListChecker as PLC

from oceantracker.util.parameter_checking import ParamValueChecker as PVC, ParameterTimeChecker as PTC

from oceantracker.util import class_importer_util, basic_util
from oceantracker.util.message_logger import  MessageLogger

from time import  perf_counter
from oceantracker.util import cord_transforms
import numpy as np
from oceantracker.util.scheduler import Scheduler
from oceantracker.particle_properties import particle_operations

# useful utility classes to enable auto complete
class _Object(object):  pass
class _SharedStruct():
    '''
    holds variables as clas attributes to enable auto complete hints
    and give iterators over these variables

    allows both  instance.backtracking and i['instance.backtracking']

    '''
    def __init__(self):
        # add  class variables in ._class_.__dict__
        # to instance __dict__ by adding attributes
        for key, item in self.__class__.__dict__.items():
            if not key.startswith('_'):
                setattr(self, key, item)
    def as_dict(self):
        d = dict()
        for key, item in self.__dict__.items():
            if not key.startswith('_'): d[key] = item
        return d

    def possible_values(self):
        d = []
        for key in self.__dict__.keys():
            if not key.startswith('_'): d.append(key) 
        return d
    def items(self): return self.as_dict().items()

    def __getitem__(self, name:str):
        return getattr(self,name)
    def __setitem__(self, name:str, value):
        setattr(self,name, value)

# default settings structure

class _DefaultSettings(_SharedStruct):

    root_output_dir=  PVC('root_output_dir', str, doc_str='base dir for all output files')
    add_date_to_run_output_dir =  PVC(False,bool, doc_str='Append the date to the output dir. name to help in keeping output from different runs separate' )
    output_file_base =    PVC('output_file_base', str,
                doc_str= 'The start/base of all output files and name of sub-dir of "root_output_dir" where output will be written' )
    time_step = PVC(3600., float, min=0.001, units='sec',doc_str='Time step in seconds for all cases' )
    screen_output_time_interval = PVC(3600., float, doc_str='Time in seconds between writing progress to the screen/log file' )
    backtracking =   PVC(False, bool, doc_str='Run model backwards in time')
    regrid_z_to_uniform_sigma_levels = PVC(True, bool,
                doc_str='much faster 3D runs by re-griding hydo-model fields in the z to uniform sigma levels on read, based on sigma most curve z_level profile. Some hydo-model are already uniform sigma, so this param is ignored, eg ROMS' )
    display_grid_at_start = PVC(False, bool,
                doc_str='Pause during strat up to plot the grid for checking using matplotlib, clicking om image will print a coord' )
    dev_debug_plots = PVC(False, bool,expert=True, doc_str='show any debug plot generated at give dbug_level, not for general use' )
    debug = PVC(False, bool, doc_str= 'more info on errors' )
    dev_debug_opt = PVC(0, int,expert=True,doc_str= 'does extra checks given by integer, not for general use' )
    minimum_total_water_depth = PVC(0.25, float, min=0.0, units='m', doc_str='Min. water depth used to decide if stranded by tide and which are dry cells to block particles from entering' )
                #'write_output_files =     PVC(True,  bool, doc_str='Set to False if no output files are to be written, eg. for output sent to web' )
    write_dry_cell_flag = PVC(True, bool,
                doc_str='Write dry cell flag to all cells when writing particle tracks, which can be used to show dry cells on plots,may create large grid file, currently cannot be used with nested grids ' )
    max_run_duration = PVC(definitions.max_timedelta_in_seconds, float,min=.00001,units='sec',
                           doc_str='Useful in testing setup with shorter runs, as normally run duration is determined from release groups. This  limits the maximum duration in seconds of model runs.' )  # limit all cases to this duration
    max_particles = PVC(10**10, int, min=1, doc_str='Maximum number of particles to release, useful to restrict if splitting particles' )  # limit all cases to this number
    processors = PVC(None, int, min=1,doc_str='number of processors used, if > 1 then cases in the case_list to run in parallel' )
    max_warnings = PVC(50,    int, min=0,doc_str='Number of warnings stored and written to output, useful in reducing file size when there are warnings at many time steps' )  # dont record more that this number of warnings, to keep caseInfo.json finite
    use_random_seed = PVC(False,  bool,doc_str='Makes results reproducible, only use for testing developments give the same results!', expert=True )
    NUMBA_function_cache_size = PVC(4048, int, min=128, expert=True,
                                    doc_str='Size of memory cache for compiled numba functions in kB' )
    NUMBA_cache_code = PVC(False, bool, expert=True,
                           doc_str='Speeds start-up by caching complied Numba code on disk in root output dir. Can ignore warning/bug from numba "UserWarning: Inspection disabled for cached code..."' )
    multiprocessing_case_start_delay = PVC(0., float, min=0.,expert=True,
                                           doc_str='Delay start of each sucessive case run parallel, to reduce congestion reading first hydo-model file' )  # which large numbers of case, sometimes locks up at start al reading same file, so ad delay
    EPSG_code_meters_grid = PVC(None, int,
                doc_str='If hydro-model has lon_lat coords, then grid is converted to this meters system. For codes see https://epsg.io/. eg EPSG for NZ Transverse Mercator use 2193. Default grid is UTM' )
    write_tracks = PVC(True, bool, doc_str='Flag if "True" will write particle tracks to disk. For large runs and statistics done on the fly, is normally set to False to reduce output volumes' )
    user_note = PVC('No user note', str, doc_str='Any run note to store in case info file' )
    z0 = PVC(0.005, float, units='m', doc_str='Bottom roughness, used for tolerance and log layer calcs. ', min=0.0001 )  # default bottom roughness
    water_density = PVC(1025., float, units='kg/m^3', doc_str='Water density , default is seawater, an example of use is in calculating friction velocity from bottom stress, ', min=900. )
    use_open_boundary = PVC(False, bool, doc_str='Allow particles to leave open boundary, only works if open boundary nodes  can be read or inferred from hydro-model, current schism using hgrid file, and inferred for structed grids like ROMS ' )
    block_dry_cells = PVC(True, bool, doc_str='Block particles moving from wet to dry cells, ie. treat dry cells as if they are part of the lateral boundary' )
    use_geographic_coords = PVC(False, bool,
                          doc_str='Used geographic coordniated for inputs and outputs ( lon, lat_), normally auto detected based in hindcast coords (if True and hindcast already geographic coords, then reader must have EPGS code',
                                expert=True)
    use_A_Z_profile = PVC(True, bool,
                doc_str='Use the hydro-model bottom_stress variable for friction velocity calculation , where it is needed for resuspension, if variable is in hindcast files')
    use_bottom_stress = PVC(True, bool,
                doc_str='Use hydro models bottom_stress variable for friction velocity calculation, if mapped variable is in files. Friction velocity is used in resuspension')
    use_dispersion = PVC(True, bool,
                doc_str='Include random walk, allows it to be turned off if needed for applications like Lagrangian coherent structures')
    use_resuspension = PVC(True, bool,
                doc_str='Allow particles to resuspend')
    NCDF_time_chunk = PVC(24, int, min=1,expert=True,
                 doc_str='Used when writing time series to netcdf output, is number of time steps per time chunk in the netcdf file')

    particle_buffer_initial_size= PVC(500_000, int, min=1, doc_str='Starting size of particle property memory buffer. This expands by particle_buffer_chunk_size as needed', expert=True)
    particle_buffer_chunk_size = PVC(500_000, int, min=1, doc_str='How much particle property memory buffer sizes are increased by when they are full',
                                     expert=True)
        #  #'loops_over_hindcast =  PVC(0, int, min=0 )  #, not implemented yet,  artifically extend run by rerun from hindcast from start, given number of times
        # profiler = PVC('oceantracker', str, possible_values=available_profile_types,
        #                 doc_str='in development- Default oceantracker profiler, writes timings of decorated methods/functions to run/case_info file use of other profilers in development and requires additional installed modules ' )
        # 'debug_level =               PVC(0, int,min=0, max=10, doc_str='Gives  diferent levels of debug, in development' )


# blocks that make up parts of shared info
class _ClassRoles(_SharedStruct):
    release_groups =[]
    fields = []  # user fields calculated from other fields  on reading
    particle_properties =  [] # user added particle properties, eg DistanceTraveled
    velocity_modifiers = []  # user added velocity effects, eg TerminalVelocity
    trajectory_modifiers = []  # change particle paths, eg. re-suspension
    particle_statistics = []  # heat map inside polygon statistics calculated on the fly
    particle_concentrations = []  # writes concentration of particles and other properties calculated on the fly.   files ,eg PolygonEntryExit
    nested_readers = []
    event_loggers =  [] # writes events files ,eg PolygonEntryExit
    time_varying_info = [] # particle info,eg. time,or  tide at at tide gauge, core example is particle time

class _CoreClassRoles(_SharedStruct):
    reader = None
    interpolator = None
    #todo below beter as None
    solver = {}
    field_group_manager = {}
    particle_group_manager = {}
    tracks_writer = {}
    dispersion = {}
    tidal_stranding = {}
    resuspension = {}
    integrated_model = None # this is here as there can be only one at a time
from dataclasses import dataclass

class _ParticleStatusFlags(_SharedStruct):
    '''Particle status flags mapped to integer values'''
    unknown  = -20
    bad_cord = -16
    cell_search_failed=  -15
    notReleased = -10
    dead = -5
    outside_open_boundary =-2
    stationary = 0
    stranded_by_tide = 3
    on_bottom = 6
    moving =  10


class _VerticalGridTypes(_SharedStruct):
    '''Particle status flags mapped to integer values'''
    Slayer  = 'Slayer'
    LSC = 'LSC'
    Sigma = 'Sigma'
    Zfixed = 'Zfixed'

class _RunInfo(_SharedStruct):
    is3D_run = None
    backtracking =None
    vector_components = None
    model_direction = None
    free_wheeling = None
    start_time = None
    end_time = None
    current_model_time = None
    current_model_date = None
    current_model_time_step = 0
    times = None
    duration = None
    run_output_dir = None
    output_file_base = None
    caseID = None
    time_of_nominal_first_occurrence = None
    total_alive_particles = 0
    time_steps_completed = 0
    hindcast_start_time = None
    hindcast_end_time = None

class _UseFullInfo(_SharedStruct):
    # default reader classes used by auto-detection of file type
    large_float = 1.0E50

# Shared class, build using the above
#------------------------------------------------------------
class _SharedInfoClass():
    """Gives access to shared variables, classes and some utility methods. Attributes are:

        eg SharedInfo.particle_properties is a dictionary of named particle property instances

    """
    settings = _DefaultSettings() # will be overwritten with actual values by case runner
    class_roles = _ClassRoles()
    core_class_roles = _CoreClassRoles()
    particle_operations = particle_operations
    default_settings = _DefaultSettings()
    particle_status_flags = _ParticleStatusFlags() # need to be instances to allow particle_status_flags[key] form
    vertical_grid_types = _VerticalGridTypes()
    run_info  = _RunInfo()
    hindcast_info = None
    msg_logger = MessageLogger()
    block_timers={}
    classes = {}  # todo deprecated
    info = _UseFullInfo

    # list of params only setable in base case
    base_case_only_params = ['root_output_dir', 'output_file_base', 'processors', 'max_warnings',
                                  'multi_processing_method', 'multiprocessing_case_start_delay',
                                  'backtracking', 'add_date_to_run_output_dir', 'debug', 'use_random_seed']
    base_case_only_params += [key for key in default_settings.possible_values() if 'numba' in key.lower()]


    def __init__(self):

        self.default_polygon_dict_params = {'user_polygonID': PVC(0, int, min=0),
                                       'name': PVC(None, str),
                                       'points': PCC(None, is_required=True, doc_str='Points making up the polygon as, N by 2 or 3 list of locations where particles are released. eg for 2D ``[[25,10],[23,2],....]``, must be convertible into N by 2 or 3 numpy array')
                                       }
        pass
    def add_settings(self, settings):
        for key in self.settings.possible_values():
            setattr(self.settings, key, settings[key])

    def _setup(self):
        # this allows shared info to make a class importer when needed
        self.msg_logger.set_screen_tag('Prelim')
        self.msg_logger.reset()
        self._class_importer = class_importer_util.ClassImporter(self.msg_logger)

        # empty out roles and core roles in case of rerunning and shared info import only happens once
        for role in self.core_class_roles.as_dict().keys():
            setattr(self.core_class_roles, role, None)
        for role in self. class_roles.as_dict().keys():
            setattr(self.class_roles, role, {})


    def add_class(self,class_role,params={}, default_classID=None,caller=None,crumbs ='', initialize=False,check_for_unknown_keys=True, **kwargs):
        #todo get rid in initialize
        ml = self.msg_logger
        crumbs += f'Adding class {class_role}>'

        if class_role=='fields':
            ml.msg('Cannot use si.add_class() method to add fields',
                   hint='Use add_reader_field(name,  params) or si.add_custom_field(name, params, default_classID=None)',
                   fatal_error=True)


        if params is None: params ={}
        if type(params) != dict :
            ml.msg(f'Params must be a dictionary', hint= f'Got type {str(type(params))}',
                        error=True, crumbs=crumbs,
                        check_for_unknown_keys=check_for_unknown_keys, caller=caller)
            return None

        params= dict(params,**kwargs) # join params and kwargs

        if class_role in self.core_class_roles.possible_values():
            #core  roles
            params['name'] = None
            i =  self._class_importer.make_class_instance_from_params(class_role, params, default_classID=default_classID, crumbs=crumbs,
                                          check_for_unknown_keys=check_for_unknown_keys, caller=caller, initialize=initialize)
            i.info['instanceID'] = 0
            self.core_class_roles[class_role] = i

        elif class_role in self.class_roles.possible_values():
            #other roles
            instanceID= len(self.class_roles[class_role])
            i = self._class_importer.make_class_instance_from_params(class_role, params, default_classID=default_classID,
                                                                     crumbs=crumbs, caller=caller,initialize=initialize)
            i.info['instanceID'] = instanceID
            if params['name'] is None:
                # if no name in params or default param
                params['name'] = f'{class_role}_{instanceID:04d}'

            self.class_roles[class_role][params['name']] = i

        else:
            ml.msg(f'Unknown class role {class_role}', hint=f'Must be one of core_class_roles {str(self.core_class_roles.possible_values())} or other roles {str(self.class_roles.possible_values())}',
                   error=True, crumbs=crumbs, caller=caller)
            return None

        # add classes required by this class
        i.add_required_classes_and_settings(self.settings, self.run_builder['reader_builder'], self.msg_logger)
        i.si = self # for alternative acess to shared info
        return i

    # wrapers for adding fields
    def add_reader_field(self, name, params):
        return self.core_class_roles.field_group_manager.add_reader_field(name, params)
    def add_custom_field(self, name, params, default_classID=None):
        return self.core_class_roles.field_group_manager.add_custom_field(name, params, default_classID=default_classID)


    def _all_class_instance_pointers_iterator(self):
        # build list of all points for iteration, eg in calling all close methods
        p = []
        for role in self.core_class_roles.possible_values():
            i = getattr(self.core_class_roles, role)
            if i is not None: p.append(i)

        for role in self.class_roles.possible_values():
            r = getattr(self.class_roles, role)
            if r is not None:
                for key, i in r.items():
                    if i is not None:  p.append(i)
        return p
    def block_timer(self,name,t0):
        b = self.block_timers
        if name not in b:
            b[name] = dict(time=0.,calls=0)
        b[name]['time'] += perf_counter()-t0
        b[name]['calls'] += 1


# make the instance used throughout code
shared_info = _SharedInfoClass()
pass