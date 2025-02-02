import numpy as np
from oceantracker.util.parameter_base_class import ParameterBaseClass
from oceantracker.util import time_util
from oceantracker.util.parameter_checking import ParamValueChecker as PVC, ParameterListChecker as PLC, ParameterCoordsChecker as PCC
from numba import njit
from oceantracker.util.numba_util import njitOT
from oceantracker.common_info_default_param_dict_templates import large_float

class PointRelease(ParameterBaseClass):
    # releases particles at fixed points, inside optional radius
    # add checks to see if points inside domain and dry if released in a radius
    #todo make a parent release base class

    def __init__(self):
        # set up info/attributes
        super().__init__()
        self.add_default_params({
                                 'points':          PCC(None,is_required=True, doc_str='A N by 2 or 3 list or numpy array of locations where particles are released. eg for 2D [[25,10],[23,2],....] '),
                                 'release_radius':  PVC(0., float, min= 0., doc_str= 'Particles are released from random locations in circle of given radius around each point.'),
                                 'pulse_size' :     PVC(1, int, min=1, doc_str= 'Number of particles released in a single pulse, this number is released every release_interval.'),
                                 'release_interval':PVC(0., float, min =0.,units='sec', doc_str= 'Time interval between released pulses. To release at only one time use release_interval=0.'),
                                 'release_start_date': PVC(None, 'iso8601date', doc_str='Must be an ISO date as string eg. "2017-01-01T00:30:00" '),
                                   # to do add ability to release on set dates/times 'release_dates': PLC([], 'iso8601date'),
                                 'release_duration': PVC(None, float,min=0.,
                                                    doc_str='Time in seconds particles are released for after they start being released, ie releases stop this time after first release.,an alternative to using "release_end_date"' ),
                                'release_end_date': PVC(None, 'iso8601date', doc_str='Date to stop releasing particles, ignored if release_duration give, must be an ISO date as string eg. "2017-01-01T00:30:00" '),
                                'max_age': PVC(None,float,min=1.,
                                                    doc_str='Particles older than this age in seconds are culled,ie. status=dead, and removed from computation, very useful in reducing run time'),
                                 'user_release_groupID' : PVC(0,int, doc_str= 'User given ID number for this group, held by each particle. This may differ from internally uses release_group_ID.'),
                                 'user_release_group_name' : PVC(None,str,doc_str= 'User given name/label to attached to this release groups to make it easier to distinguish.'),
                                 'allow_release_in_dry_cells': PVC(False, bool,
                                              doc_str='Allow releases in cells which are currently dry, ie. either permanently dry or temporarily dry due to the tide.'),
                                 'z_range': PLC([],[float, int], min_length=2, obsolete='use z_min and/or z_max'),

                                'z_min': PVC(-large_float, [float, int],doc_str='min/ deepest z value to release for to randomly release in 3D, overrides any given release z value'),
                                'z_max': PVC( large_float, [float, int], doc_str='max/ highest z vale release for to randomly release in 3D, overrides any given release z value'),
                                'release_offset_above_bottom': PVC(False, [float, int], min= 0., doc_str=' 3D release particles at fixed give height above the bottom at the release location ', units='m'),
                                #'water_depth_min': PVC(None, float,doc_str='min water depth to release in, useful for releases with a depth rage, eg larvae from inter-tidal shellfish', units='m'),
                                #'water_depth_max': PVC(None, float, doc_str='max water depth to release in', units='m'),

            #Todo implement release group particle with different parameters, eg { 'oxygen' : {'decay_rate: 0.01, 'initial_value': 5.}
                                'max_cycles_to_find_release_points': PVC(200, int, min=100, doc_str='Maximum number of cycles to search for acceptable release points, ie. inside domain, polygon etc '),
                                 })
        self.class_doc(description= 'Release particles at 1 or more given locations. Pulse_size particles are released every release_interval. All these particles are tagged as a single release_group.')


    def initial_setup(self):
        # must be called after unpack_x0
        # tidy up parameters to make them numpy arrays with first dimension equal to number of locations

        info=self.info
        info['points'] =   np.array(self.params['points']).astype(np.float64)

        info['number_released'] = 0 # count of particles released in this group
        info['pulse_count'] = 0
        info['release_type'] = 'point'

    def set_up_release_times(self):
        # get release times based on release_start_date, duration
        params = self.params
        info = self.info
        si = self.shared_info
        ml = self.msg_logger

        hindcast_start, hindcast_end  =  si.classes['field_group_manager'].get_hindcast_start_end_times()

        model_time_step = si.settings['time_step']

        self.info['release_info'] ={'first_release_date': None, 'last_release_date':None,
                                    'last_time_alive':None}
        # short cut
        release_info =self.info['release_info']

        if params['release_start_date'] is None:
            # no user start date so use  model runs' start date
            time_start= hindcast_start if not si.backtracking else hindcast_end
        else:
            # user given start date
            time_start = time_util.isostr_to_seconds(params['release_start_date'])

        # now check if start in range
        n_groups_so_far =len(si.classes['release_groups'])
        if not hindcast_start <= time_start <= hindcast_end:
            si.msg_logger.msg('Release group= ' + str(n_groups_so_far + 1) + ', name= ' + self.info['name'] + ',  parameter release_start_time is ' +
                                    time_util.seconds_to_isostr(time_start)
                              + '  is outside hindcast range ' + time_util.seconds_to_isostr(hindcast_start)
                                    + ' to ' + time_util.seconds_to_isostr(hindcast_end), warning=True)

        # set max age of particles
        release_interval = model_time_step if params['release_interval'] is None else params['release_interval']

        # world out release times
        if release_interval == 0.:
            time_end = time_start
        elif self.params['release_duration'] is not None:
            time_end = time_start + si.model_direction*self.params['release_duration']

        elif self.params['release_end_date'] is not None:
            time_end = time_util.isostr_to_seconds(self.params['release_end_date'])
        else:
            # default is limit of hindcast
            time_end = hindcast_start if si.backtracking else hindcast_end

        # get time steps for release in a dow safe way
        model_time_step = si.settings['time_step']


        # get release times within the hindcast
        if abs(time_end-time_start) < model_time_step:
            # have only one release
            release_info['release_times'] = np.asarray(time_start)
        else:
            release_info['release_times'] = time_start + np.arange(0., abs(time_end-time_start),release_interval )*si.model_direction

        # trim releases to be within hindcast
        sel = np.logical_and( release_info['release_times'] >= hindcast_start,  release_info['release_times']  <= hindcast_end)
        release_info['release_times'] = release_info['release_times'][sel]

        if release_info['release_times'].size ==0:
            ml.msg(f'No release times in range of hydro-model for release_group {info["instanceID"]:2d}, ',
                   fatal_error=True,
                   hint=' Check hydro-model date range and release dates  ')

        # get time steps when released, used to determine when to release
        release_info['release_time_steps'] =  np.round(( release_info['release_times']- hindcast_start)/model_time_step).astype(np.int32)

        # find last time partiles alive
        max_age = 1.0E30 if params['max_age'] is None else params['max_age']
        release_info['last_time_alive'] =  release_info['release_times'][-1] + si.model_direction*max_age
        release_info['last_time_alive'] =  min(max(hindcast_start,release_info['last_time_alive']),hindcast_end) # trim to limits of hind cast

        # useful info
        release_info['first_release_time'] = release_info['release_times'][0]
        release_info['last_release_time'] = release_info['release_times'][-1]

        release_info['release_dates'] = time_util.seconds_to_datetime64(release_info['release_times'])
        release_info['first_release_date'] = time_util.seconds_to_datetime64(release_info['first_release_time'])
        release_info['last_release_date'] = time_util.seconds_to_datetime64(release_info['last_release_time'])

        # index of release the  times to be released next
        release_info['index_of_next_release'] =  0

        if not   hindcast_start <= release_info['first_release_time'] <= hindcast_end :
            ml.msg(f'Release group "{info["name"]}" >  start time {time_util.seconds_to_isostr(release_info["first_release_time"])}  is outside the range of hydro-model times for release_group instance #{info["instanceID"]:2d}, ',
                   fatal_error=True,hint=f' Check release start time is in hydro-model  range of  {time_util.seconds_to_isostr(hindcast_start)}  to {time_util.seconds_to_isostr(hindcast_start)} ')

    def get_release_locations(self, time_sec):
        # set up full set of release locations inside  polygons
        si = self.shared_info
        info= self.info
        params=self.params

        n_required = self.get_number_required()

        x0           = np.full((0, info['points'].shape[1]), 0.,dtype=np.float64, order='C')
        n_cell0 = np.full((0,), -1, dtype=np.int32)
        hydro_model_gridID0 = np.full((0,), 0, dtype=np.int8)
        bc_cords0 = np.full((0,3), 0, dtype=np.float64)

        count = 0
        n_found = 0

        while x0.shape[0] < n_required:
            # get 2D release candidates
            x_guess = self.get_release_location_candidates()

            x_guess, n_cell_guess, bc_cords, hydro_model_gridID = self.check_potential_release_locations_in_bounds(x_guess)

            if x_guess.shape[0] > 0:
                is_ok = self.filter_release_points(x_guess,  time_sec= time_sec, n_cell=n_cell_guess,
                                                   bc_cords=bc_cords, hydro_model_gridID=hydro_model_gridID)
                # if any ok then add to list
                n_found += x_guess.shape[0]
                x0          = np.concatenate((x0, x_guess[is_ok,:]), axis =0)
                n_cell0 = np.concatenate((n_cell0, n_cell_guess[is_ok],))
                bc_cords0 = np.concatenate((bc_cords0, bc_cords[is_ok,:],))
                hydro_model_gridID0= np.concatenate((hydro_model_gridID0, hydro_model_gridID[is_ok],))

            # allow max_cycles_to_find_release_points cycles to find points
            count += 1
            if count > self.params["max_cycles_to_find_release_points"]: break

        if n_found < n_required:
            si.msg_logger.msg(f'Release group-"{self.info["name"]}", only found {n_found} of {n_required} required points inside domain after {self.params["max_cycles_to_find_release_points"]} cycles',
                              warning=True,
                           hint=f'Maybe, release points outside the domain?, or hydro-model grid and release points use different coordinate systems?? or increase parameter  max_cycles_to_find_release_points, current value = {self.params["max_cycles_to_find_release_points"]:3}' )
            n_required = n_found #


        # trim initial location, cell  etc to required number
        x0 = x0[:n_required, :]
        n_cell0 = n_cell0[:n_required]
        bc_cords0 = bc_cords0[:n_required,:]
        hydro_model_gridID0= hydro_model_gridID0[:n_required]

        n = x0.shape[0]
        IDrelease_group = self.info['instanceID']
        IDpulse = info['pulse_count']
        info['pulse_count'] += 1
        user_release_groupID = self.params['user_release_groupID']

        info['number_released'] += n  # count number released in this group

        if si.is3D_run:

            if params['z_min'] >= params['z_max']:
                si.msg_logger.msg(f'Release group-"{self.info["name"]}", zmin >= zmax, (zmin,zmax) =({info["z_range"][0]:.3e}, {info["z_range"][1]:.3e}) ',fatal_error=True)

            fgm = si.classes['field_group_manager']
            water_depth = fgm.interp_named_field_at_given_locations_and_time('water_depth', x0, time_sec=None,     n_cell=n_cell0, bc_cords=bc_cords0, hydro_model_gridID=hydro_model_gridID)
            tide        = fgm.interp_named_field_at_given_locations_and_time('tide',        x0, time_sec=time_sec, n_cell=n_cell0, bc_cords=bc_cords0, hydro_model_gridID=hydro_model_gridID)
            if x0.shape[1] == 2:
                # expand x0 to 3D if needed
                x0 = np.concatenate((x0, np.zeros((x0.shape[0],1), dtype=x0.dtype)), axis=1)

            x0 = self.get_z_release_in_depth_range(x0,params['z_min'],params['z_max'],  water_depth, tide)

        return x0, IDrelease_group, IDpulse, user_release_groupID, n_cell0, bc_cords0, hydro_model_gridID0

    @staticmethod
    @njitOT
    def get_z_release_in_depth_range(x,z_min, z_max, water_depth,tide):
        # get random release within zrange within bounds of water depth and tide

        zr =  np.full((2,),0.)

        for n in range(x.shape[0]):
            # get mean depth of triangle by summing

            z1  = max(water_depth[n], z_min)
            z2 = min(tide[n], z_max)
            a = np.random.uniform(z1, z2, size=1)[0]
            x[n, 2] = np.random.uniform(z1, z2, size=1)[0]

        return x




    def get_number_required(self):
        return self.params['pulse_size']*self.info['points'].shape[0]

    def get_release_location_candidates(self):
        si = self.shared_info
        x = np.repeat(self.info['points'], self.params['pulse_size'], axis=0)

        if self.params['release_radius']> 0.:
            rr = abs(float(self.params['release_radius']))
            n = x.shape[0]
            rr = np.repeat(rr, n, axis=0)
            r = np.random.random((n,)) * rr * np.exp(1.0j * np.random.random((n,)) * 2.0 * np.pi)
            r = r.reshape((-1, 1))
            x[:, :2] += np.hstack((np.real(r), np.imag(r)))

        return x


    def filter_release_points(self, x, time_sec= None, n_cell=None,bc_cords=None, hydro_model_gridID=None):
        # user can filter release points by inheritance of this class and overriding this method
        return np.full((x.shape[0],), True)

    def check_potential_release_locations_in_bounds(self, x):
        si= self.shared_info
        # use KD tree to find points those outside model domain
        fgm = si.classes['field_group_manager']
        is_inside, n_cell, bc, hydro_model_gridID  = fgm.are_points_inside_domain(x, self.params['allow_release_in_dry_cells'])

        # keep those inside domain
        x = x[is_inside, :]
        n_cell = n_cell[is_inside]
        bc = bc[is_inside,:]
        hydro_model_gridID= hydro_model_gridID[is_inside]
        return x, n_cell, bc, hydro_model_gridID

