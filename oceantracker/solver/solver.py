from time import perf_counter
from numba import njit

from oceantracker.util import time_util
from datetime import datetime

from oceantracker.particle_properties.util import particle_operations_util
from oceantracker.util.parameter_base_class import ParameterBaseClass
from oceantracker.util.parameter_checking import ParamDictValueChecker as PVC



class Solver(ParameterBaseClass):
    #  does particle tracking solution as class to allow multi processing

    def __init__(self):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults

        self.add_default_params({ 'screen_output_step_count':   PVC(1, int),
                                  'RK_order':                   PVC(4, int, possible_values=[1, 2, 4]),
                                  'n_sub_steps':                PVC(1, int, min=1),
                                  'name':                       PVC('solver',str) })

    def initialize(self):

        self.code_timer.start('solver_initialization')
        si = self.shared_info
        si.classes['particle_group_manager'].create_particle_property('manual_update', dict(name='v_temp', vector_dim=si.classes['particle_properties']['x'].num_vector_dimensions(), write=False))

        # set up working space for RK stesp to impriove L3 cache performance
        self.code_timer.stop('solver_initialization')

    def check_requirements(self):

        msg_list = self.check_class_required_fields_properties_grid_vars_and_3D(
            required_fields=['water_velocity'],
            required_props=['x','status', 'x_last_good', 'particle_velocity', 'v_temp'],
            required_grid_vars=[])


        return msg_list

    def initialize_run(self):
        si = self.shared_info
        info = self.info
        # set up particle velocity working space for solver
        info['n_time_steps_completed'] = 0
        info['total_num_particles_moving'] = 0

    def solve_for_data_in_buffer(self, nb0, num_in_buffer, nt0):
        # solve fro dat in buffer
        si = self.shared_info
        info = self.info
        p   = si.classes['particle_group_manager']
        computation_started = datetime.now()


        for nb in range(nb0,nb0 + num_in_buffer-1): # one less step as last step is initial condition for next block

            t_hindcast = si.grid['time'][nb]  # make time exactly that of grid

            # do sub steps with hindcast model step
            for ns in range(self.params['n_sub_steps']):
                # round start time to nearest hindcast step
                t0_step = perf_counter()
                t1 = t_hindcast + ns*si.model_substep_timestep*si.model_direction

                p.release_particles(nb, t1)

                # post release write time varying part properties and statistics
                if si.write_output_files and si.write_tracks and info['n_time_steps_completed'] % si.classes['tracks_writer'].params['output_step_count'] == 0:
                    p.write_time_varying_prop_and_data()

                # update and write stats
                self.code_timer.start('on_the_fly_statistics')
                for s in si.class_list_interators['particle_statistics']['all'].values():
                    if abs(t1 - s.info['time_last_stats_recorded']) >= s.params['calculation_interval']:
                        s.update(time=t1)
                self.code_timer.stop('on_the_fly_statistics')

                # update triangle concentrations
                self.code_timer.start('particle_concentrations')
                for s in si.class_list_interators['particle_concentrations']['all'].values():
                    if abs(t1 - s.info['time_last_stats_recorded']) >= s.params['calculation_interval']:
                        s.update(nb,t1)
                self.code_timer.stop('particle_concentrations')

                # write events
                self.code_timer.start('event_logging')
                for e in  si.class_list_interators['event_loggers']['all'].values():
                    e.update(time=t1)
                self.code_timer.stop('event_logging')

                # print screen data
                if (info['n_time_steps_completed']  + ns) % self.params['screen_output_step_count'] == 0:
                    self.screen_output(info['n_time_steps_completed'] , nt0, nb0, nb, ns, t1, t0_step,computation_started)

                #  Main integration step
                #  --------------------------------------
                self.code_timer.start('integration_step')
                #  --------------------------------------

                moving = self.integration_step(nb,t1)
                #--------------------------------------
                self.code_timer.stop('integration_step')
                info['total_num_particles_moving'] += moving.shape[0]

                t2 = t1 + si.model_substep_timestep*si.model_direction

                si.classes['dispersion'].update(nb, t2, moving)

                self.post_step_bookeeping(si.model_start_time,nb,info['n_time_steps_completed'], t2, moving)

                p.kill_old_particles(t2)

                if si.compact_mode: p.remove_dead_particles_from_memory()

                info['n_time_steps_completed']  += 1

                if abs(t2 - si.model_start_time) > si.model_duration:  break

        return info['n_time_steps_completed'], t2

    def integration_step(self, nb, t):
        # single step in particle tracking, t is time in seconds, is_moving are indcies of moving particles
        # this is done inplace directly operation on the particle properties
        # nb is buffer offset
        si = self.shared_info
        RK_order =self.params['RK_order']
        f = si.classes['field_group_manager']
        part_prop =  si.classes['particle_properties']

        # note here subStep_time_step has sign of forwards/backwards
        dt = si.model_substep_timestep*si.model_direction
        # set up views of  working variable space
        x1      = part_prop['x_last_good'].data
        x2      = part_prop['x'].data
        v       = part_prop['particle_velocity'].data
        v_temp  = part_prop['v_temp'].data  # temp vel from interp at each RK substeps
        is_moving = part_prop['status'].compare_all_to_a_value('eq', si.particle_status_flags['moving'], out = self.get_particle_index_buffer())

        # this makes x1, ['x_last_good']  at start of new integration step for moving particles, allowing updates to x2 ['x']
        particle_operations_util.copy(x1, x2, is_moving)

        #  step 1 from current location and time
        f.setup_interp_time_step(nb, t,  x1, is_moving)

        if RK_order==1:
            self.update_particle_velocity(t,v, is_moving) # put vel into permanent place
            self.euler_substep(x2, x1, v, dt, is_moving)
            return is_moving

        self.update_particle_velocity(t,v_temp, is_moving)  # velocity in temp place

        # do first half step location from RK1 to update values
        self.euler_substep(x2, x1, v_temp, dt / 2., is_moving)

        # accumulate RK velocity to reduce space taken by temporary working variables
        particle_operations_util.copy(v, v_temp, is_moving, scale=1.0 / 6.0) # vel at start of step

        # step 2, get improved half step velocity
        t2=t + 0.5*dt
        f.setup_interp_time_step(nb, t2, x2, is_moving)

        if RK_order==2:
            self.update_particle_velocity(t2,v, is_moving)
            self.euler_substep(x2, x1, v, dt, is_moving)
            return is_moving

        self.update_particle_velocity(t2, v_temp, is_moving)

        # step 3, a second half step
        self.euler_substep(x2, x1, v_temp, dt / 2., is_moving)  # improve half step position
        particle_operations_util.add_to(v, v_temp, is_moving, scale=2.0 / 6.0)  # next accumulation of velocity step 2

        t2 = t + 0.5 * dt
        f.setup_interp_time_step(nb, t2, x2, is_moving)
        self.update_particle_velocity(t2,v_temp , is_moving)  # v3, better velocity at half step

        self.euler_substep(x2, x1, v_temp, dt, is_moving)  # improve half step position values
        particle_operations_util.add_to(v, v_temp, is_moving, scale=2.0 / 6.0)  # next accumulation of velocity from step 3

        # step 4, full step
        t2 = t + dt
        f.setup_interp_time_step(nb, t2,  x2, is_moving)
        self.update_particle_velocity( t2, v_temp, is_moving)  # full step for v4

        particle_operations_util.add_to(v, v_temp, is_moving, scale=1.0 / 6.0)  # last accumulation of velocity for v4

        # final location using  accumulation in "v"
        # below is emulated by accumulation above of
        #  v = (v1 + 2.0 * (v2 + v3) + v4) /6
        #  x2 = x1 + v*dt
        self.euler_substep(x2, x1, v, dt, is_moving)  # set final location directly to particle x property

        return  is_moving

    def update_particle_velocity(self, t, v, active):
        # get water velocity plus any particle affects on velocity in particle_velocity, eg terminal velocity
        si= self.shared_info
        si.classes['field_group_manager'].interp_named_field_at_particle_locations('water_velocity', active, output=v)

        # any effects of super-imposable modifiers of the velocity, eg terminal  velocity in 3D to be added to w component
        for key, vm in si.classes['velocity_modifiers'].items():
            v = vm.modify_velocity(v, t, active)

    @staticmethod
    @njit
    def euler_substep(xnew, xold, velocity, dt, active):
        # do euler substep, xnew = xold+velocity*dt for active particles
        for n in active:
            for m in range(xold.shape[1]):
                xnew[n, m] = xold[n, m] + velocity[n, m] * dt
                    
    def post_step_bookeeping(self,t, nb, n_steps, t2, moving):
        # do dispersion, modify trajectories,
        # do strandings etc to change particle status
        # update part prop, eg  interp mapped reader fields to particle locations
        self.code_timer.start('post_step_bookeeping')
        si = self.shared_info
        pm = si.classes['particle_group_manager'] # internal short cuts
        fgm = si.classes['field_group_manager']
        part_prop  =  si.classes['particle_properties']

        # after dispersion some may be outside, update search cell status to see which are now outside domain etc
        fgm.setup_interp_time_step(nb, t2, part_prop['x'].data, moving)

        # user particle movements, eg resupension for all particles
        # re-find alive particles after above movements
        sel = part_prop['status'].compare_all_to_a_value('gteq', si.particle_status_flags['frozen'], out=self.get_particle_index_buffer())
        for i in si.class_list_interators['trajectory_modifiers']['all'].values():
            i.update(nb, t2, sel)

        # after moves, update search cell status to see which are now outside domain etc
        fgm.setup_interp_time_step(nb, t2, part_prop['x'].data, sel)

        # now all  particle movements complete after trajectory changes, move backs, update cell and bc cords for latest locations, update particle properties
        sel = part_prop['status'].compare_all_to_a_value('gteq', si.particle_status_flags['frozen'], out=self.get_particle_index_buffer())
        # now update part prop with good cell and bc cords, eg  interp water_depth
        pm.update_PartProp(t2, sel)

        # do any status only changes,
        # eg total water depth used for tidal stranding must be up to date
        # if 'total_water_depth' in part_prop :
        #     self.tidal_stranding_from_total_water_depth(part_prop['total_water_depth'].data,
        #                                             si.particle_status_flags['frozen'],
        #                                             si.particle_status_flags['stranded_by_tide'],
        #                                             si.particle_status_flags['moving'],
        #                                             si.minimum_total_water_depth,
        #                                             sel,
        #                                             part_prop['status'].data)

        # temp insert to test and include old stranding
        if  si.grid['is_dry_cell'] is not None:
            # un-strand those now stranded, later reset if still stranded
            currently_stranded = part_prop['status'].find_subset_where(
                sel, 'eq', si.particle_status_flags['stranded_by_tide'],
                out = self.get_particle_subset_buffer())
            part_prop['status'].set_values(si.particle_status_flags['moving'], currently_stranded)

            # now move back and strand all currently in dry cells, block if cell dry at this or next hindcast time step
            is_dry = self.is_in_dry_cell(
                nb, si.grid['is_dry_cell'], part_prop['n_cell'].dataInBufferPtr(), 
                sel, self.get_particle_subset_buffer())

            if is_dry.shape[0] > 0:
                part_prop['status'].set_values(si.particle_status_flags['stranded_by_tide'], is_dry)
                self.return_to_last_position(nb, t2, is_dry)
        
        self.code_timer.stop('post_step_bookeeping')

    @staticmethod
    @njit
    def tidal_stranding_from_total_water_depth(total_water_depth, status_frozen, status_stranded,status_moving, min_water_depth, sel, status):
        # look at all particles in buffer to check total water depth < min_water_depth
        for n in sel:
            if status[n] >= status_frozen:
                if total_water_depth[n] < min_water_depth:
                    status[n] = status_stranded
                elif status[n] == status_stranded:
                    # unstrand if already stranded, if status is on bottom,  remains as is
                    status[n] = status_moving
    
    
    @staticmethod
    @njit
    def is_in_dry_cell(nb,dry_cell_buffer,n_cell, active, out):
        # return a view of out with indices of those in dry cells
        nfound =0
        for n in active:
            nc = n_cell[n]
            # in dry cell if cell dry at this or next time step in buffer
            # todo should this be time interpolation between time steps and dry if > 0.5
            if dry_cell_buffer[nb, nc] ==1 or dry_cell_buffer[nb+1, nc] ==1:
                out[nfound]= n
                nfound += 1
        return out[:nfound]

    def return_to_last_position(self,nb,t2,sel):
        si = self.shared_info
        p_prop = si.classes['particle_properties']
        # return to last good position
        p_prop['x'].copy_prop('x_last_good', sel)

        # ensure bc cords and cell are right for any returned particles
        # particle
        si.classes['field_group_manager'].setup_interp_time_step(nb, t2, p_prop['x'].dataInBufferPtr(), sel)


    def screen_output(self,n_steps, nt0, nb0,  nb,  ns, t1, t0_step,computation_started):

        si= self.shared_info
        fraction_done= abs((t1 -si.model_start_time) / si.model_duration)
        s = '%02.0f%%:' % (100* fraction_done)
        s += '%06.0f:' % n_steps + 'h%06.0f:' % (nt0+nb-nb0) + 's%02.0f:' % ns + 'b%03.0f:' % nb

        t = abs( t1-si.model_start_time)
        s += 'Day ' +  ('-' if si.backtracking else '+')
        s += time_util.day_hms(t)
        s += ' ' + time_util.seconds_to_pretty_str(t1) + ':'
        s +=    si.classes['particle_group_manager'].screen_info()
        s += ' Finishes: ' + (datetime.now() +(1.-fraction_done)*(datetime.now()-computation_started)).strftime('%y-%m-%d %H:%M')
        timePerStep = perf_counter() - t0_step
        s +=  ' Step-%4.0f ms' % (timePerStep * 1000.)

        si.case_log.write_msg(s)

    def close(self):
        a=1



