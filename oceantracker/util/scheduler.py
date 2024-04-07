from inspect import signature
import numpy as np
from oceantracker.util import  time_util
class Scheduler(object):
    # set up event shedule based on times since 1/1/1970
    # rounds starts, times and intervals to model time steps,
    # uses times given, otherwise start and interval
    # all times in seconds
    def __init__(self,settings, run_info,hindcast_info,
                 start=None, end=None, duration=None,
                 interval = None, times=None,cancel_when_done=True):


        self.cancel_when_done = cancel_when_done
        md = run_info.model_direction
        dt = settings.time_step

        if times is None:
            # make from start time and interval
            times,interval = self._start_end_from_interval(settings, run_info, start, end, duration, interval)
        else:
            # use times given, but round
            n = (times - run_info.start_time)/dt
            times = run_info.start_time + np.round(n) * dt
            times = md* np.sort(md*times) # ensure they are in right order for backwards/forwards
            interval = None

        start_time_outside_run_times = times[0]*md < run_info.start_time*md or times[-1]*md > run_info.end_time*md

        # trim to fit inside the run
        sel = np.logical_and( times * md  >= run_info.start_time * md, times * md  <= run_info.end_time * md )
        self.scheduled_times = times[sel]

        # make a task flag for each time step of the model run
        self.task_flag = np.full_like(run_info.times,False, dtype=bool)
        nt_task = (np.abs(self.scheduled_times - run_info.start_time) / settings.time_step).astype(np.int32)
        self.task_flag[nt_task] = True

        # flag times steps scheduler is active, ie start to end
        self.active_flag = np.full_like(run_info.times, False, dtype=bool)
        self.active_flag[nt_task[0]:nt_task[-1]+1] = True

        # record info
        duration=abs(self.scheduled_times[-1] - self.scheduled_times[0]),
        self.info= dict(start_time=self.scheduled_times[0], interval=interval,
                        end_time=self.scheduled_times[-1],
                        duration = duration,
                        duration_str = time_util.seconds_to_pretty_duration_string(duration),
                        interval_str = time_util.seconds_to_pretty_duration_string(interval),
                        start_date=time_util.seconds_to_isostr(self.scheduled_times[0]),
                        end_date=time_util.seconds_to_isostr(self.scheduled_times[-1]),
                        number_scheduled_times = self.scheduled_times.size,
                        cancel_when_done=cancel_when_done,
                        start_time_outside_run_times =start_time_outside_run_times,
                        )
        i = self.info
        b = f'{12*" "} Scheduler{15*" "}Run:{14*" "}Hindcast\n'
        b += f'Start- {i["start_date"]} {hindcast_info["start_date"]}  {time_util.seconds_to_isostr(run_info.times[0])}  \n'
        b += f'Ends - {i["end_date"]} {hindcast_info["end_date"]}  {time_util.seconds_to_isostr(run_info.times[-1])}]\n'
        b += f'{10*" "}interval = {i["interval_str"]}, backtracking={settings.backtracking}'
        i['bounds_table']= b
        pass

    def _start_end_from_interval(self,settings,run_info, start,end, duration, interval):

        md = run_info.model_direction
        dt = settings.time_step
        if start is None:
            start = run_info.start_time
        else:
            # use given start rounded to time step
            n = (start - run_info.start_time) / dt  # number of model steps since the start
            start = run_info.start_time + round(n) *  dt

        if duration is not None:
            # use duration for end if given
            end = start + md * duration
        elif end is None:
            end = run_info.end_time

        if interval is None:
            interval=  dt
        elif interval > 0.:
            # round interval to time step, but not less than one per time step
            interval=  max(round(interval/dt)*dt, dt)
        else:
            interval = 0.

        duration = abs(start-end)
        if interval ==0:
            # only one event
            times =np.asarray([start])
        else:
            times = start + md * np.arange(0,duration+settings.time_step,interval)

        return  times, interval

    def do_task(self, n_time_step):
        # check if task flag is set
        do_it = self.task_flag[n_time_step]

        if self.cancel_when_done and do_it:
            # ensure task is not repeated by another operation at the same time step
            self.task_flag[n_time_step] = False
        return do_it

    def cancel_task(self, n_time_step):
        # check if task flag is set
         self.task_flag[n_time_step] = False

    def see_task_flag(self, n_time_step):
        # returns if task is happening  without any cancellation of task when done
        return  self.task_flag[n_time_step]

    def is_active(self, n_time_step):
        # check if task is between start and end from active_flag is set
        return self.active_flag[n_time_step]