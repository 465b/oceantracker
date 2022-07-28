from os import path, mkdir
import argparse
import glob
import numpy as np

from oceantracker.util import rUtil
import oceantracker.user_post_processing.loadOutputFiles as loadOutputFiles
from oceantracker.oceantrackersim import OceanTrackerRunner
import  oceantracker.user_post_processing.particlePlot as particlePlot
from copy import deepcopy
import  oceantracker.oceantrackersim

from matplotlib import animation

import yaml
#todo why does bouyant particle not follow free surface?


poly_points=[[1597682.1237, 5489972.7479],
                        [1598604.1667, 5490275.5488],
                        [1598886.4247, 5489464.0424],
                        [1597917.3387, 5489000],
                        [1597300, 5489000], [1597682.1237, 5489972.7479]]
poly_points_large=[[1597682.1237, 5489972.7479],
                        [1598604.1667, 5490275.5488],
                        [1598886.4247, 5489464.0424],
                        [1597917.3387, 5487000],
                        [1597300, 5487000],
                       [1597682.1237, 5489972.7479]]
demo_base_params=\
{'shared_params' :{'input_dir': 'demohindcast','output_file_base' :'demo01_plot_tracks',
                   'root_output_dir': 'output', 'debug': True},
 'reader':  {   'file_mask': 'demoHindcast2D*.nc',
                'water_velocity_map': {'u' : 'east_vel', 'v' : 'north_vel'},
                'dimension_map': {'time': 'time'},
                'grid_map': {   'time': 'time_sec', 'x': 'east','y': 'north',
                                'triangles': 'tri', 'dry_cells': 'dry_cell_flag'},
                'field_map' :{'water_depth' : 'depth'},
                'time_buffer_size': 24, 'isodate_of_hindcast_time_zero': '1970-01-01'} ,
 'base_case_params' : {
    'run_params' : {'dispersion_on' : True},
    'dispersion': {'A_H': 0.1},
    'backtracking': False,
    'solver': {'n_sub_steps': 2},
        'particle_release_groups': [ {'points': [[1594500, 5483000]], 'pulse_size': 200, 'release_interval': 0}],
        'user_particle_properties': [
                               #     {'name': 'Oxygen', 'class_name': 'oceantracker.user_particle_properties.ageDecay.AgeDecay', 'decay_time_scale': 1. * 3600 * 24,'initial_value' : 20.},
                                        {'class_name': 'oceantracker.user_particle_properties.distanceTravelled.DistanceTravelled'}
                                        ]},
                        }

p1= deepcopy(demo_base_params)
p1['base_case_params']['run_params']['backtracking']= True


params = []

# schsim 3D
s2 = deepcopy(schsim_base_params)
s2['shared_params'].update({'output_file_base' :'demo55_SCHISM_3D_fall_velocity' })
s2['reader'].update({ 'depth_average': False})

s2['base_case_params']['particle_release_groups']=[{'points': [[1594500, 5487000, -1]], 'pulse_size': 10, 'release_interval': 3600},
                                   {'class_name': 'oceantracker.particle_release_groups.polygonRelease.PolygonRelease',
                                    'points': poly_points,'z_min' : -2,'z_max' : -4,
                                    'pulse_size': 10, 'release_interval':  3600}
                                   ]
s2['base_case_params']['dispersion'].update( {'A_H':0.1, 'A_V':0.001})
s2['base_case_params']['user_velocity_modifiers']= [
       {'class_name' : 'oceantracker.particle_velocity.terminal_velocity.AddTerminalVelocity', 'mean': -0.002}
]

s2['base_case_params']['user_onfly_particle_statistics']=[
                  {   'class_name': 'oceantracker.user_onfly_particle_statistics.GriddedStatistics.GriddedStats2D_timeBased',
                      'calculation_interval': 3600, 'particle_property_list': ['water_depth'],
                      'grid_size': [120, 121]}]
params.append(s2)


def build_demos(root_folder, root_output_dir_base):

    print('building demos in ' + root_output_dir_base)


    JSONdir =path.join(root_folder, 'JSONinputFiles')
    if not path.isdir(JSONdir):
        mkdir(JSONdir)

    YAMLdir = path.join(root_folder, 'YAMLinputFiles')
    if not path.isdir(YAMLdir):
        mkdir(YAMLdir)

    input_folder = path.join(root_folder, 'demohindcast')
    print('demo hindcast in ', input_folder)
    for demo in params:
        demo['shared_params']['debug'] = True


        demo['shared_params']['input_dir'] = input_folder
        rUtil.write_JSON( path.join(JSONdir,demo['shared_params']['output_file_base'] + '.json'), demo)
        rUtil.write_YAML( path.join(YAMLdir, demo['shared_params']['output_file_base'] + '.yaml'), demo)


if __name__ == "__main__":

    root_folder = path.normpath(path.abspath(path.dirname(path.dirname(oceantracker.oceantrackersim.__file__))))
    root_folder = path.join(root_folder,'experiments')

    print('output dir ',root_folder, root_folder)

    build_demos(root_folder,root_folder)

    ot = OceanTrackerRunner()
    runInfo_file_name= ot.run(params)

