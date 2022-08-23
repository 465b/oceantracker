#%%
from concurrent.futures import process
from oceantracker import main
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 

import os
import numpy as np

#%%
#-----------------------------------------------
run_name = '22_08_19_hpc_testing_stranding_no_resus_v00'
#-----------------------------------------------

input_dir = "/work/uh0296/u301513/hzg_data/"
output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60
processors = 10

max_time = 3600*24*20

max_particle = 2000
max_splitting_ratio = 0.0001
threshold_to_cull = 10
fraction_to_cull = 0.01

sa_resolution = 11
replicates = 1

output_step_multiplier = 8 # hours between track recorded


params={
    'shared_params' :{
	    "output_file_base": run_name,
        "root_output_dir": output_dir,
        "compact_mode": True,
        "processors": processors,
        "replicates": replicates,
    },
	"reader": {
		"class_name": "oceantracker.reader.schism_reader.SCHSIMreaderNCDF",
        "input_dir": input_dir,
		"file_mask": "schout_*.nc",
		"field_variables": {
			"water_depth": "depth",
			"water_salinity": "salt"
		},
		"time_buffer_size": 24,
        "field_variables_to_depth_average":['water_velocity']
	},
	"base_case_params": {
		"run_params": {
			"write_tracks": True,
            "duration": max_time,
			"particle_buffer_size": max_particle,
            "open_boundary_type": 0,
            "block_dry_cells": True
		},
        "solver": {
            "RK_order": 2,
            "n_sub_steps": n_sub_steps,
            "screen_output_step_count": int(n_sub_steps*output_step_multiplier)
        },
        "tracks_writer": {
            "output_step_count": int(n_sub_steps*output_step_multiplier)
        },
        "dispersion": {
                "A_H": 0.1,
                "A_V": 0.01
            },
        "particle_properties": [
                {'class_name': 'oceantracker.particle_properties.friction_velocity.FrictionVelocity'}
        ], 
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.resuspension.BasicResuspension",
                "critical_friction_velocity": 0.0
            }
        ],
		"particle_release_groups": [
			{
				"class_name": "oceantracker.particle_release_groups.polygon_release.PolygonRelease",
				"points": [
					[
						593043,
						5919010
					],
					[
						592988,
						5919053
					],
					[
						592924,
						5918944
					],
					[
						592981,
						5918899
					]
				],
				"pulse_size": int(max_particle/2),
				"z_min": -1,
				"z_max": 1,
				"release_interval": 0
			}
        ]
    }
}

#%% run the model
runInfo = main.run(params)
runInfo_file_name = runInfo[0]
print(runInfo_file_name)

#%% load and draw statistical overview
path_to_dir = os.path.join(params['shared_params']['root_output_dir'],params['shared_params']['output_file_base'])
print(f'path_to_dir: {path_to_dir}')
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

# draw overview of all cases split up into repetitions and migration
# stats_plot.plot_sa_total_polycount(cases,output_path=path_to_dir)

#%% draw animations
stats_plot.animate_cases(cases)
