#%%
from concurrent.futures import process
from oceantracker import main
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 

import os
import numpy as np

#%%
#-----------------------------------------------
run_name = '23_01_12_mini_test_set_w_rep'
#-----------------------------------------------

input_dir = "/scratch/local1/hzg2/"
output_dir = "/scratch/local1/output/"

# tweeked parameters:
n_sub_steps = 2
processors = 10

max_time = 3600*24*1

max_particle = 100
max_splitting_ratio = 0.0001
threshold_to_cull = 10
fraction_to_cull = 0.01

sa_resolution = 3
replicates = 2

output_step_multiplier = 1 # hours between track recorded

sinking_parameters = np.linspace(-1e-2,+1e-2,sa_resolution)
# pop'ing the zero verticle velocity parameters because all =0 cases are the same
# as no verticle migration and just waste comp time
sinking_parameters = sinking_parameters[sinking_parameters != 0]
splitting_ratios = np.linspace(0,max_splitting_ratio,sa_resolution)


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
                "A_H": 1.0,
                "A_V": 0.001
            },
        "particle_properties": [
                {'class_name': 'oceantracker.particle_properties.stranded_dryout.stranded_dryout',
                 'max_time_stranded': 3600*24*7},
                {'class_name': 'oceantracker.particle_properties.light_limitation.light_limitation',
                 'max_time_wo_light': 3600*24*14},
                {'class_name' : 'oceantracker.particle_properties.distance_travelled.DistanceTravelled'}
        ], 
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.resuspension.BasicResuspension",
                "critical_friction_velocity": 0.005
            },
            {
                "class_name": "oceantracker.trajectory_modifiers.cull_particles.ParticleCullConcentration",
                "cull_interval": 60,
                "cull_status_greater_than": 'frozen',
                "threshold_to_cull": threshold_to_cull,
                "probability_of_culling": fraction_to_cull,
                "concentration_field": "water_salinity"
            }
        ],
        "fields": [
                {'class_name': 'oceantracker.fields.friction_velocity.FrictionVelocity'},
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
				"z_range": [-5,5],
				"release_interval": 0
			}
        ]
    }
}


#%% run the model
runInfo = main.run(params)
runInfo_file_name = runInfo[0]
