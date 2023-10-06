#%% 
import os
import numpy as np

from oceantracker import main
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import plot_statistics

#%%
#-----------------------------------------------    
run_name = '23_08_17_bath_jump_traj_v01'
#-----------------------------------------------

# v01
# - added sinking velocity

"""
The intention of this experiment is to create a dataset to observe particle
trajectories around the bathymetric jump in the harbor
"""

input_dir = "/scratch/local1/hzg2/"
output_dir = "/scratch/local1/output/"
# input_dir = "/work/uh0296/u301513/hzg_data/"
# output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60*6
processors = 10

max_time = 3600*24*28*1

max_particle = 1e4

pulse_size = 1
pulse_interval = 10*60 # every minute

threshold_to_cull = 10
fraction_to_cull = 1

sa_resolution = 0
replicates = 1

output_step_multiplier = 0.001 # hours between track recorded

release_polygon = [
    {
        'points': np.array([
					[593043,5919010],
                    [592988,5919053],
					[592924,5918944],
					[592981,5918899]
        ])
    }
]

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
			"water_salinity": "salt",
            "A_Z": "diffusivity"
		},
        'field_variables_to_depth_average': ['water_velocity'],
		"time_buffer_size": 24,
	},
	"base_case_params": {
		"run_params": {
			"write_tracks": True,
            "duration": max_time,
			"particle_buffer_size": int(max_particle),
            "open_boundary_type": 0,
            "block_dry_cells": True
		},
        "solver": {
            "RK_order": 2,
            "n_sub_steps": n_sub_steps,
            "screen_output_step_count": n_sub_steps
        },
        "tracks_writer": {
            "output_step_count": int(24*n_sub_steps*output_step_multiplier)
        },
        "dispersion": {
                'class_name': 'oceantracker.dispersion.random_walk_varyingAz.RandomWalkVaryingAZ',
                'A_H': 0.1,
            },
        "particle_properties": [
            {
                'class_name': 'oceantracker.particle_properties.total_water_depth.TotalWaterDepth'
            }
        ],
        "velocity_modifiers": [
            {
                "class_name": "oceantracker.velocity_modifiers.terminal_velocity.TerminalVelocity",
                "mean": -0.001,
            } 
        ],
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.resuspension.BasicResuspension",
                "critical_friction_velocity": 0.0
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
                {
                    'class_name': 'oceantracker.fields.friction_velocity.FrictionVelocity'
                },
                {
                    'class_name': 'oceantracker.fields.field_vertical_gradient.VerticalGradient',
                    'name_of_field': 'A_Z',
                    'name': 'A_Z_vertical_gradient'
                }                
        ], 
		"particle_release_groups": [
			{
				"class_name": "oceantracker.particle_release_groups.polygon_release.PolygonRelease",
				"points": release_polygon[0]['points'],
				"pulse_size": int(pulse_size),
				"release_interval": pulse_interval
			}
        ],
        "particle_statistics": [
        ],
        'particle_concentrations': [
        ],
    }
}

#%%
runInfo = main.run(params)

#%%
path_to_dir = os.path.join(params['shared_params']['root_output_dir'],params['shared_params']['output_file_base'])
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

