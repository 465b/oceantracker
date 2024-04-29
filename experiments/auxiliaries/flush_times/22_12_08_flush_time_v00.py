#%% 
import os
import numpy as np

from oceantracker import main

from oceantracker.post_processing.read_output_files import load_output_files 

from oceantracker.post_processing.plotting import stats_plot 
from oceantracker.post_processing.plotting import plot_utilities
from oceantracker.post_processing.plotting import plot_transects
from oceantracker.post_processing.plotting import plot_tracks
from oceantracker.post_processing.plotting import plot_statistics

from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_path_in_vertical_section
from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_relative_height

from oceantracker.post_processing.transect_projection.transect import Transect

#%%
#-----------------------------------------------
run_name = '22_12_08_flush_time_v00'
#-----------------------------------------------w

#v07
# added distance traveled 
# testing new data post processing

#v09
# added dynamic verticle dispersion

#v11
# reduced growth rate

#v12
# fixed buggy stats 

#v14
# fixed broken diel migration

#v15
# changed init light starvation to 0

#v21
# based on v15 - increased salinity threshold to 20

#v22
# set critical friction velocity to 0.009
# increased max particle count to 20000

#v23 & v24 & v25
# increased growth rate bound

# based on retention v25
# particles always resuspend get removed after stranded for a day
# no splitting




input_dir = "/work/uh0296/u301513/hzg_data/"
output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60
processors = 256

max_time = 3600*24*365

max_particle = 200000
threshold_to_cull = 20
fraction_to_cull = 1


output_step_multiplier = 12 # hours between track recorded


params={
    'shared_params' :{
	    "output_file_base": run_name,
        "root_output_dir": output_dir,
        "compact_mode": True,
        "processors": processors,
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
		"time_buffer_size": 24,
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
                "class_name": "oceantracker.dispersion.random_walk_varyingAz.RandomWalkVaryingAZ",
                "A_H": 0.1
        },
        "particle_properties": [
                {
                'class_name': 'oceantracker.particle_properties.total_water_depth.TotalWaterDepth'
                },
                {'class_name': 'oceantracker.particle_properties.stranded_dryout.stranded_dryout',
                 'max_time_stranded': 3600*24*1},
                {'class_name' : 'oceantracker.particle_properties.distance_travelled.DistanceTravelled'}
        ], 
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.resuspension.BasicResuspension",
                "critical_friction_velocity": 0.000
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
                "class_name": "oceantracker.fields.friction_velocity.FrictionVelocity"
            },
            {
                'class_name': 'oceantracker.fields.field_vertical_gradient.VerticalGradient','name_of_field': 'A_Z','name':'A_Z_vertical_gradient'
            },
        ],
		"particle_release_groups": [
			{
				"class_name": "oceantracker.particle_release_groups.polygon_release.PolygonRelease",
				"points": [
                    [ 600900, 5911947],
                    [ 600900, 6012889],
                    [ 435445, 6012889],
                    [ 435445, 5911947]
                ],
				"pulse_size": int(max_particle/2),
                # release every month
				"release_interval": 3600*24*28,
                "maximum_age": 3600*24*28
            },
        ],
    }
}

#%% run the model
runInfo = main.run(params)