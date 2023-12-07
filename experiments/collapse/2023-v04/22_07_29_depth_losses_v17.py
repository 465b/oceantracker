#%% 
import os
import numpy as np

from oceantracker import main


#-----------------------------------------------    
run_name = '22_11_01_depth_losses_v17'
#-----------------------------------------------


# v10
# upgraded from ot v03 to v04

# v11
# reduced release pulse interval as buffer is overflowing

# v12
# full year instead of 3 months and buffer increase to find equilibrium pop.

# v13
# small demo version (1 month) with hourly output

# v14
# trying to add turbidity to the model      
# reducing particles and model time to speed up testing

# v15
# adding average illumination to the model

# v16
# addingTurbidityInducedPrecipitationLikeSinking
# pleanty of debugging and getting everything to work

# v17
# adding illumination based culling

input_dir = "/scratch/local1/hzg2/"
output_dir = "/scratch/local1/output/"
# input_dir = "/work/uh0296/u301513/hzg_data/"
# output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
max_time = 3600*24*30*3
max_particle = 1e6

pulse_size = 1
pulse_interval = 3600

threshold_to_cull = 10
fraction_to_cull = 1


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
    "root_output_dir": output_dir,
    "output_file_base": run_name,
    "debug": False,
    # "compact_mode": True,
    # "processors": processors,
    # "replicates": replicates,

    # "run_params":
    "write_tracks": True,
    "max_run_duration": max_time,
    "max_particles": int(max_particle),
    "open_boundary_type": 0,
    "block_dry_cells": True,
    "time_step": 60,
    "screen_output_time_interval": 24*3600,

	"reader": {
		"input_dir": input_dir,
		"file_mask": "schout_*.nc",
		"field_variables": {
			"water_depth": "depth",
			"salinity": "salt",
            "turbidity": "turbidity",
            # "A_Z": "diffusivity"
		},
	},

    "solver": {
        "RK_order": 2,
    },

    "dispersion": {
        # 'class_name': 'oceantracker.dispersion.random_walk_varyingAz.RandomWalkVaryingAZ',
        "A_H": 0.1,
        "A_V": 0.001
    },

    "particle_properties": {  
        "total_water_depth": {
            "class_name": "oceantracker.particle_properties.total_water_depth.TotalWaterDepth"
        },
        "dryout": {
            "class_name": "oceantracker.particle_properties.stranded_dryout.StrandedDryout",
            "max_time_stranded": 3600*24*1
        },
        "illumination": {
            'class_name': 'oceantracker.particle_properties.illumination.AverageIllumination',
            'name_of_turbidity_field': 'turbidity',
            'name_of_irradiance_field': 'irradiance',
            'c': 5,
            'time_to_average': 24*3600
        },

    },
    "fields": {
        "irradiance": {
            'class_name': 'oceantracker.fields.irradiance.Irradiance',
            # 'name_of_field': 'turbidity',
            'latitude': 53.5,
            'longitude': 9.9,
            'timezone': 'UTC',
            'albedo': 0.1
        }
    #     "VerticalGradient": {
    #         'class_name': 'oceantracker.fields.field_vertical_gradient.VerticalGradient',
    #         'name_of_field': 'A_Z',
    #         'name': 'A_Z_vertical_gradient'
    #     }
    }, 

    "resuspension": {
        "critical_friction_velocity": 0.000
    },
    
    "velocity_modifiers": {
        # "terminal_velocity": {
        #     "class_name": "oceantracker.velocity_modifiers.terminal_velocity.TurbidityInducedPrecipitationLikeSinking",
        #     "turbidity_field": "turbidity",
        #     "diameter_of_particles": 0.001,
        #     "diameter_of_dust": 0.001,
        #     "density_dust": 1000,
        #     "max_sinking_velocity": -0.005
        # }
    },

    "trajectory_modifiers": {
        "salinity_induced_mortality": {
            "class_name": "oceantracker.trajectory_modifiers.cull_particles.ParticleCullConcentration",
            "cull_interval": 3600*24,
            "cull_status_greater_than": "dead",
            "threshold_to_cull": threshold_to_cull,
            "probability_of_culling": fraction_to_cull,
            "concentration_field": "salinity"
        },
        "light_starvation_induced_mortality": {
            "class_name": "oceantracker.trajectory_modifiers.cull_particles.IlluminationBasedLightLimitation",
            "cull_interval": 60,
            "cull_status_greater_than": "dead",
            "required_illumination": 30, #W/m2
            "probability_of_culling": 3.56e-05,
        }
    },

    "release_groups": {
        "poly1": {
            "class_name": "oceantracker.release_groups.polygon_release.PolygonRelease",
            "points": release_polygon[0]['points'],
            "pulse_size": int(pulse_size),
            "release_interval": pulse_interval
        }
    },

    "tracks_writer": {
        "update_interval": int(24*3600),
        "write_dry_cell_index": True
    },

    # 'particle_concentrations': {
    #         "top_layer": {
    #             "class_name": 'oceantracker.particle_concentrations.particle_concentrations.ParticleConcentrationsDepthLayer',
    #             "update_interval": 24*3600
    #         },
    #         "full_coloumn": {
    #             "class_name": 'oceantracker.particle_concentrations.particle_concentrations.ParticleConcentrations2D',
    #             "update_interval": 24*3600
    #         }
    # },
}

#%%
runInfo = main.run(params)
