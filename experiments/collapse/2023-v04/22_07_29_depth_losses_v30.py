#%% 

import numpy as np
from oceantracker import main


#-----------------------------------------------    
run_name = '22_11_01_depth_losses_v30'
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

# v18
# adding particle radius, density to calculate
# a terminal velocity based on stokes law
# with density and radius being changed by particles
# colliding with SPM particles and sticking to them

# v20
# merged dev041 to fix vertical diffusion bug

# v21
# using hzg3 containing montly averaged spm fields
# 1month run with 60s output

# v22
# testing new multi case run
# SA for initial particle size and stickiness

# v23
# small fixing and debugging run

# v28
# based on v23 after identifying a bug in regridding

# v29
# added the grid workaround again. has not been fixed as assumed
# also removed emtpy case from case list 
# removed unnecessary fields from the output

# v30
# adding "cause_of_death"


input_dir = "/scratch/local1/hzg3/"
output_dir = "/scratch/local1/output/"
# input_dir = "/work/uh0296/u301513/hzg_data/"
# output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:

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
    # 'numba_caching': False,
    # "compact_mode": True,
    "processors": 1,
    "replicates": 1,

    "regrid_z_to_uniform_sigma_levels": False,

    # "run_params":tm
    "write_tracks": True,
    "max_run_duration": 3600*24*14,
    "max_particles": int(1e6),
    "open_boundary_type": 0,
    "block_dry_cells": True,
    "time_step": 60,
    "screen_output_time_interval": 24*3600,
    "write_dry_cell_flag": True,

	"reader": {
		"input_dir": input_dir,
		"file_mask": "schout_*.nc",
        "load_fields": [
            "turbidity",
            "spm_sum_of_all_classes",
            "spm_very_fine_silt",
            "spm_fine_silt",
            "spm_medium_silt",
            "spm_coarse_silt",
            "spm_very_fine_sand",
            "salinity"
        ],
		"field_variable_map": {
			"water_depth": "depth",
			"salinity": "salt",
            "turbidity": "spm_sum_of_all_classes",
            "spm_sum_of_all_classes": "spm_sum_of_all_classes",
            "spm_very_fine_silt": "spm_very_fine_silt",
            "spm_fine_silt": "spm_fine_silt",
            "spm_medium_silt": "spm_medium_silt",
            "spm_coarse_silt": "spm_coarse_silt",
            "spm_very_fine_sand": "spm_very_fine_sand",
            "A_Z_profile": "diffusivity"
		},
	},

    "solver": {
        "RK_order": 2,
    },

    "dispersion": {
        'class_name': 'oceantracker.dispersion.random_walk_varyingAz.RandomWalkVaryingAZ',
        "A_H": 0.1,
        # "A_V": 0.001
    },

    "particle_properties": {  

        "dryout": {
            "class_name": "oceantracker.particle_properties.stranded_dryout.StrandedDryout",
            "max_time_stranded": 3600*24*1
        },
        "illumination": {
            'class_name': 'oceantracker.particle_properties.illumination.AverageIllumination',
            'name_of_turbidity_field': 'spm_sum_of_all_classes',
            'name_of_irradiance_field': 'irradiance',
            'c': 5,
            'time_to_average': 24*3600
        },
        "cause_of_death": {
            "class_name": "oceantracker.particle_properties.cause_of_death.CauseOfDeath",
        },
        # "density": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.Density",
        #     "initial_value": 1000
        # },
        # "radius": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.Radius",
        #     # "initial_value": 0.0005
        #     "initial_value": 0.5e-4
        # },
        # "buoyancy": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.Buoyancy",
        #     "gravity": 9.81,
        #     "mu": 1e-6  # kinematic viscosity of the water
        # },
        # collision with particle classes:
        # "collision_very_fine_silt": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
        #     "stickyness": stickiness,
        #     "spm_field": "spm_very_fine_silt",
        #     "spm_radius": 6e-6/2,
        #     "spm_density": 2650.
        # },
        # "collision_fine_silt": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
        #     "stickyness": stickiness,
        #     "spm_field": "spm_fine_silt",
        #     "spm_radius": 12e-6/2,
        #     "spm_density": 2650.
        # },
        # "collision_medium_silt": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
        #     "stickyness": stickiness,
        #     "spm_field": "spm_medium_silt",
        #     "spm_radius": 24e-6/2,
        #     "spm_density": 2650.
        # },
        # "collision_coarse_silt": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
        #     "stickyness": stickiness,
        #     "spm_field": "spm_coarse_silt",
        #     "spm_radius": 47e-6/2,
        #     "spm_density": 2650.
        # },
        # "collision_very_fine_sand": {
        #     "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
        #     "stickyness": stickiness,
        #     "spm_field": "spm_very_fine_sand",
        #     "spm_radius": 94e-6/2,
        #     "spm_density": 2650.
        # },
    },

    "fields": {
        "irradiance": {
            'class_name': 'oceantracker.fields.irradiance.Irradiance',
            'name_of_field': 'spm_sum_of_all_classes',
            'latitude': 53.5,
            'longitude': 9.9,
            'timezone': 'UTC',
            'albedo': 0.1
        },
    }, 

    "resuspension": {
        "critical_friction_velocity": 0.009
    },
    
    # "velocity_modifiers": {
    #     "buyoancy_based_terminal_velocity": {
    #         "class_name": "oceantracker.velocity_modifiers.stokes_based_buoyancy.StokesBasedBuoyancy"
    #     }
    # },

    "trajectory_modifiers": {
        "salinity_induced_mortality": {
            "class_name": "oceantracker.trajectory_modifiers.cull_particles.ParticleCullConcentration",
            "cull_interval": 3600*24,
            "cull_status_greater_than": "dead",
            "threshold_to_cull": 20,
            "probability_of_culling": 0.005,
            "concentration_field": "salinity"
        },
        # "light_starvation_induced_mortality": {
        #     "class_name": "oceantracker.trajectory_modifiers.cull_particles.IlluminationBasedLightLimitation",
        #     "cull_interval": 60,
        #     "cull_status_greater_than": "dead",
        #     "required_illumination": 30, #W/m2
        #     "probability_of_culling": 3.56e-05,
        # }
    },

    "release_groups": {
        "poly1": {
            "class_name": "oceantracker.release_groups.polygon_release.PolygonRelease",
            "points": release_polygon[0]['points'],
            "pulse_size": 1,
            "release_interval": 60*30
        }
    },
    "tracks_writer": {
        "update_interval": int(3600*1),
        "turn_off_write_particle_properties_list": [
                    "water_velocity",
                    "particle_velocity",
                    "velocity_modifier"

                    # "dry_cell_index",
                    # "particles_written_per_time_step",
                    # "particle_ID",
                    # "write_step_index",
                    # "time_step_range",
                    # "time",
                    # "num_part_released_so_far",
                    # "x",
                    # "status",
                    # "age",
                    # "ID",
                    # "IDrelease_group",
                    # "user_release_groupID",
                    # "IDpulse",
                    # "time_released",
                    # "x_last_good",
                    "turbidity",
                    "spm_very_fine_sand",
                    "spm_very_fine_silt",
                    # "tide",
                    "water_depth",
                    "spm_medium_silt",
                    # "salinity",
                    "spm_coarse_silt",
                    "spm_fine_silt",
                    "spm_sum_of_all_classes",
                    "irradiance",
                    "density",
                    "buoyancy",
                    "radius",
                    "collision_very_fine_silt",
                    "collision_fine_silt",
                    "collision_medium_silt",
                    "collision_coarse_silt",
                    "collision_very_fine_sand",
                    # "dryout",
                    # "illumination",
                    "release_points",
                    "number_of_release_points",
                    "is_polygon_release",                    
                ],
        # "time_steps_per_per_file": 24,
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

initial_size_list = np.linspace(1e-6,1e-4,3)
stickiness_list = np.linspace(1e-4,1e-2,3)

#prepand 0 stickiness and initial size of 1e-5
initial_size_list = np.insert(initial_size_list, 0, 1e-5)
stickiness_list = np.insert(stickiness_list, 0, 0)

initial_size_list = [1e-5]
stickiness_list = [0]

case_list = []

# reference case - no collision, neutral buoyancy

for initial_size in initial_size_list:
    for stickiness in stickiness_list:
        case_list.append({
            "particle_properties": {
                "density": {
                    "class_name": "oceantracker.particle_properties.buoyancy.Density",
                    "initial_value": 1000
                },
                "buoyancy": {
                    "class_name": "oceantracker.particle_properties.buoyancy.Buoyancy",
                    "gravity": 9.81,
                    "mu": 1e-6  # kinematic viscosity of the water
                },
                "radius": {
                    "class_name": "oceantracker.particle_properties.buoyancy.Radius",
                    "initial_value": initial_size
                },
                "collision_very_fine_silt": {
                    "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
                    "stickyness": stickiness,
                    "spm_field": "spm_very_fine_silt",
                    "spm_radius": 6e-6/2,
                    "spm_density": 2650.
                },
                "collision_fine_silt": {
                    "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
                    "stickyness": stickiness,
                    "spm_field": "spm_fine_silt",
                    "spm_radius": 12e-6/2,
                    "spm_density": 2650.
                },
                "collision_medium_silt": {
                    "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
                    "stickyness": stickiness,
                    "spm_field": "spm_medium_silt",
                    "spm_radius": 24e-6/2,
                    "spm_density": 2650.
                },
                "collision_coarse_silt": {
                    "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
                    "stickyness": stickiness,
                    "spm_field": "spm_coarse_silt",
                    "spm_radius": 47e-6/2,
                    "spm_density": 2650.
                },
                "collision_very_fine_sand": {
                    "class_name": "oceantracker.particle_properties.buoyancy.ParticleCollision",
                    "stickyness": stickiness,
                    "spm_field": "spm_very_fine_sand",
                    "spm_radius": 94e-6/2,
                    "spm_density": 2650.
                },

            },
            "velocity_modifiers": {
                "buyoancy_based_terminal_velocity": {
                "class_name": "oceantracker.velocity_modifiers.stokes_based_buoyancy.StokesBasedBuoyancy"
                }
            }
        })




#%%
case_info_files= main.run_parallel(params, case_list)

        
# quick fix to yet another bug in ross's HEAD

import os

def modify_json_files(folder_path):
    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.json'):
            filepath = os.path.join(folder_path, filename)

            # Read the content of the file
            with open(filepath, 'r') as file:
                filedata = file.read()

            # Replace the target string
            target_string = '"solver": null\n    },'
            replacement_string = '"solver": null,\n        "grid": "22_11_01_depth_losses_v29_C000_grid.nc",\n        "grid_outline": "22_11_01_depth_losses_v29_C000_grid_outline.json"\n    },'
            if target_string in filedata:
                filedata = filedata.replace(target_string, replacement_string)

                # Write the modified content back to the file
                with open(filepath, 'w') as file:
                    file.write(filedata)
            else:
                print(f"No target string found in {filename}")

# Replace 'your_folder_path' with the path to the folder containing your JSON files
modify_json_files(os.path.join(params['root_output_dir'],params['output_file_base']))


