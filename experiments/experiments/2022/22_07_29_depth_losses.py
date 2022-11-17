#%% from oceantracker import main
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 
from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_path_in_vertical_section
from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_relative_height


from oceantracker.post_processing.transect_projection.transect import Transect
from oceantracker.post_processing.plotting import plot_transects

import os
import numpy as np

#%%

#-----------------------------------------------
run_name = '22_11_01_depth_losses_v01'
#-----------------------------------------------

input_dir = "/scratch/local1/hzg2/"
output_dir = "/scratch/local1/output/"
# input_dir = "/work/uh0296/u301513/hzg_data/"
# output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60
processors = 1 

max_time = 3600*24*7*4

max_particle = 100000

pulse_size = 5
pulse_interval = 60

max_splitting_ratio = 0.0001

threshold_to_cull = 10
fraction_to_cull = 0.01

sa_resolution = 0
replicates = 1

output_step_multiplier = 1 # hours between track recorded

sinking_parameters = np.linspace(-1e-2,+1e-2,sa_resolution)
# pop'ing the zero verticle velocity parameters because all =0 cases are the same
# as no verticle migration and just waste comp time
sinking_parameters[sinking_parameters != 0]
splitting_ratios = np.linspace(0,max_splitting_ratio,sa_resolution)

polygone = [
    {'points':
        np.array([
            [569979,5930284],
            [570239,5930483],
            [569147,5931534],
            [568979,5931251]
        ])
    },
    {'points':
        np.array([
            [568988,5931254],
            [569147,5931537],
            [567622,5932312],
            [567502,5932043]
        ])
    },
    {'points':
        np.array([
            [567584,5931996],
            [567704,5932265],
            [565924,5932656],
            [565900,5932289]
        ])
    },
    {'points':
        np.array([
            [565900,5932282],
            [565943,5932663],
            [564086,5933401],
            [563971,5933043]
        ])
    },
    {'points':
        np.array([
            [564000,5933064],
            [564096,5933417],
            [562586,5933178],
            [562643,5932757]
        ])
    },
    {'points':
        np.array([
            [562638,5932757],
            [562614,5933141],
            [561359,5933128],
            [561359,5932690]
        ])
    }
]

planes = [
    {'points':
        np.array([
            [570181,5930237],
            [568926,5931456]
        ])
    },
    {'points':
        np.array([
            [569123,5931325],
            [567444,5932218]
        ])
    },
    {'points':
        np.array([
            [567771,5932080],
            [565804,5932477]
        ])
    },
    {'points':
        np.array([
            [565977,5932423],
            [563961,5933259]
        ])
    },
    {'points':
        np.array([
            [564207,5933266],
            [562446,5932925]
        ])
    },
    {'points':
        np.array([
            [562704,5932919],
            [561174,5932947]
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
                {'class_name': 'oceantracker.particle_properties.depth_losses.depth_losses',
                 'initial_value': 0},
                {'class_name': 'oceantracker.particle_properties.friction_velocity.FrictionVelocity'}
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
        "derived_fields": [
           {"class_name": "oceantracker.derived_fields.total_water_depth.TotalWaterDepth"}
        ],
		"particle_release_groups": [
			{
				"class_name": "oceantracker.particle_release_groups.polygon_release.PolygonRelease",
				"points": polygone[0]['points'],
				"pulse_size": pulse_size,
				# "z_min": -1,
				# "z_max": 1,
				"release_interval": pulse_interval
			}
        ],
        "particle_statistics": [
            # {
            #     "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
            #     "calculation_interval": 60,
            #     "count_status_equal_to": 'moving',
            #     "particle_property_list": [],
            #     "polygon_list": statistical_polygon_list
            # },
            # {
            #     "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
            #     "calculation_interval": 60,
            #     "count_status_equal_to": 'stranded_by_tide',
            #     "particle_property_list": [],
            #     "polygon_list": statistical_polygon_list
            # },
            # {
            #     "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
            #     "calculation_interval": 60,
            #     "count_status_equal_to": 'on_bottom',
            #     "particle_property_list": [],
            #     "polygon_list": statistical_polygon_list
            # }
        ]
    }
}

cases = []

# drifting aka no migration
drifting_cases = []
sinking = 0
for splitting in splitting_ratios:
    drifting_cases.append({
        "velocity_modifiers": [
            {
                "class_name": "oceantracker.velocity_modifiers.terminal_velocity.AddTerminalVelocity",
                "mean": sinking
            }
        ],
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.split_particles.SplitParticles",
                "splitting_interval": 60,
                "split_status_greater_than": 'frozen',
                "probability_of_splitting": splitting
            }
        ]
    })

# Monotonic migration
monotonic_cases = []
for sinking in sinking_parameters:
    for splitting in splitting_ratios:
        monotonic_cases.append({
            "velocity_modifiers": [
                {
                    "class_name": "oceantracker.velocity_modifiers.terminal_velocity.AddTerminalVelocity",
                    "mean": sinking
                }
            ],
            "trajectory_modifiers": [
                {
                    "class_name": "oceantracker.trajectory_modifiers.split_particles.SplitParticles",
                    "splitting_interval": 60,
                    "split_status_greater_than": 'frozen',
                    "probability_of_splitting": splitting
                }
            ]
        })

# Diel Migration
diel_cases = []
for sinking in sinking_parameters:
    for splitting in splitting_ratios:
        diel_cases.append({
            "velocity_modifiers": [
                {
                    "class_name": "oceantracker.velocity_modifiers.terminal_velocity.AddDielVelocity",
                    "mean": sinking
                }
            ],
            "trajectory_modifiers": [
                {
                    "class_name": "oceantracker.trajectory_modifiers.split_particles.SplitParticles",
                    "splitting_interval": 60,
                    "split_status_greater_than": 'frozen',
                    "probability_of_splitting": splitting
                }
            ]
        })


cases = []
cases += drifting_cases
cases += monotonic_cases
cases += diel_cases
params['case_list'] = cases

print(f'Total amount of cases: {len(cases)}')

# runInfo = main.run(params)
# runInfo_file_name = runInfo[0]

#%%
path_to_dir = os.path.join(params['shared_params']['root_output_dir'],params['shared_params']['output_file_base'])
print(f'path_to_dir: {path_to_dir}')
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

#%%
stats_plot.animate_cases(cases)

#%%
transect = [{'polygon_vertices': poly['points'],'plane_vertices': plane['points']} for poly,plane in zip(polygone,planes)]
norderelbe_transect = Transect(transect)
norderelbe_transect.project_track_data(path_to_dir,var_list=['water_depth','age','tide','depth_losses'],n_case=0)

#%%
for ii in range(20):
    plot_path_in_vertical_section(norderelbe_transect.track_data,particleID=ii)

#%%
starving = norderelbe_transect.track_data['depth_losses']
starving = np.max(starving,axis=0)
for days_starving in range(8):
    print(f'Assuming phytoplankton can has an {days_starving} energy budget')
    print(f'Then {len(starving[starving > (3600*24*days_starving)])/len(starving)*100}% are surviving')
print('Maybe they are only surviving because they are not accounted for anymore (out of transect) and therefore not dying')

#%%
# for nt in range(0,100,1):
# # nt = 30
#     axis_lims = [560500,572500,5.923e6,5.935e6]
#     # plot_transects.plot_transect_map(norderelbe_transect,nt=nt,axis_lims=axis_lims,plot_file_name='transect_map.png')
#     plot_transects.plot_projected_horizontal_tracks(norderelbe_transect,nt=nt)
#     plot_transects.plot_projected_verticle_tracks(norderelbe_transect,nt=nt)




# %%
# threshold = -1
# for nt in range(0,120,10):
#     depth = norderelbe_transect.track_data['x'][nt,:,2]
#     depth = depth[~np.isnan(depth)]
#     print(f'len depth {len(depth)}')

#     n_above = len(depth[depth >= threshold])
#     n_below = len(depth[depth < threshold])
#     print(n_above,n_below,n_above/len(depth))