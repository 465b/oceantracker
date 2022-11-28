#%%
from concurrent.futures import process
from oceantracker import main
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 

import os
import numpy as np

#%%
#-----------------------------------------------
run_name = '22_07_27_retention_v06'
#-----------------------------------------------

input_dir = "/work/uh0296/u301513/hzg_data/"
output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60
processors = 256  

max_time = 3600*24*365

max_particle = 10000
max_splitting_ratio = 2e-5
threshold_to_cull = 10
fraction_to_cull = 0.01

sa_resolution = 11
replicates = 1

output_step_multiplier = 12 # hours between track recorded

sinking_parameters = np.linspace(-1e-2,+1e-2,sa_resolution)
# pop'ing the zero verticle velocity parameters because all =0 cases are the same
# as no verticle migration and just waste comp time
sinking_parameters[sinking_parameters != 0]
splitting_ratios = np.linspace(0,max_splitting_ratio,sa_resolution)

statistical_polygon_list = [
    {
        "user_polygon_name": "geesthacht",
        "points": [
            [
                585700.0,
                5921404.0
            ],
            [
                586412.0,
                5920095.0
            ],
            [
                594187.0,
                5917410.0
            ],
            [
                594346.0,
                5919695.0
            ],
            [
                587003.0,
                5922492.0
            ]
        ]
    },
    {
        "user_polygon_name": "neuengamme",
        "points": [
            [
                577905.0,
                5916461.0
            ],
            [
                582811.0,
                5915980.0
            ],
            [
                586412.0,
                5920095.0
            ],
            [
                585700.0,
                5921404.0
            ],
            [
                583464.0,
                5921191.0
            ],
            [
                578488.0,
                5917969.0
            ]
        ]
    },
    {
        "user_polygon_name": "kirchwerder",
        "points": [
            [
                577905.0,
                5916461.0
            ],
            [
                578488.0,
                5917969.0
            ],
            [
                574562.0,
                5927524.0
            ],
            [
                568765.0,
                5921992.0
            ]
        ]
    },
    {
        "user_polygon_name": "wilhelmsburg",
        "points": [
            [
                568765.0,
                5921992.0
            ],
            [
                574562.0,
                5927524.0
            ],
            [
                574027.0,
                5931853.0
            ],
            [
                569497.0,
                5933189.0
            ],
            [
                568607.0,
                5931906.0
            ],
            [
                568330.0,
                5931212.0
            ],
            [
                569596.0,
                5930624.0
            ],
            [
                569893.0,
                5930303.0
            ],
            [
                570051.0,
                5929822.0
            ],
            [
                569774.0,
                5929582.0
            ],
            [
                567993.0,
                5930437.0
            ],
            [
                564729.0,
                5924023.0
            ]
        ]
    },
    {
        "user_polygon_name": "harbor",
        "points": [
            [
                567380.0,
                5936476.0
            ],
            [
                554995.0,
                5934712.0
            ],
            [
                554600.0,
                5932895.0
            ],
            [
                564729.0,
                5924023.0
            ],
            [
                567993.0,
                5930437.0
            ],
            [
                569774.0,
                5929582.0
            ],
            [
                570051.0,
                5929822.0
            ],
            [
                569893.0,
                5930303.0
            ],
            [
                569596.0,
                5930624.0
            ],
            [
                568330.0,
                5931212.0
            ],
            [
                568607.0,
                5931906.0
            ],
            [
                569497.0,
                5933189.0
            ]
        ]
    },
    {
        "user_polygon_name": "schulau",
        "points": [
            [
                554600.0,
                5932895.0
            ],
            [
                554995.0,
                5934712.0
            ],
            [
                539772.0,
                5943332.0
            ],
            [
                535644.0,
                5936297.0
            ],
            [
                554909.0,
                5929980.0
            ],
            [
                555597.0,
                5931806.0
            ]
        ]
    },
    {
        "user_polygon_name": "stade",
        "points": [
            [
                535644.0,
                5936297.0
            ],
            [
                539772.0,
                5943332.0
            ],
            [
                533518.0,
                5959826.0
            ],
            [
                522176.0,
                5952740.0
            ]
        ]
    },
    {
        "user_polygon_name": "glückstadt",
        "points": [
            [
                522176.0,
                5952740.0
            ],
            [
                533518.0,
                5959826.0
            ],
            [
                525679.0,
                5969373.0
            ],
            [
                517339.0,
                5963566.0
            ]
        ]
    },
    {
        "user_polygon_name": "freiburg",
        "points": [
            [
                517339.0,
                5963566.0
            ],
            [
                525679.0,
                5969373.0
            ],
            [
                513670.0,
                5975376.0
            ],
            [
                513586.0,
                5965633.0
            ]
        ]
    },
    {
        "user_polygon_name": "brunsbüttel",
        "points": [
            [
                513586.0,
                5965633.0
            ],
            [
                513670.0,
                5975376.0
            ],
            [
                498742.0,
                5980888.0
            ],
            [
                491403.0,
                5958744.0
            ]
        ]
    },
    {
        "user_polygon_name": "ottendorf",
        "points": [
            [
                491403.0,
                5958744.0
            ],
            [
                498742.0,
                5980888.0
            ],
            [
                490569.0,
                5986990.0
            ],
            [
                471722.0,
                5956775.0
            ]
        ]
    },
    {
        "user_polygon_name": "cuxhafen",
        "points": [
            [
                471722.0,
                5956775.0
            ],
            [
                490569.0,
                5986990.0
            ],
            [
                494155.0,
                5985612.0
            ],
            [
                501411.0,
                5986005.0
            ],
            [
                501577.0,
                5991222.0
            ],
            [
                496073.0,
                6002441.0
            ],
            [
                485732.0,
                5999095.0
            ],
            [
                475475.0,
                5999390.0
            ],
            [
                468219.0,
                5995158.0
            ],
            [
                463883.0,
                5989056.0
            ],
            [
                460130.0,
                5982561.0
            ],
            [
                457461.0,
                5975475.0
            ],
            [
                459630.0,
                5968782.0
            ],
            [
                462799.0,
                5961499.0
            ]
        ]
    },
    {
        "user_polygon_name": "north_sea",
        "points": [
            [
                471722.0,
                5956775.0
            ],
            [
                462799.0,
                5961499.0
            ],
            [
                459630.0,
                5968782.0
            ],
            [
                457461.0,
                5975475.0
            ],
            [
                463883.0,
                5989056.0
            ],
            [
                468219.0,
                5995158.0
            ],
            [
                475475.0,
                5999390.0
            ],
            [
                485732.0,
                5999095.0
            ],
            [
                496073.0,
                6002441.0
            ],
            [
                480479.0,
                6011414.0
            ],
            [
                451123.0,
                6010595.0
            ],
            [
                438114.0,
                5987654.0
            ],
            [
                438114.0,
                5957612.0
            ],
            [
                468803.0,
                5943956.0
            ]
        ]
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
                {'class_name': 'oceantracker.particle_properties.light_limitation.light_limitation',
                 'max_time_wo_light': 3600*24*14},
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
        ],
        "particle_statistics": [
            {
                "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
                "calculation_interval": 60,
                "count_status_equal_to": 'moving',
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            },
            {
                "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
                "calculation_interval": 60,
                "count_status_equal_to": 'stranded_by_tide',
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            },
            {
                "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
                "calculation_interval": 60,
                "count_status_equal_to": 'on_bottom',
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            }
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

#%% run the model
runInfo = main.run(params)
runInfo_file_name = runInfo[0]
print(runInfo_file_name)

#%% load and draw statistical overview
path_to_dir = os.path.join(params['shared_params']['root_output_dir'],params['shared_params']['output_file_base'])
print(f'path_to_dir: {path_to_dir}')
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

# draw overview of all cases split up into repetitions and migration
stats_plot.plot_sa_total_polycount(cases,output_path=path_to_dir)

#%% draw animations
stats_plot.animate_cases(cases)
