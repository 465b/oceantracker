from oceantracker.oceanTrackerMain import run_oceantracker
from oceantracker.user_post_processing import loadOutputFiles
from oceantracker.user_post_processing import statsPlot
import numpy as np


statistical_polygon_list = [
    {
        "__comment": "geesthacht",
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
        "__comment": "neuengamme",
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
        "__comment": "kirchwerder",
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
        "__comment": "wilhelmsburg",
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
        "__comment": "harbor",
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
        "__comment": "schulau",
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
        "__comment": "stade",
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
        "__comment": "glückstadt",
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
        "__comment": "freiburg",
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
        "__comment": "brunsbüttel",
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
        "__comment": "ottendorf",
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
        "__comment": "cuxhafen",
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
        "__comment": "north_sea",
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
                    "input_dir": "/scratch/local1/hzg2/",
		            "output_file_base": "22_04_22_light_dryness_frictionresusp_tests_v04",
		            "root_output_dir": "/scratch/local1/output/",
		            "compact_mode": True,
		            "processors": 10,
		            "replicates": 1,
                    "multiprocessing_start_method_spawn": False,
		            "write_to_screen": True,
                    "debug": True
                    },
	"reader": {
		"class_name": "oceantracker.readers.SchismNCDFreader.SCHSIMreaderNCDF",
		"file_mask": "schout_*.nc",
		"field_variables": {
			"water_depth": "depth",
			"salinity": "salt"
		},
		"time_buffer_size": 24,
        "field_variables_to_depth_average":['water_velocity']
	},
	"base_case_params": {
		"run_params": {
			"write_tracks": True,
            "duration": 3600*24*365,
			"particle_buffer_size": 11000,
		},
		"solver": {
			"n_sub_steps": 360,
            "screen_output_step_count": 360
		},
        "tracks_writer": {
            "output_step_count": 3600
        },
        "user_particle_properties": [
                {'class_name': 'oceantracker.user_particle_properties.stranded_dryout.stranded_dryout'},
                {'class_name': 'oceantracker.user_particle_properties.distanceTravelled.DistanceTravelled'},
                {'class_name': 'oceantracker.user_particle_properties.light_limitation.light_limitation'}
        ],
		"particle_release_groups": [
			{
				"class_name": "oceantracker.particle_release_groups.polygonRelease.PolygonRelease",
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
				"pulse_size": 10000,
				"z_min": -1,
				"z_max": 1,
				"release_interval": 0
			}
        ],
        "user_onfly_particle_statistics": [
            {
                "class_name": "oceantracker.user_onfly_particle_statistics.PolygonStatistics.PolygonStats2D_timeBased",
                "calculation_interval": 100,
                "count_status_equal_to": 10,
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            },
            {
                "class_name": "oceantracker.user_onfly_particle_statistics.PolygonStatistics.PolygonStats2D_timeBased",
                "calculation_interval": 100,
                "count_status_equal_to": 5,
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            },
            {
                "class_name": "oceantracker.user_onfly_particle_statistics.PolygonStatistics.PolygonStats2D_timeBased",
                "calculation_interval": 100,
                "count_status_equal_to": 4,
                "particle_property_list": [],
                "polygon_list": statistical_polygon_list
            }
        ]
    }
}

sinking_parameters = np.linspace(-1e-2,+1e-2,5)
splitting_ratios = np.linspace(0,0.00025,5)

cases = []
# Monotonic Migration
for sinking in sinking_parameters:
    for splitting in splitting_ratios:
        cases.append({
            "user_velocity_modifiers": [
                {
                    "class_name": "oceantracker.particle_velocity.terminal_velocity.AddTerminalVelocity",
                    "mean": sinking
                }
            ],
            "user_trajectory_modifiers": [
                {
                    "class_name": "oceantracker.user_trajectory_modifiers.resuspension.BasicResuspension",
                    "critical_friction_velocity": 0.01
                },
                {
                    "class_name": "oceantracker.user_trajectory_modifiers.splitParticles.ParticleSplit",
                    "splitting_interval": 100,
                    "split_status_greater_than": 4,
                    "probability_of_splitting": splitting
                },
                # mortality for particles reaching high salinity water
                {
                    "class_name": "oceantracker.user_trajectory_modifiers.cullParticles.ParticleCullConcentration",
                    "cull_interval": 100,
                    "status_to_cull": 0,
                    "threshold_to_cull": 1.0,
                    "fraction_to_cull": 0.1,
                    "concentration_field": "salinity"
                }
            ]
        })
        

params['case_list'] = cases

runInfo_file_name = run_oceantracker(params)
print('done')

runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(runInfo_file_name)
statsPlot.plot_sa_total_polycount(runCaseInfo,[5,5])
statsPlot.animate_cases(runCaseInfo)