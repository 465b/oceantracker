from oceantracker.oceanTrackerMain import run_oceantracker

params={
    'shared_params' :{
                    "input_dir": "/scratch/local1/hzg2/",
		            "output_file_base": "22_04_04_testing_load_depth_averaged",
		            "root_output_dir": "/scratch/local1/output/",
                    "debug": True,
		            "write_to_screen": True,
                    },
	"reader": {
		"class_name": "oceantracker.readers.SchismNCDFreader.SCHSIMreaderNCDF",
		"file_mask": "schout_*.nc",
		"depth_average": False,
        "load_depth_aver_vel": True
	},
	"base_case_params": {
		"run_params": {
            "duration": 3600*24*1,
			"particle_buffer_size": 1100,
		},
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
				"pulse_size": 1000,
				"z_min": -1,
				"z_max": 1,
				"release_interval": 0
			}
        ]
    }
}

runInfo_file_name = run_oceantracker(params)
