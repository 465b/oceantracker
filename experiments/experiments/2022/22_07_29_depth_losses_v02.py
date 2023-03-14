#%% 
import os
import numpy as np

from oceantracker import main

from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.read_output_files import load_output_files

from oceantracker.post_processing.transect_projection.transect import Transect

from oceantracker.post_processing.plotting import stats_plot 
from oceantracker.post_processing.plotting import plot_utilities
from oceantracker.post_processing.plotting import plot_transects
from oceantracker.post_processing.plotting import plot_tracks

from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_path_in_vertical_section
from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_relative_height


#%%

#-----------------------------------------------
run_name = '22_11_01_depth_losses_v02'
#-----------------------------------------------

input_dir = "/scratch/local1/hzg2/"
output_dir = "/scratch/local1/output/"
# input_dir = "/work/uh0296/u301513/hzg_data/"
# output_dir = "/work/uh0296/u301513/ot_output/"

# tweeked parameters:
n_sub_steps = 60
processors = 1 

max_time = 3600*24*14

max_particle = 1e6

pulse_size = 1000
pulse_interval = 60 # every minute

threshold_to_cull = 0.001
fraction_to_cull = 1

sa_resolution = 0
replicates = 1

output_step_multiplier = 1 # hours between track recorded

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


observational_polygone = [
    {
        'name': 'pre-splitting',
        'points': np.array([
           [572085,5922830],
           [572169,5922719],
           [572246,5922848],
           [572162,5922948] 
        ])
    },
    # norderelbe
    ## low-depth
    {
        'name': 'norder-bunthaus',
        'points': np.array([
            [570470,5925498],
            [570537,5925521],
            [570508,5925681],
            [570433,5925668]
        ])
    },
    {
        'name': 'norder-kaltehofe',
        'points': np.array([
            [569021,5931319],
            [569084,5931274],
            [569139,5931362],
            [569077,5931410]
        ])
    },
    ## medium-depth
    {
        'name': 'norder-elphi',
        'points': np.array([
            [565158,5932675],
            [565229,5932795],
            [565151,5932838],
            [565069,5932716]
        ])
    },
    ## high-depth
    {
        'name': 'norder-dockland',
        'points': np.array([
            [562387,5933002],
            [562390,5932860],
            [562513,5932872],
            [562507,5933012]
        ])
    },
    # suederelbe
    ## low-depth
    {
        'name': 'sueder-bunthaus',
        'points': np.array([
            [570405,5923818],
            [570454,5923894],
            [570397,5923946],
            [570345,5923871]
        ])
    },
    {
        'name': 'sueder-eurpabruecke',
        'points': np.array([
            [565993,5925470],
            [566099,5925433],
            [566083,5925340],
            [565977,5925374]
        ])
    },
    ## medium depth
    {
        'name': 'sueder-nynas',
        'points': np.array([
            [563601,5926690],
            [563675,5926768],
            [563630,5926837],
            [563552,5926761]
        ])
    },
    ## high depth
    {
        'name': 'sueder-kohlbrandt',
        'points': np.array([
            [562176,5931052],
            [562290,5931042],
            [562293,5931190],
            [562180,5931192]
        ])
    },
    # post sueder-norder unification
    {
        'name': 'post-unification',
        'points': np.array([
            [560484,5932946],
            [560485,5932765],
            [560640,5932764],
            [560637,5932946]
        ])
    }
]

transect_polygone = [
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

transect_planes = [
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
                "A_H": 0.01,
                "A_V": 0.001
            },
        "particle_properties": [
        ], 
        "trajectory_modifiers": [
            {
                "class_name": "oceantracker.trajectory_modifiers.resuspension.BasicResuspension",
                "critical_friction_velocity": 0.001
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
				"points": release_polygon[0]['points'],
				"pulse_size": pulse_size,
				"release_interval": pulse_interval
			}
        ],
        "particle_statistics": [
            {
                "class_name": "oceantracker.particle_statistics.polygon_statistics.PolygonStats2D_timeBased",
                "calculation_interval": 60,
                "count_status_equal_to": 'moving',
                "particle_property_list": [],
                "polygon_list": observational_polygone
            }
        ]
    }
}


#%%

runInfo = main.run(params)

 
#%%
path_to_dir = os.path.join(params['shared_params']['root_output_dir'],params['shared_params']['output_file_base'])
print(f'path_to_dir: {path_to_dir}')
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

m = load_output_files.load_stats_file(cases[0], nsequence=1)

#%%

tracks = load_output_files.load_particle_track_vars(cases[0],var_list=['water_depth'])
# plot_tracks.plot_tracks(tracks,plot_file_name='polygon_test.png',axis_lims=tracks['axis_lim'],polygon_list_to_plot=m['polygon_list'])

#%%
# plot_tracks.animate_particles(tracks,axis_lims=tracks['axis_lim'],movie_file='polygon_test.mp4')

import matplotlib.pyplot as plt
import shapely

area = [Polygon(item['points']).area for item in observational_polygone]

fig,ax = plt.subplots(2,1,figsize=(12,12))
fig.set_tight_layout(tight=True)
for jj in range(2):
    for ii in range(4):
        index = jj*4+ii+1
        print(index)
        print(area[index])
        print(m['count'][:,0,index])
        ax[jj].plot(m['count'][:,0,index]/area[index],label=m['polygon_list'][index]['name'])
        ax[jj].set_xlim(5000,10000)
        ax[jj].set_ylim(0,0.03)
fig.legend()
fig.savefig('stats_poly_test.png')


    
# %%


# %%
