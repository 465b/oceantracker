from oceantracker.run_oceantracker import main
from oceantracker.util import yaml_util, json_util
from os import path
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-norun', action='store_true')
parser.add_argument('-hpc', action='store_true')
args = parser.parse_args()

# https://tidesandcurrents.noaa.gov/ofs/lsofs/lsofs.html
# https://www.ncei.noaa.gov/thredds/catalog/model-lsofs-files/catalog.html
output_file_base='FVCOM_Lake_Superior_test'


if args.hpc:
    input_dir='/hpcfreenas/hindcast/LakeSuperior/'
    root_output_dir = '/hpcfreenas/ross/oceanTrackerOutput/LakeSuperior/'
else:
    input_dir=r'F:\Hindcasts\colaborations\LakeSuperior\historical_sample\2022'
    root_output_dir = r'F:\OceanTrackerOutput'

print(args,input_dir)
file_mask ='nos.lsofs.fields.n000*.nc'

points = [[256203.6793068961, 5193002.88896844, -10],
           [416692.1617094234, 5216000.828769726, -10],
           [666422.5233426465, 5189371.635315605, -10],
           [429178.67979108455, 5417656.448290474, -10],
           [439094.44415005075, 5265627.962025132, -10]]

# just use one point
points= [[439094.44415005075, 5265627.962025132, -10]]

params={'output_file_base' : output_file_base,
        'regrid_z_to_uniform_sigma_levels': False,
        'root_output_dir':root_output_dir,
        'time_step' : 20*60,
        'reader': {"class_name": 'oceantracker.reader.FVCOM_reader.unstructured_FVCOM',
                'input_dir': input_dir,
                'file_mask': file_mask},
        'user_note':'test of notes',
        'release_groups': {'mytest_points': {'points': points, 'pulse_size': 250, 'release_interval': 7200}} ,
        'particle_statistics':{'example_gridded_stat' :
                  {'class_name': 'oceantracker.particle_statistics.gridded_statistics.GriddedStats2D_timeBased',
                      'update_interval': 72000,   'grid_size': [320, 321],'grid_span':[ 250000,250000],'grid_center':points[0]}
                    }

}

#yaml_util.write_YAML(output_file_base+'.yaml',params)
#json_util.write_JSON(output_file_base+'.json',params)

if args.norun:
    # infer run file name
    caseInfo_file_name = path.join(params['root_output_dir'], params['output_file_base'], params['output_file_base'] + '_caseInfo.json')
else:
    # run oceantracker
    caseInfo_file_name= main.run(params)


from oceantracker.post_processing.read_output_files import load_output_files
from oceantracker.post_processing.plotting import plot_utilities
from oceantracker.post_processing.plotting.plot_tracks import animate_particles, plot_tracks
from oceantracker.post_processing.plotting.plot_statistics import animate_heat_map, plot_heat_map

grid= load_output_files.load_grid(caseInfo_file_name)

#plot_utilities.display_grid(grid, ginput=6)

track_data = load_output_files.load_track_data(caseInfo_file_name,fraction_to_read=.1)

animate_particles(track_data,  show_grid=True,axis_lims=None,
                  heading='FVCOM reader test',show_dry_cells=False,
                  release_group=None,  fps=15,size=6,
                  back_ground_depth=True, interval=20)
plot_tracks(track_data)

# heat maps from on the fly counts
stats_data = load_output_files.load_stats_data(caseInfo_file_name)

animate_heat_map(stats_data,'mytest_points',  heading=output_file_base + ' particle count heat map',
                 vmax=100.)
plot_heat_map(stats_data,'mytest_points',  heading=output_file_base + ' particle count heat map', vmax=100.)