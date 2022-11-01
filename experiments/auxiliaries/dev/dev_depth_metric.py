#%%
from re import A
import numpy as np
import matplotlib.pyplot as plt

from oceantracker.post_processing.read_output_files.load_output_files import get_case_info_files_from_dir, load_particle_track_vars, get_case_info_file_from_run_file
from oceantracker.post_processing.plotting.plot_vertical_tracks import plot_path_in_vertical_section, plot_relative_height
from oceantracker.post_processing.read_output_files import read_ncdf_output_files

#%%
def flat_n_clean(x):
    x = x.flatten()
    x = x[~np.isnan(x)]
    return x
#%%
path_to_dir = '/scratch/local1/output/22_08_19_testing_stranding_no_resus_v01'
case_info_file_name = get_case_info_files_from_dir(path_to_dir)[0]

#%%
tracks_file = '/scratch/local1/output/22_07_27_retention_v06/22_07_27_retention_v06_C143_tracks.nc'
var_list = ['x', 'time','status', 'age', 'tide', 'water_depth']
tracks = read_ncdf_output_files.read_particle_tracks_file(tracks_file, var_list)
# track_data = load_particle_track_vars(case_info_file_name, var_list=['tide', 'water_depth','age'])

min_status = 0
x = tracks['x']
sel = tracks['status'][:, :] < min_status
x[sel] = np.nan

below_free_surface = x[:,:,2] - tracks['tide']
# %%
month = 60*60*24*28
# %%
long_lived =  tracks['age'][-1] > 3*month

bfs_long = flat_n_clean(below_free_surface[:,long_lived])
bfs_short = flat_n_clean(below_free_surface[:,~long_lived])

plt.boxplot([bfs_long,bfs_short],whis=(0,100))

# %%
plt.hist(bfs_long,bins=100)

# %%
status = tracks['status'][:,long_lived]
status = status.flatten()
status = status[status > 0]
stranded_long = len(status[status == 3])/len(status)
print(f'{stranded_long} stranded ratio for long ')

status = tracks['status'][:,~long_lived]
status = status.flatten()
status = status[status > 0]
stranded_short = len(status[status == 3])/len(status)
print(f'{stranded_short} stranded ratio for short ')

status = tracks['status'][:,:5000]
status = status.flatten()
status = status[status > 0]
stranded_initial = len(status[status == 3])/len(status)
print(f'{stranded_initial} stranded ratio for short ')

# %%
# add distance traveled as a metric
