import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

from oceantracker.post_processing.read_output_files.load_output_files import load_grid
from oceantracker.post_processing.read_output_files.load_output_files import get_case_info_files_from_dir
from oceantracker.post_processing.read_output_files.load_output_files import read_case_info_file
from oceantracker.post_processing.read_output_files.load_output_files import load_track_data


def slice_tracks_by_iso8601(tracks, lower_threshold_iso8601, upper_threshold_iso8601):
    """
    Assumes two iso8601 strings as inputs
    """

    # generate datetime from threshold
    lower_threshold_datetime = datetime.datetime.fromisoformat(lower_threshold_iso8601)
    upper_threshold_datetime = datetime.datetime.fromisoformat(upper_threshold_iso8601)

    # tracks['time'] from posix to datetime 
    model_time_datetime = tracks['time'].astype('datetime64[s]').astype(datetime.datetime)

    # slice the tracks by datetime
    slice = np.logical_and(model_time_datetime >= lower_threshold_datetime, model_time_datetime <= upper_threshold_datetime)

    return slice


def average_part_prop_along_transect(transect,part_prop_name, nbins=10,time_slice=None):

    x = transect.track_data['x']
    
    part_prop = transect.track_data[part_prop_name]
    
    distance_downstream = x[:,:,0]

    alive = transect.track_data['status'] > 0
    particles_in_transect = ~np.isnan(x[:,:,0])
    in_and_alive = np.logical_and(particles_in_transect, alive)

    if time_slice is not None:
        time_slice = slice_tracks_by_iso8601(transect.track_data,time_slice[0],time_slice[1])
        in_and_alive = np.logical_and(time_slice[:,np.newaxis], in_and_alive)
        
    distance_downstream = distance_downstream[in_and_alive]
    part_prop = part_prop[in_and_alive]
    
    # bin the data along the transect and calculate the mean depth in each bin

    num_bins = nbins  # Example: 10 bins. Adjust this to your specific needs.
    bins = np.linspace(distance_downstream.min(), distance_downstream.max(), num_bins + 1)
    
    # Digitize the data to find out which bin each data point belongs to.
    bin_indices = np.digitize(distance_downstream, bins) - 1
    
    # Calculate the mean depth in each bin.
    mean_distance = np.array([distance_downstream[bin_indices == i].mean() for i in range(num_bins)])

    
    # print coount of particles in bin
    print([np.sum(bin_indices == i) for i in range(num_bins)])

    min_depths = np.zeros(num_bins)
    max_depths = np.zeros(num_bins)
    mean_depths = np.zeros(num_bins)
    std_depths = np.zeros(num_bins)

    for ii in range(num_bins):
        # print(f'bin {ii} has {np.sum(bin_indices == ii)} particles')
        print(np.sum(bin_indices == ii))
        if np.sum(bin_indices == ii) == 0:
            min_depths[ii] = np.nan
            max_depths[ii] = np.nan
            mean_depths[ii] = np.nan
            std_depths[ii] = np.nan
        else:
            min_depths[ii] = np.min(part_prop[bin_indices == ii])
            max_depths[ii] = np.max(part_prop[bin_indices == ii])
            mean_depths[ii] = np.mean(part_prop[bin_indices == ii])
            std_depths[ii] = np.std(part_prop[bin_indices == ii])
    
    # mean_depths = np.nan_to_num(mean_depths, nan=0.0)
    # std_depths = np.nan_to_num(std_depths, nan=0.0)
    
    # Return the mean depths. Alternatively, return bins and mean_depths for more detailed analysis.
    return mean_distance, (min_depths,max_depths), (mean_depths,std_depths)


def plot_part_prop_along_transect(transect, case, part_prop_name, part_prop_unit = '', time_slice=None):


    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = [ part_prop_name ]
    
    track_data = load_track_data(case,var_list)
    # print('Loaded data')

    transect.project_track_data(track_data)
    print('Projected data')
    mean_distance, (min_depths,max_depths), (mean_depths,std_depths) = average_part_prop_along_transect(transect,part_prop_name,nbins=100)
    print('Binned data')

    fig,ax = plt.subplots()

    ax.plot(mean_distance*1e-3, -max_depths,label=f'max. depth at transect',color='k')
    ax.plot(mean_distance*1e-3, -mean_depths,label=f'{initial_radius*1e6}um, {stickiness}')
    ax.fill_between(mean_distance*1e-3, -mean_depths-std_depths, -mean_depths+std_depths, alpha=0.3)

    ax.set_xlim(0,140)
    # ax.set_ylim(-30,2)
    ax.set_xlabel('distance downstream [km]')
    ax.set_ylabel(f'{part_prop_name} [{part_prop_unit}]')
    ax.set_title(f'average {part_prop_name} along transect')
    
    # legend to the right of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def average_depth_along_transect(transect,nbins=10,time_slice=None):

    x = transect.track_data['x']
    
    water_depth = transect.track_data['water_depth']
    tide = transect.track_data['tide']
    
    distance_downstream = x[:,:,0]
    depth_below_surface = tide - x[:,:,2]

    alive = transect.track_data['status'] > 0
    particles_in_transect = ~np.isnan(x[:,:,0])
    in_and_alive = np.logical_and(particles_in_transect, alive)

    if time_slice is not None:
        time_slice = slice_tracks_by_iso8601(transect.track_data,time_slice[0],time_slice[1])
        in_and_alive = np.logical_and(time_slice[:,np.newaxis], in_and_alive)
        
    distance_downstream = distance_downstream[in_and_alive]
    depth_below_surface = depth_below_surface[in_and_alive]
    
    # bin the data along the transect and calculate the mean depth in each bin

    num_bins = nbins  # Example: 10 bins. Adjust this to your specific needs.
    bins = np.linspace(distance_downstream.min(), distance_downstream.max(), num_bins + 1)
    
    # Digitize the data to find out which bin each data point belongs to.
    bin_indices = np.digitize(distance_downstream, bins) - 1
    
    # Calculate the mean depth in each bin.
    mean_distance = np.array([distance_downstream[bin_indices == i].mean() for i in range(num_bins)])

    
    # print coount of particles in bin
    print([np.sum(bin_indices == i) for i in range(num_bins)])

    # min_depths = np.array([depth_below_surface[bin_indices == i].min() for i in range(num_bins)])
    # max_depths = np.array([depth_below_surface[bin_indices == i].max() for i in range(num_bins)])
    # mean_depths = np.array([depth_below_surface[bin_indices == i].mean() for i in range(num_bins)])
    # std_depths = np.array([depth_below_surface[bin_indices == i].std() for i in range(num_bins)])

    min_depths = np.zeros(num_bins)
    max_depths = np.zeros(num_bins)
    mean_depths = np.zeros(num_bins)
    std_depths = np.zeros(num_bins)

    for ii in range(num_bins):
        # print(f'bin {ii} has {np.sum(bin_indices == ii)} particles')

        if np.sum(bin_indices == ii) == 0:
            min_depths[ii] = np.nan
            max_depths[ii] = np.nan
            mean_depths[ii] = np.nan
            std_depths[ii] = np.nan
        else:
            min_depths[ii] = np.min(depth_below_surface[bin_indices == ii])
            max_depths[ii] = np.max(depth_below_surface[bin_indices == ii])
            mean_depths[ii] = np.mean(depth_below_surface[bin_indices == ii])
            std_depths[ii] = np.std(depth_below_surface[bin_indices == ii])
    
    # mean_depths = np.nan_to_num(mean_depths, nan=0.0)
    # std_depths = np.nan_to_num(std_depths, nan=0.0)
    
    # Return the mean depths. Alternatively, return bins and mean_depths for more detailed analysis.
    return mean_distance, (min_depths,max_depths), (mean_depths,std_depths)


def average_relative_depth_along_transect(transect,nbins=10,time_slice=None):

    x = transect.track_data['x']
    
    water_depth = transect.track_data['water_depth']
    tide = transect.track_data['tide']
    
    distance_downstream = x[:,:,0]
    depth_below_surface = tide - x[:,:,2]
    total_water_depth = water_depth + tide

    alive = transect.track_data['status'] > 0
    particles_in_transect = ~np.isnan(x[:,:,0])
    in_and_alive = np.logical_and(particles_in_transect, alive)

    if time_slice is not None:
        time_slice = slice_tracks_by_iso8601(transect.track_data,time_slice[0],time_slice[1])
        in_and_alive = np.logical_and(time_slice[:,np.newaxis], in_and_alive)

    distance_downstream = distance_downstream[in_and_alive]
    depth_below_surface = depth_below_surface[in_and_alive]
    total_water_depth = total_water_depth[in_and_alive]
    relative_water_depth = depth_below_surface / total_water_depth
    # relative_water_depth should range from 0 to 1 with 0 being the surface and 1 being the bottom
    
    # bin the data along the transect and calculate the mean depth in each bin

    num_bins = nbins  # Example: 10 bins. Adjust this to your specific needs.
    bins = np.linspace(distance_downstream.min(), distance_downstream.max(), num_bins + 1)
    
    # Digitize the data to find out which bin each data point belongs to.
    bin_indices = np.digitize(distance_downstream, bins) - 1
    
    # Calculate the mean depth in each bin.
    mean_distance = np.array([distance_downstream[bin_indices == i].mean() for i in range(num_bins)])

    
    # print coount of particles in bin
    print([np.sum(bin_indices == i) for i in range(num_bins)])

    # min_depths = np.array([depth_below_surface[bin_indices == i].min() for i in range(num_bins)])
    # max_depths = np.array([depth_below_surface[bin_indices == i].max() for i in range(num_bins)])
    # mean_depths = np.array([depth_below_surface[bin_indices == i].mean() for i in range(num_bins)])
    # std_depths = np.array([depth_below_surface[bin_indices == i].std() for i in range(num_bins)])

    min_depths = np.zeros(num_bins)
    max_depths = np.zeros(num_bins)
    mean_depths = np.zeros(num_bins)
    std_depths = np.zeros(num_bins)

    for ii in range(num_bins):
        # print(f'bin {ii} has {np.sum(bin_indices == ii)} particles')

        if np.sum(bin_indices == ii) == 0:
            min_depths[ii] = np.nan
            max_depths[ii] = np.nan
            mean_depths[ii] = np.nan
            std_depths[ii] = np.nan
        else:
            min_depths[ii] = np.min(relative_water_depth[bin_indices == ii])
            max_depths[ii] = np.max(relative_water_depth[bin_indices == ii])
            mean_depths[ii] = np.mean(relative_water_depth[bin_indices == ii])
            std_depths[ii] = np.std(relative_water_depth[bin_indices == ii])
    
    # mean_depths = np.nan_to_num(mean_depths, nan=0.0)
    # std_depths = np.nan_to_num(std_depths, nan=0.0)
    
    # Return the mean depths. Alternatively, return bins and mean_depths for more detailed analysis.
    return mean_distance, (min_depths,max_depths), (mean_depths,std_depths)

def average_illumination_along_transect(transect,nbins=10,time_slice=None):
    x = transect.track_data['x']
    illumination = transect.track_data['illumination']
    
    distance_downstream = x[:,:,0]
    
    alive = transect.track_data['status'] > 0
    particles_in_transect = ~np.isnan(x[:,:,0])
    in_and_alive = np.logical_and(particles_in_transect, alive)

    if time_slice is not None:
        time_slice = slice_tracks_by_iso8601(transect.track_data,time_slice[0],time_slice[1])
        in_and_alive = np.logical_and(time_slice[:,np.newaxis], in_and_alive)

    distance_downstream = distance_downstream[in_and_alive]
    illumination = illumination[in_and_alive]
    
    # bin the data along the transect and calculate the mean depth in each bin

    num_bins = nbins  # Example: 10 bins. Adjust this to your specific needs.
    bins = np.linspace(distance_downstream.min(), distance_downstream.max(), num_bins + 1)
    
    # Digitize the data to find out which bin each data point belongs to.
    bin_indices = np.digitize(distance_downstream, bins) - 1
    
    # Calculate the mean depth in each bin.
    mean_distance = np.array([distance_downstream[bin_indices == i].mean() for i in range(num_bins)])
    mean_illu = np.array([illumination[bin_indices == i].mean() for i in range(num_bins)])
    std_illu = np.array([illumination[bin_indices == i].std() for i in range(num_bins)])

    # print coount of particles in bin
    print([np.sum(bin_indices == i) for i in range(num_bins)])
    
    mean_illu = np.nan_to_num(mean_illu, nan=0.0)
    std_illu = np.nan_to_num(std_illu, nan=0.0)
    
    # Return the mean depths. Alternatively, return bins and mean_depths for more detailed analysis.
    return mean_distance, (mean_illu,std_illu)    


def plot_depth_along_transect(transect, case):


    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = [ 'illumination',
                    'buoyancy',
                    'density',
                    'radius_spherical',
                    'tide',
                    'turbidity',
                    'water_depth']
    
    track_data = load_track_data(case,var_list)
    print('Loaded data')

    print('Sliced data')

    transect.project_track_data(track_data)
    print('Projected data')
    mean_distance, (min_depths,max_depths), (mean_depths,std_depths) = average_depth_along_transect(transect,nbins=100)
    print('Binned data')

    fig,ax = plt.subplots()

    ax.plot(mean_distance*1e-3, -max_depths,label=f'max. depth at transect',color='k')
    ax.plot(mean_distance*1e-3, -mean_depths,label=f'{initial_radius*1e6}um, {stickiness}')
    ax.fill_between(mean_distance*1e-3, -mean_depths-std_depths, -mean_depths+std_depths, alpha=0.3)

    ax.set_xlim(0,140)
    ax.set_ylim(-30,2)
    ax.set_xlabel('distance downstream [km]')
    ax.set_ylabel('depth [m]')
    ax.set_title(f'average depth along transect')
    
    # legend to the right of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def plot_relative_depth_along_transect(transect, case):


    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = [ 'illumination',
                    'buoyancy',
                    'density',
                    'radius_spherical',
                    'tide',
                    'turbidity',
                    'water_depth']
    
    track_data = load_track_data(case,var_list)
    print('Loaded data')

    print('Sliced data')

    transect.project_track_data(track_data)
    print('Projected data')
    mean_distance, (min_depths,max_depths), (mean_depths,std_depths) = average_relative_depth_along_transect(transect,nbins=100)
    print('Binned data')

    fig,ax = plt.subplots()

    # ax.plot(mean_distance*1e-3, -max_depths,label=f'max. depth at transect',color='k')
    ax.plot(mean_distance*1e-3, -mean_depths,label=f'{initial_radius*1e6}um, {stickiness}')
    ax.fill_between(mean_distance*1e-3, -mean_depths-std_depths, -mean_depths+std_depths, alpha=0.3)

    ax.set_xlim(0,140)
    ax.set_ylim(0,-1)
    ax.set_xlabel('distance downstream [km]')
    ax.set_ylabel('depth [m]')
    ax.set_title(f'average depth along transect')
    
    # legend to the right of plot
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


def rotate_points(x, y, label=None, angle = 32 * np.pi / 180, point = [5.489e5,5.9312e6]):

    x_rotated = (x - point[0]) * np.cos(angle) - (y - point[1]) * np.sin(angle) + point[0]
    y_rotated = (x - point[0]) * np.sin(angle) + (y - point[1]) * np.cos(angle) + point[1]

    if label is not None:
        for key, value in label.items():
            x_rotated_label = (value[0] - point[0]) * np.cos(angle) - (value[1] - point[1]) * np.sin(angle) + point[0]
            y_rotated_label = (value[0] - point[0]) * np.sin(angle) + (value[1] - point[1]) * np.cos(angle) + point[1]
            label[key] = [x_rotated_label, y_rotated_label]

        return x_rotated, y_rotated, label

    return x_rotated, y_rotated

def hexmap_of_death(case, property_name_and_unit='dead', labels_dict=None,
                            xlim=[4.80e5, 5.75e5], ylim=[5.924e6, 5.950e6], vmin=0, vmax=1000,
                            num_bins=30, fraction=1, save_path=None):

    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = ['cause_of_death']

    track_data = load_track_data(case,var_list)    
    
    x,y = track_data['x'][:,:,0],track_data['x'][:,:,1]
    dead = track_data['status'][-1] < 0
    particle_property = track_data['cause_of_death']

    # go thru all dead particles and find their last non nan location
    dead_particles = np.where(dead)[0]
    # get the last non nan location of each dead particle
    last_non_nan = np.zeros(len(dead_particles),dtype=int)
    for i,dead_particle in enumerate(dead_particles):
        # print(i, dead_particle)

        last_idx = np.where(np.isnan(track_data['x'][:,dead_particle]) == False)[0]
        if len(last_idx) == 0:
            last_non_nan[i] = -1
        else:
            last_non_nan[i] = last_idx[-1]

    x = track_data['x'][last_non_nan,dead_particles][:,0]
    y = track_data['x'][last_non_nan,dead_particles][:,1]
    particle_property = particle_property[last_non_nan,dead_particles]

    if labels_dict is not None:
        x, y, labels_rotated = rotate_points(x, y, labels_dict)
    else:
        x, y = rotate_points(x, y)

    # flatten
    x = x.flatten()
    y = y.flatten()
    particle_property = particle_property.flatten()

    # drop all nan values
    is_nan = np.isnan(x) | np.isnan(y) | np.isnan(particle_property)
    x = x[~is_nan]
    y = y[~is_nan]
    particle_property = particle_property[~is_nan]

    # select those in between axis limits
    in_window = np.logical_and(x >= xlim[0], x <= xlim[1])
    in_window = np.logical_and(in_window, y >= ylim[0])
    in_window = np.logical_and(in_window, y <= ylim[1])

    x = x[in_window]
    y = y[in_window]
    particle_property = particle_property[in_window]

    aspect_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))

    # slite down to a fraction of particles
    # create a random selection of indices to fit the length
    indices_subset = np.random.choice(x.shape[0], int(x.shape[0] * fraction), replace=False)

    x = x[indices_subset]
    y = y[indices_subset]
    particle_property = particle_property[indices_subset]

    # -- PLOTTING --

    # 12cm in width
    fig, ax = plt.subplots(figsize=(12, 12/(aspect_ratio*0.85)), dpi=300)

    # set title
    ax.set_title(f'{property_name_and_unit} of {initial_radius*1e6}um, stickiness: {stickiness}')

    ax.set_facecolor('lightgray')

    hb = ax.hexbin(x, y, 
                   gridsize=[int(num_bins*aspect_ratio), num_bins],
                   cmap=cm.get_cmap('inferno', 6), vmin=vmin, vmax=vmax,
                   mincnt=1,
                   extent=[xlim[0], xlim[1], ylim[0], ylim[1]], zorder=3)
    cb = fig.colorbar(hb, ax=ax, label=property_name_and_unit, pad=0.01)
    # add cbar label
    cb.ax.set_ylabel(property_name_and_unit, rotation=90, labelpad=10, fontsize=14)

    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # add an arrow pointing 32 degree north in the top center
    arrow_pos = [0.95, 0.1]
    ax.arrow(arrow_pos[0], arrow_pos[1], -np.sin(32 * np.pi / 180)*0.15, np.cos(32 * np.pi / 180)*0.15, 
             transform=ax.transAxes, length_includes_head=True, overhang=0.3,
             head_width=0.03, head_length=0.04, color='k', zorder=4)

    # add N to arrow
    ax.text(arrow_pos[0]+0.015, arrow_pos[1]-0.07, 'N', transform=ax.transAxes, fontsize=14, zorder=4)

    # add labels
    if labels_dict is not None:
        for key, value in labels_rotated.items():
            ax.text(value[0], value[1], key, ha='center', va='center', zorder=4, fontsize=14)

    # fix aspect ratio
    ax.set_aspect('equal', 'box')

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def hexmap_of_particle_prop(case, particle_property, property_name_and_unit='', labels_dict=None,
                            xlim=[4.80e5, 5.75e5], ylim=[5.924e6, 5.950e6], vmin=0, vmax=100,
                            num_bins=30, fraction=1, save_path=None):

    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = ['illumination','tide','turbidity', 'water_depth']

    track_data = load_track_data(case,var_list)
    
    x,y = track_data['x'][:,:,0],track_data['x'][:,:,1]
    particle_property = track_data[particle_property]

    if labels_dict is not None:
        x, y, labels_rotated = rotate_points(x, y, labels_dict)
    else:
        x, y = rotate_points(x, y)

    # flatten
    x = x.flatten()
    y = y.flatten()
    particle_property = particle_property.flatten()

    # drop all nan values
    is_nan = np.isnan(x) | np.isnan(y) | np.isnan(particle_property)
    x = x[~is_nan]
    y = y[~is_nan]
    particle_property = particle_property[~is_nan]

    # select those in between axis limits
    in_window = np.logical_and(x >= xlim[0], x <= xlim[1])
    in_window = np.logical_and(in_window, y >= ylim[0])
    in_window = np.logical_and(in_window, y <= ylim[1])

    x = x[in_window]
    y = y[in_window]
    particle_property = particle_property[in_window]

    aspect_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))

    # slite down to a fraction of particles
    # create a random selection of indices to fit the length
    indices_subset = np.random.choice(x.shape[0], int(x.shape[0] * fraction), replace=False)

    x = x[indices_subset]
    y = y[indices_subset]
    particle_property = particle_property[indices_subset]

    # -- PLOTTING --

    # 12cm in width
    fig, ax = plt.subplots(figsize=(12, 12/(aspect_ratio*0.85)), dpi=300)

    # set title
    ax.set_title(f'{property_name_and_unit} of {initial_radius*1e6}um, stickiness: {stickiness}')

    ax.set_facecolor('lightgray')

    hb = ax.hexbin(x, y, C=particle_property, reduce_C_function=np.mean,
                   gridsize=[int(num_bins*aspect_ratio), num_bins],
                   cmap=cm.get_cmap('viridis', 6), vmin=vmin, vmax=vmax,
                   extent=[xlim[0], xlim[1], ylim[0], ylim[1]], zorder=3)
    cb = fig.colorbar(hb, ax=ax, label=property_name_and_unit, pad=0.01)
    # add cbar label
    cb.ax.set_ylabel(property_name_and_unit, rotation=90, labelpad=10, fontsize=14)

    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # add an arrow pointing 32 degree north in the top center
    arrow_pos = [0.95, 0.1]
    ax.arrow(arrow_pos[0], arrow_pos[1], -np.sin(32 * np.pi / 180)*0.15, np.cos(32 * np.pi / 180)*0.15, 
             transform=ax.transAxes, length_includes_head=True, overhang=0.3,
             head_width=0.03, head_length=0.04, color='k', zorder=4)

    # add N to arrow
    ax.text(arrow_pos[0]+0.015, arrow_pos[1]-0.07, 'N', transform=ax.transAxes, fontsize=14, zorder=4)

    # add labels
    if labels_dict is not None:
        for key, value in labels_rotated.items():
            ax.text(value[0], value[1], key, ha='center', va='center', zorder=4, fontsize=14)

    # fix aspect ratio
    ax.set_aspect('equal', 'box')

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def hexmap_of_particle_prop(case, particle_property, property_name_and_unit='', labels_dict=None,
                            xlim=[4.80e5, 5.75e5], ylim=[5.924e6, 5.950e6], vmin=0, vmax=100,
                            num_bins=30, fraction=1, save_path=None):

    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = ['illumination','tide','turbidity', 'water_depth']

    track_data = load_track_data(case,var_list)
    
    x,y = track_data['x'][:,:,0],track_data['x'][:,:,1]
    particle_property = track_data[particle_property]

    if labels_dict is not None:
        x, y, labels_rotated = rotate_points(x, y, labels_dict)
    else:
        x, y = rotate_points(x, y)

    # flatten
    x = x.flatten()
    y = y.flatten()
    particle_property = particle_property.flatten()

    # drop all nan values
    is_nan = np.isnan(x) | np.isnan(y) | np.isnan(particle_property)
    x = x[~is_nan]
    y = y[~is_nan]
    particle_property = particle_property[~is_nan]

    # select those in between axis limits
    in_window = np.logical_and(x >= xlim[0], x <= xlim[1])
    in_window = np.logical_and(in_window, y >= ylim[0])
    in_window = np.logical_and(in_window, y <= ylim[1])

    x = x[in_window]
    y = y[in_window]
    particle_property = particle_property[in_window]

    aspect_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))

    # slite down to a fraction of particles
    # create a random selection of indices to fit the length
    indices_subset = np.random.choice(x.shape[0], int(x.shape[0] * fraction), replace=False)

    x = x[indices_subset]
    y = y[indices_subset]
    particle_property = particle_property[indices_subset]

    # -- PLOTTING --

    # 12cm in width
    fig, ax = plt.subplots(figsize=(12, 12/(aspect_ratio*0.85)), dpi=300)

    # set title
    ax.set_title(f'{property_name_and_unit} of {initial_radius*1e6}um, stickiness: {stickiness}')

    ax.set_facecolor('lightgray')

    hb = ax.hexbin(x, y, C=particle_property, reduce_C_function=np.mean,
                   gridsize=[int(num_bins*aspect_ratio), num_bins],
                   cmap=cm.get_cmap('viridis', 6), vmin=vmin, vmax=vmax,
                   extent=[xlim[0], xlim[1], ylim[0], ylim[1]], zorder=3)
    cb = fig.colorbar(hb, ax=ax, label=property_name_and_unit, pad=0.01)
    # add cbar label
    cb.ax.set_ylabel(property_name_and_unit, rotation=90, labelpad=10, fontsize=14)

    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # add an arrow pointing 32 degree north in the top center
    arrow_pos = [0.95, 0.1]
    ax.arrow(arrow_pos[0], arrow_pos[1], -np.sin(32 * np.pi / 180)*0.15, np.cos(32 * np.pi / 180)*0.15, 
             transform=ax.transAxes, length_includes_head=True, overhang=0.3,
             head_width=0.03, head_length=0.04, color='k', zorder=4)

    # add N to arrow
    ax.text(arrow_pos[0]+0.015, arrow_pos[1]-0.07, 'N', transform=ax.transAxes, fontsize=14, zorder=4)

    # add labels
    if labels_dict is not None:
        for key, value in labels_rotated.items():
            ax.text(value[0], value[1], key, ha='center', va='center', zorder=4, fontsize=14)

    # fix aspect ratio
    ax.set_aspect('equal', 'box')

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def hexmap_of_fraction_above_threshold(case,particle_property,
                            property_name_and_unit='',
                            threshold = 30,
                            labels_dict=None,
                            xlim=[4.80e5, 5.75e5], ylim=[5.924e6, 5.950e6],
                            vmin=0, vmax=1,
                            num_bins=30, fraction=1,
                            save_path=None):

    case_info = read_case_info_file(case)
    initial_radius = case_info['full_case_params']['class_dicts']['particle_properties']['radius_spherical']['initial_value']
    stickiness = case_info['full_case_params']['class_dicts']['particle_properties']['collision_very_fine_silt']['stickyness']

    var_list = ['illumination','tide','turbidity', 'water_depth']

    track_data = load_track_data(case,var_list)
    
    x,y = track_data['x'][:,:,0],track_data['x'][:,:,1]
    particle_property = track_data[particle_property]

    if labels_dict is not None:
        x, y, labels_rotated = rotate_points(x, y, labels_dict)
    else:
        x, y = rotate_points(x, y)

    # flatten
    x = x.flatten()
    y = y.flatten()
    particle_property = particle_property.flatten()

    # drop all nan values
    is_nan = np.isnan(x) | np.isnan(y) | np.isnan(particle_property)
    x = x[~is_nan]
    y = y[~is_nan]
    particle_property = particle_property[~is_nan]

    # select those in between axis limits
    in_window = np.logical_and(x >= xlim[0], x <= xlim[1])
    in_window = np.logical_and(in_window, y >= ylim[0])
    in_window = np.logical_and(in_window, y <= ylim[1])

    x = x[in_window]
    y = y[in_window]
    particle_property = particle_property[in_window]


    aspect_ratio = (np.max(x) - np.min(x)) / (np.max(y) - np.min(y))

    # slite down to a fraction of particles
    # create a random slection of indices to fit the length
    indices_subset = np.random.choice(x.shape[0], int(x.shape[0] * fraction), replace=False)

    x = x[indices_subset]
    y = y[indices_subset]
    particle_property = particle_property[indices_subset]

    # fraction of particles with particle_property value above threshold
    above_threshold = particle_property > threshold

    # create reduce_C_function to color hexbins by fraction of particles above threshold
    def reduce_C_function_above_threshold(values):
        # print(type(values))
        # print(values)
        # print( np.sum(np.array(values) > threshold),
        #        np.array(values).shape[0],
        #        np.sum(np.array(values) > threshold) / np.array(values).shape[0]
        # )
        return np.sum(values) / len(values)
    


    # -- PLOTTING --

    # 12cm in width
    fig,ax = plt.subplots(figsize=(12 , 12/(aspect_ratio*0.85)), dpi=300)

    # set title
    ax.set_title(f'{property_name_and_unit} of {initial_radius*1e6}um, stickiness: {stickiness}')

    ax.set_facecolor('lightgray')


    hb = ax.hexbin(x, y, C=above_threshold, reduce_C_function=reduce_C_function_above_threshold,
               gridsize=[int(num_bins*aspect_ratio),num_bins],
               cmap=cm.get_cmap('viridis',6), vmin=vmin, vmax=vmax,
               extent=[xlim[0], xlim[1], ylim[0], ylim[1]], zorder=3)
    cb = fig.colorbar(hb, ax=ax, label=property_name_and_unit,pad=0.01)
# add cbar label
    cb.ax.set_ylabel(property_name_and_unit, rotation=90, labelpad=10,fontsize=14)


# Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

# Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

# add an arrow pointing 32 degree north in the top center
    arrow_pos = [0.95, 0.1]
    ax.arrow(arrow_pos[0],arrow_pos[1], -np.sin(32 * np.pi / 180)*0.15, np.cos(32 * np.pi / 180)*0.15, 
         transform=ax.transAxes, length_includes_head=True, overhang=0.3,
         head_width=0.03, head_length=0.04, color='k', zorder=4)

# add N to arrow
    ax.text(arrow_pos[0]+0.015,arrow_pos[1]-0.07, 'N', transform=ax.transAxes, fontsize=14, zorder=4)

# add labels
    if labels_dict is not None:
        for key, value in labels_rotated.items():
            ax.text(value[0], value[1], key, ha='center', va='center', zorder=4, fontsize=14)

# fix aspect ratio
    ax.set_aspect('equal', 'box')


    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')


def animate_full_transect(tracks, transect, t_min=0, t_duration=None, plot_file_name=None, axis_lims_map=None, axis_lims_transect_h=None, axis_lims_transect_v=None):
    """
    Creates three subplots:
    - horizontal tracks on the unprojected model domain
    - vertical tracks on the projected model domain
    - horizontal tracks on the projected model domain
    """

    def update(nt):
        ax[0].clear()
        ax[1].clear()

        x = transect.track_data['x']

        alive = transect.track_data['status'][nt] > 0
        particles_in_transect = ~np.isnan(x[nt, :, 0])
        in_and_alive = np.logical_and(particles_in_transect, alive)

        x = x[nt, in_and_alive, :]
        depth = transect.track_data['water_depth'][nt, particles_in_transect]
        water_level = transect.track_data['tide'][nt, particles_in_transect]


        ## MAP
        tri01 = ax[0].tripcolor(transect.track_data['grid']['x'][:,0],
                    transect.track_data['grid']['x'][:,1],
                    transect.track_data['grid']['water_depth'],
                    triangles=transect.track_data['grid']['triangles'],
                    shading='gouraud', cmap='Blues', edgecolors='none', zorder=1)

        scatter01 = ax[0].scatter(tracks['x'][nt, particles_in_transect, 0], tracks['x'][nt, particles_in_transect, 1], c=np.arange(tracks['x'][:,particles_in_transect,:].shape[1]) % 20)


        ## VERTICAL TRANSECT
        # particle
        scatter21 = ax[1].scatter(x[:, 0], x[:, 2], c=np.arange(x.shape[0]) % 20)
        # bottom at particle location
        scatter22 = ax[1].scatter(x[:, 0], -depth, marker='_', c=np.arange(x.shape[0]) % 20)
        # water level at particle location
        scatter23 = ax[1].scatter(x[:, 0], water_level, marker='_', c=np.arange(x.shape[0]) % 20)
        # line from water level to bottom set behing the particle scatter plot
        lines20 = ax[1].vlines(x[:, 0], -depth, water_level, alpha=0.5)
        # ax[0].set_title(f'verticle tracks at timestep {transect.track_data["time"][nt].astype("datetime64[s]")}')

        ## AXIS LIMS   
        ax[0].set_title(f'map at timestep {transect.track_data["time"][nt].astype("datetime64[s]")}')

        ax[0].set_title(f'position in transect at timestep {transect.track_data["time"][nt].astype("datetime64[s]")}')
        ax[1].set_xlabel('distance downstream [m]')
        ax[1].set_ylabel('depth [m]')
        # swap orientation of x axis
        
        ax[1].set_xlim((0, transect.transect[-1]['length'] + transect.transect[-1]['distance_downstream']))
        ax[1].invert_xaxis()

        if axis_lims_map is not None:
            ax[0].set_xlim(axis_lims_map[:2])
            ax[0].set_ylim(axis_lims_map[2:])
        if axis_lims_transect_v is not None:
            ax[1].set_ylim(axis_lims_transect_v[2:])   

        return scatter01, scatter21, scatter22, scatter23, lines20
    

    fig, ax = plt.subplots(2,1,figsize=(10,15),dpi=100)

    ani = FuncAnimation(fig, update,
                        frames=range(t_min, t_min + (t_duration if t_duration else len(transect.track_data['x']) - t_min), 1),
                        blit=True)

    if plot_file_name is not None:
        ani.save(plot_file_name, writer='ffmpeg', fps=10)
    else:
        plt.show()
