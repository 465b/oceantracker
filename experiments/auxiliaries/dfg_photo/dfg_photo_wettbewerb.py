#%%
from oceantracker.post_processing.read_output_files.load_output_files import load_particle_track_vars
from oceantracker.post_processing.read_output_files.load_output_files import get_case_info_files_from_dir
from oceantracker.post_processing.read_output_files.load_output_files import read_case_info_file
from oceantracker.post_processing.read_output_files.load_output_files import load_grid
from oceantracker.post_processing.read_output_files.load_output_files import _extract_useful_params
from oceantracker.post_processing.plotting import plot_utilities
from oceantracker.post_processing.read_output_files import read_ncdf_output_files

from matplotlib import colors, animation
from oceantracker.util import time_util
from oceantracker.util import cord_transforms
import pyproj

import numpy as np
import matplotlib.pyplot as plt

#%%
path_to_dir = '/scratch/local1/output/22_08_19_testing_stranding_no_resus_v01'
cases = get_case_info_files_from_dir(path_to_dir)
case_info = read_case_info_file(cases[0])
# track_data = load_particle_track_vars(cases[0])

tracks_file = '/scratch/local1/output/22_07_27_retention_v06/22_07_27_retention_v06_C143_tracks.nc'
var_list = list(set(['time', 'x','status'])) # default vars
tracks = read_ncdf_output_files.read_particle_tracks_file(tracks_file, var_list)

grid_file = '/scratch/local1/output/22_07_27_retention_v06/22_07_27_retention_v02_grid.nc'
d = read_ncdf_output_files.read_grid_file(grid_file)

# load  grid outline and convert outline to numpy arrays
grid_outline_file = '/scratch/local1/output/22_07_27_retention_v06/22_07_27_retention_v02_grid_outline.json'
d['grid_outline'] = read_ncdf_output_files.read_grid_outline_file(grid_outline_file)

# make outline list np arrays for plotting
for key in d['grid_outline']['domain']:
    d['grid_outline']['domain'][key] = np.asarray(d['grid_outline']['domain'][key])
for n in range(len(d['grid_outline']['islands'])):
    for key in  d['grid_outline']['islands'][n]:
        d['grid_outline']['islands'][n][key]= np.asarray( d['grid_outline']['islands'][n][key])
    
tracks['grid'] = d

tracks.update({'particle_status_flags': case_info['particles']['particle_status_flags'],
            'particle_release_groups': case_info['particle_release_groups']})
# tracks['full_params'] = case_info['full_params']


track_data = tracks

# tracks= _extract_useful_params(case_info, tracks)
# %%

color_palette={'land': (np.asarray([146, 179, 140,0])/256).tolist(), 'land_edge': [.6, .6, .6]}

def plot_field(grid, field_vals, ax=plt.gca(), color_map=None, vmin=None, vmax=None, zorder=3):
    # use tri surf to color map feild in 3D, defaul view is 2D from above

    ax.tripcolor(grid['x'][:,0], grid['x'][:,1], field_vals,
                 triangles=grid['triangles'], alpha=0.6,
                  shading='gouraud', cmap=color_map, edgecolors='none',
                  vmin=vmin, vmax=vmax, zorder=zorder)

def plot_coloured_depth(grid, ax=plt.gca(), color_map=None, zorder=3):
    # find depth range inside axes to set max  and min depth

    # plot colored depth, but dilute deepest colour but setting vmax 20% larger, to set colormap limits based on nodes inside axies
    sel =  np.logical_and(grid['x'][:,0] >= ax.get_xlim()[0],  grid['x'][:,0] <= ax.get_xlim()[1])
    sel = np.logical_and(sel, grid['x'][:, 1] >= ax.get_ylim()[0])
    sel = np.logical_and(sel, grid['x'][:, 1] <= ax.get_ylim()[1])
    depth = grid['water_depth']
    vmax = np.nanmax(depth[sel])

    # blue particle can be lost in deep water, cant get alphaand  not to draw grid edges, so outscale deepest colour
    if color_map is None:
        color_map = 'Blues'
        vmax=1.3*vmax

    plot_field(grid, depth, ax=ax, color_map=color_map, vmin= 0., vmax=vmax, zorder=zorder)

def draw_base_map(grid, osm_img_path,
                  ax=plt.gca(), axis_lims=None, background_lims=None, back_ground_depth=True,
                  show_grid=False, back_ground_color_map='Blues', title=None, text1=None, credit=None):

    # get grid bounds to fill a recgtangle
    bounds= [np.min(grid['x'][:, 0]), np.max(grid['x'][:, 0]), np.min(grid['x'][:, 1]), np.max(grid['x'][:, 1])]
    dx,dy = bounds[1]- bounds[0], bounds[3]- bounds[2]
    f= 0.05
    bounds =np.asarray([ [bounds[0]-f*dx, bounds[1]+f*dx], [bounds[2]-f*dy,  bounds[3]+f*dy]]) # l
    b = np.asarray([bounds[0,:], [bounds[1,0], bounds[0,1] ], bounds[1, :], [bounds[0,0],bounds[1,1] ], bounds[0,:] ] )

    if axis_lims is None: axis_lims= bounds.flatten().tolist()
    ax.set_xlim(axis_lims[:2])
    ax.set_ylim(axis_lims[2:])

    print(axis_lims)
    # print(cord_transforms.WGS84_to_ETRS89(np.array([5.92e6,4.4e5])))
    
    # fill background land retangle
    osm_map = plt.imread(osm_img_path)
    if background_lims is not None:
        ax.imshow(osm_map, zorder=0, extent=background_lims, aspect= 'equal',alpha=0.7)
    else:
        ax.imshow(osm_map, zorder=0, aspect= 'equal')



    # fill domain as white
    # ax.fill(grid['grid_outline']['domain']['points'][:,0], grid['grid_outline']['domain']['points'][:,1],
    #         edgecolor= None, facecolor=(1., 1., 1.), linewidth=.5, zorder=0)

    # plot islands from outline
    for g in grid['grid_outline']['islands']:
        ax.fill(g['points'][:, 0], g['points'][:, 1], edgecolor=color_palette['land_edge'],
                facecolor=color_palette['land'], linewidth= 1, zorder= 3)
    ax.plot(grid['grid_outline']['domain']['points'][:, 0],
            grid['grid_outline']['domain']['points'][:, 1], c=color_palette['land_edge'], 
            linewidth=2, zorder=3)

    if  back_ground_depth:
        plot_coloured_depth(grid, ax=ax,color_map= back_ground_color_map,zorder=1)

    if show_grid:
        ax.triplot(grid['x'][:, 0], grid['x'][:, 1], grid['triangles'], color=(0.8, 0.8, 0.8), linewidth=.5, zorder=1)
    for o in grid['grid_outline']['open_boundary_nodes']:
        plt.plot(grid['x'][o, 0], grid['x'][o, 1], '--','red')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis="both", direction="in", right=True, top=True)

    if title is not None:  ax.set_title(title)
    if text1 is not None:  text_norm(.4, .1, text1, fontsize=8)
    plot_utilities.add_credit(credit)
    plot_utilities.add_map_scale_bar(axis_lims, ax=ax)

    return grid

def draw_dfg_img(track_data, osm_img, nt,
                 axis_lims=None, background_lims=None, colour_using_data= None, show_grid=False, title=None, max_duration=None,
                 movie_file= None,  dpi=300, size=8,
                 min_status=0, back_ground_depth=True, back_ground_color_map = None, credit=None, heading= None,
                 size_using_data= None,  part_color_map=None,
                 vmin=None, vmax=None,
                 release_group=None):


    fig = plt.figure(
        figsize=(23.3,10),
        dpi=300,
        frameon=False)

    ax = plt.gca()
    draw_base_map(
        track_data['grid'],osm_img, ax=ax, axis_lims=axis_lims,background_lims=background_lims,
        show_grid=show_grid, title=title, credit=credit, back_ground_depth=back_ground_depth,
        back_ground_color_map=back_ground_color_map)
        
    s0 = size
        
    # colour by status
    colour_using_data = np.full_like(track_data['status'], -127)

    stat_types = track_data['particle_status_flags']
    status_list = [stat_types['outside_open_boundary'],stat_types['dead'],  stat_types['frozen'], stat_types['stranded_by_tide'], stat_types['on_bottom'], stat_types['moving']]
    for n, val in enumerate(status_list):
        colour_using_data[track_data['status']==val] = n  #replace status with range(status_list)

    colour_using_data = colour_using_data.astype(np.float64)
    status_colour_map = np.asarray([[0.9961,    0.8906,    0.7070],[0, 0., 0.],[.8, 0, 0.],  [0, .5, 0.], [0.5, 0.5, .5], [0, 0, 1.] ])
    cmap = colors.ListedColormap(status_colour_map)

    sc = ax.scatter(track_data['x'][nt, :, 0], track_data['x'][nt, :, 1], c=colour_using_data[nt, :], vmin=0,
                    vmax=len(status_list), s=s0, edgecolors=None, cmap=cmap, zorder=5)
    ax.scatter([],[],color='blue',label='drifting plankton')                    
    ax.scatter([],[],color='green',label='stranded plankton')                    

    # only plot alive particles
    x = track_data['x'][nt, :, :2].copy() # copy so as not to change original data
    sel = track_data['status'][nt, :] < min_status # get rid of dead particles
    x[sel,:] = np.nan
    sc.set_offsets(x)
    sc.set_array(colour_using_data[nt, :].astype(np.float64))
    sc.set_zorder(5)
    if size_using_data is not None: sc.set_sizes(scaled_marker_size[nt, :])
    time_text = plt.text(.05, .05, time_util.seconds_to_pretty_str(track_data['time'][0], seconds=False), transform=ax.transAxes)
    time_text.set_text(time_util.seconds_to_pretty_str(track_data['time'][nt], seconds=False))
    
    plot_utilities.add_heading(heading)
    plot_utilities.show_particleNumbers(track_data['x'].shape[1])
    
    # lg = ax.legend(
    #     loc='upper left',
    #     bbox_to_anchor=(0.6, 0.95),
    #     fontsize='xx-large',
    #     prop={'family':'serif','size':'xx-large'})
    # lg.set_title(
    #     """Phytoplankton Retention Experiment\nv0.6 performed at 22nd of August '22\n\nParticles represent phytoplankton modeled\nin the Elbe estuary  while the two colors\nindicate the current status of the plankton.\n""",
    #              prop={'family':'serif','size':'xx-large'})
    # lg.get_title().set_multialignment('center')
    fig.tight_layout()
    
    plt.savefig(('/home/zmaw/u301513/Documents/scr/phd/bicest/oceantracker/experiments/auxiliaries/dfg.png'))

#%%

axis_lims = [435445.1921875, 600900.9015625, 5911947.75, 6012889.25]

x_shift = -14200
y_shift = -100
x_stretch = 1.115
y_stretch = 1.05

x_left = 435445.1921875
x_right = 600900.9015625
x_width = x_right - x_left

y_bottom = 5911947.75
y_top = 6012889.25
y_hight = y_top - y_bottom

x_left += x_shift
x_right = x_left + x_width*x_stretch

y_top += y_shift
y_bottom = y_top - y_hight*y_stretch


background_lims = [x_left,x_right,y_bottom,y_top]

draw_dfg_img(track_data, '/home/zmaw/u301513/Documents/scr/phd/bicest/oceantracker/experiments/auxiliaries/dfg_photo_osm.png', 
             30,axis_lims=axis_lims,background_lims=background_lims,
             credit='Laurin Steidle & Ross Vennell 2022\nOceantracker')


# %%

# %%

# %%
