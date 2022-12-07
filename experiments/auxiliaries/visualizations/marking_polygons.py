
from oceantracker.post_processing.plotting.plot_utilities import *

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
# import matplotlib.font_manager as font_manager

# from oceantracker.util.triangle_utilities_code import convert_face_to_nodal_values

from oceantracker.post_processing.read_output_files import load_output_files
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# from matplotlib import animation
# from oceantracker.util import time_util

color_palette={'land': (np.asarray([146, 179, 140])/256).tolist(), 'land_edge': [.5, .5, .5]}

def clickable_map(grid, axis_lims=None, back_ground_depth=True,
                  show_grid=True, back_ground_color_map='Blues', title=None, text1=None, credit=None):

    # get grid bounds to fill a recgtangle
    bounds= [np.min(grid['x'][:, 0]), np.max(grid['x'][:, 0]), np.min(grid['x'][:, 1]), np.max(grid['x'][:, 1])]
    dx,dy = bounds[1]- bounds[0], bounds[3]- bounds[2]
    f= 0.05
    bounds =np.asarray([ [bounds[0]-f*dx, bounds[1]+f*dx], [bounds[2]-f*dy,  bounds[3]+f*dy]]) # l
    b = np.asarray([bounds[0,:], [bounds[1,0], bounds[0,1] ], bounds[1, :], [bounds[0,0],bounds[1,1] ], bounds[0,:] ] )


    fig, ax = plt.subplots(1, 1)

    # fill background land retangle
    ax.fill(b[:,0] , b[:, 1],  facecolor=color_palette['land'],  zorder=0)

    if axis_lims is None: axis_lims= bounds.flatten().tolist()
    ax.set_xlim(axis_lims[:2])
    ax.set_ylim(axis_lims[2:])

    # fill domain as white
    ax.fill(grid['grid_outline']['domain']['points'][:,0], grid['grid_outline']['domain']['points'][:,1],
            edgecolor= None, facecolor=(1., 1., 1.), linewidth=.5, zorder=0)

    # plot islands from outline
    for g in grid['grid_outline']['islands']:
            ax.fill(g['points'][:, 0], g['points'][:, 1], edgecolor=color_palette['land_edge'],
                    facecolor=color_palette['land'], linewidth= 0.5, zorder= 3)
    ax.plot(grid['grid_outline']['domain']['points'][:, 0], grid['grid_outline']['domain']['points'][:, 1], c=color_palette['land_edge'], linewidth=0.5, zorder=3)

    if  back_ground_depth:
        plot_coloured_depth(grid, ax=ax,color_map= back_ground_color_map,zorder=1)

    if show_grid:
        ax.triplot(grid['x'][:, 0], grid['x'][:, 1], grid['triangles'], color=(0.8, 0.8, 0.8), linewidth=.5, zorder=1)
    for o in grid['grid_outline']['open_boundary_nodes']:
        plt.plot(grid['x'][o, 0], grid['x'][o, 1], '--','red')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(axis="both", direction="in", right=True, top=True)

    if title is not None:  ax.set_title(title)
    if text1 is not None:  text_norm(.4, .1, text1, fontsize=8)
    add_credit(credit)
    add_map_scale_bar(axis_lims, ax=ax)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show()

# coords = []

def onclick(event):
    x,y = event.xdata, event.ydata
    print(f'[{round(x)},{round(y)}]')


#path_to_grid = '/scratch/local1/output/21_11_01_sa_verticle_spawning_v03/21_11_01_sa_verticle_spawning_v03_grid.nc'
# path_to_grid = '/scratch/local1/output/22_04_22_light_dryness_frictionresusp_tests_v01/22_04_22_light_dryness_frictionresusp_tests_v01_runInfo.json'
path_to_grid = '/scratch/local1/output/22_08_19_testing_stranding_no_resus_v01/22_08_19_testing_stranding_no_resus_v01_runInfo.json'

grid = load_output_files.load_grid(path_to_grid)

clickable_map(grid)