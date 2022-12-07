# %%
# from oceantracker.oceanTrackerMain import run_oceantracker
from oceantracker.user_post_processing import loadOutputFiles
from oceantracker.user_post_processing import statsPlot
from oceantracker.user_post_processing.particlePlot import * 
from oceantracker.util import otTime
import numpy as np
import os
# %%

def draw_grid(grid, show_grid=True):

    # get grid bounds to fil a recgtangle
    bounds= [np.min(grid['x'][:, 0]), np.max(grid['x'][:, 0]), np.min(grid['x'][:, 1]), np.max(grid['x'][:, 1])]
    dx,dy = bounds[1]- bounds[0], bounds[3]- bounds[2]
    f= 0.1
    b =[ [bounds[0]-f*dx, bounds[2]-f*dy], [bounds[1]+f*dx,  bounds[3]+f*dy]] # l
    b = np.asarray([b[0],[b[0][0],b[1][0] ], b[1], [b[0][0],b[1][1] ], b[0] ] )

    plt.fill(b[:,0] , b[:,1 ],  facecolor=(.9,.9,.9),  zorder=0)



    # fill domain as white
    plt.fill(grid['grid_outline']['domain']['points'][:,0], 
             grid['grid_outline']['domain']['points'][:,1],
             edgecolor=(1,1,1), facecolor=(.9,.9,.9), linewidth=.5, zorder=0)
    # plot islands from outline
    for g in grid['grid_outline']['islands']:
            plt.fill(g['points'][:, 0], g['points'][:, 1], 
            edgecolor=(0.8, 0.8, 0.8), facecolor=(.9,.9,.9), linewidth= 0.5, zorder= 0)
    if show_grid:
        plt.tripcolor(grid['x'][:,0,],grid['x'][:,1],water_depth,triangles=grid['triangles'],cmap='magma',zorder=0)
        #plt.triplot(grid['x'][:, 0], grid['x'][:, 1], grid['triangles'], color=(0.8, 0.8, 0.8), linewidth=.5, zorder=0)


    text_norm(.6, .05, 'OceanTracker- R. Vennell, 2022', fontsize=5)

def animate_particles(runInfo_fileName_or_runCaseInfoDict, axes=None, num_to_plot=10 ** 2, show_grid=True, title=None,
                      movie_file= None, fps=15, dpi=300, ncase=0):

    def draw_frame(nt):

        sel =status[nt, :] == stat_types['moving']
        sc1.set_offsets(x[nt,sel, :] )

        sel = status[nt, :] == stat_types['stranded_byTide']
        sc2.set_offsets(x[nt, sel, :])

        sel = status[nt, :] == stat_types['stranded_onBottom']
        sc3.set_offsets(x[nt, sel, :])

        sel = status[nt, :] == stat_types['frozen']
        #sc4.set_offsets(x[nt, sel, :])

        sel = status[nt, :] < stat_types['frozen']
        #sc5.set_offsets(x[nt, sel, :])

        title_text.set_text(otTime.seconds_to_pretty_str(tracks['time'][nt], seconds=False))

        return  sc1,sc2,sc3,sc4,sc5,title_text

    runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(runInfo_fileName_or_runCaseInfoDict)

    tracks = loadOutputFiles.load_particle_track_vars(runCaseInfo, ['x', 'status', 'time'], ncase=ncase)

    grid = loadOutputFiles.load_grid(runCaseInfo)

    status = tracks['status']
    stat_types = runCaseInfo['particle_status_flags']

    x = tracks['x'][:, :, :2]
    #fig = plt.gcf()
    fig = plt.figure(figsize=(7,5))
    ax = plt.gca()


    set_axes(axes, x=tracks['x'][:, :, 0], y=tracks['x'][:, :, 1])
    draw_grid(grid,show_grid=show_grid)
    plot_release_points_and_polygons(runCaseInfo, ncase= ncase)
    x0 = x[0,:,:]*np.nan # fill with NaNs to start

    #moving
    sc1 = ax.scatter(x0[:, 0], x0[:, 1], c='lime', s=2, zorder=9)
    #stranded
    sc2 = ax.scatter(x0[ :, 0], x0[ :, 1], c='green', s=2, zorder =8)
    #onbottom
    sc3 = ax.scatter(x0[:, 0], x0[:, 1], c='gray', s=2, zorder =7)
    #outside domain - i think
    sc4 = ax.scatter(x0[ :, 0], x0[ :, 1], c='r', s=2,zorder =6)
    #dead
    sc5 = ax.scatter(x0[:, 0], x0[:, 1], c='black', s=2, zorder=5)

    title_text = plt.text(.05, .05, otTime.seconds_to_pretty_str(tracks['time'][0], seconds=False), transform=ax.transAxes)


    show_particleNumbers(tracks)

    if title is not None:  ax.set_title(title)

    fig.tight_layout()

    anim = animation.FuncAnimation(fig, draw_frame, frames=np.arange(0,x.shape[0],20), interval=1, blit=True)

    animation_output(anim, movie_file, fps=fps, dpi=dpi)

    return
# %%
path = '/scratch/local1/output/22_04_22_light_dryness_frictionresusp_tests_v01/22_04_22_light_dryness_frictionresusp_tests_v01_runInfo.json'
runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path)
grid = loadOutputFiles.load_grid(runCaseInfo)

# for future use:
# make this more sophisticated and read nc automatically based on runCaseInfo
#file_mask = runCaseInfo['user_params']['reader']['file_mask']
#input_dir = runCaseInfo['user_params']['shared_params']['input_dir']
# %%
from netCDF4 import Dataset
df = Dataset('/scratch/local1/hzg2/schout_1.nc')
# %%
topo = df.variables['depth']
water_level = df.variables['elev']
water_depth = water_level[0] - topo
water_depth[water_depth > 0] = 0
water_depth[water_depth < -20] = -20

#%%
animate_particles(runCaseInfo,movie_file='animate_with_depth_coloring.mp4',axes=[520e3,590e3,5.915e6,5.97e6],ncase=3)
# %%

# %%
