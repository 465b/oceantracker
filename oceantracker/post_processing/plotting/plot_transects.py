import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from oceantracker.post_processing.plotting.plot_utilities import draw_base_map
from oceantracker.post_processing.plotting import plot_utilities 

#%%
def plot_transect_map(transect_class, nt=0, min_status = 0,
                      show_grid=False, credit=None,
                      heading =None,title=None, axis_lims=None, back_ground_depth=True,
                      back_ground_color_map= None,plot_file_name=None, 
                      polygon_list_to_plot = None):

    transect = transect_class.transect

    grid = transect_class.track_data['grid']
    track_data = transect_class.track_data

    fig = plt.figure()
    ax = plt.gca()

    fig.tight_layout()

    draw_base_map(grid, ax=ax, axis_lims=axis_lims, show_grid= show_grid, title=title, credit=credit,
                  back_ground_depth=back_ground_depth,back_ground_color_map= back_ground_color_map)

    n = len(transect)
    cmap = matplotlib.cm.get_cmap('Spectral')
    colors = cmap(np.linspace(0,1,n))

    for ii,section in enumerate(transect):
        polygon = [{'points': section['polygon_vertices']}]
        plane = [{'points': section['plane'].vertices}]

        plot_utilities.draw_polygon_list(polygon,ax=ax,color=colors[ii],label=str(ii))
        plot_utilities.draw_polygon_list(polygon,ax=ax,color='black')
        plot_utilities.draw_polygon_list(plane,ax=ax,color='black')
        
        A = section['plane'].vertices[0]
        B = section['plane'].vertices[1]
        ax.quiver(A[0],A[1],B[0]-A[0],B[1]-A[1],angles='xy',zorder=10)
        plt.scatter(A[0],A[1],color=colors[ii])
        plt.scatter(B[0],B[1],color=colors[ii])



    # if track_data is not None:
    x = track_data['x'][nt, :, :2].copy() # copy so as not to change original data
    sel = track_data['status'][nt, :] < min_status # get rid of dead particles
    x[sel,:] = np.nan

    ax.scatter(x[:,0],x[:,1],s=10)

    for ii,section in enumerate(transect):
        particles_in_poly = section['polygon'].is_inside(x)
    
        # slice track_data down to those inside poly and save their indicies
        sub_track_data = x[particles_in_poly]
        plt.scatter(sub_track_data[:,0],sub_track_data[:,1],color=colors[ii])

    if plot_file_name is not None:
        plt.savefig(plot_file_name,dpi=300)

    ax.legend()


def plot_projected_horizontal_tracks(transect,nt=0,plot_file_name=None,axis_lims=None):


    plt.figure(figsize=(20,5))


    x = transect.track_data['x']
    plt.scatter(x[nt,:,0],x[nt,:,1])
    plt.xlim((0,transect.transect[-1]['length'] + transect.transect[-1]['distance_downstream']))
    plt.xlabel('distance downstream [m]')
    plt.ylabel('distance off-center [m]')
    plt.title(f'horizontal tracks at timestep {nt}')
    
    if plot_file_name is not None:
        plt.savefig(plot_file_name,dpi=300)


def plot_projected_verticle_tracks(transect,nt=0,plot_file_name=None,axis_lims=None):
    
    x = transect.track_data['x']
   
    particles_in_transect = ~np.isnan(x[nt,:,0])

    x = x[nt,particles_in_transect] #downstream
    z = transect.track_data['z'][nt,particles_in_transect]
    depth = transect.track_data['water_depth'][nt,particles_in_transect]

    plt.figure()
    plt.scatter(x[:,0],x[:,2])
    plt.scatter(x[:,0],-depth,marker='_')
    plt.xlim((0,transect.transect[-1]['length'] + transect.transect[-1]['distance_downstream']))

    plt.title(f'verticle tracks at timestep {nt}')
    plt.xlabel('distance downstream [m]')
    plt.ylabel('depth [m]')
    
    if axis_lims is not None:
        plt.xlim(axis_lims[:2])
        plt.xlim(axis_lims[2:])

    if plot_file_name is not None:
        plt.savefig(plot_file_name,dpi=300)