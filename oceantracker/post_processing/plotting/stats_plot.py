import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from oceantracker.post_processing.read_output_files import load_output_files
from oceantracker.post_processing.plotting import plot_tracks
# %%
def load_multicase_msb_stats(cases,replicate_set=0):

    # load data from cases in dir
    
    num_cases = len(cases)
    case_info = load_output_files.read_case_info_file(cases[0]) 
    replicates = case_info['shared_params']['replicates']
    
    polygon_names = [item['user_polygon_name'] for item in case_info['full_params']['particle_statistics'][0]['polygon_list']]
   
    case_stats = {}
    for case in cases:
        moving =   load_output_files.load_stats_file(case,nsequence=1)
        stranded = load_output_files.load_stats_file(case,nsequence=2)
        bottom =   load_output_files.load_stats_file(case,nsequence=3)

        if replicates == 1:
            run_name = moving['info']['output_file'][-30:-26]
        else:
            run_name = moving['info']['output_file'][-33:-26]


        singlecase_df = {
            'name': run_name,
            'time': moving['time'],
            'm': moving['count'],
            's': stranded['count'],
            'b': bottom['count'],
            'polygon_names': polygon_names,
            'case_properties': moving['full_params']}
        
        
        case_stats[run_name] = singlecase_df

    return case_stats



# %%
def plot_sa_total_polycount(cases,title='',poly_range=(0,13),savefig=True,output_path=''):
    # todo:

    df = load_multicase_msb_stats(cases) 

    # if there are replicates
    if len(list(df)[0]) != 4:
        n_replicates = int(max([key[-2:] for key in df]))

        # split replicates into subsets
        replicate_subsets = []
        for ii in range(n_replicates):
            replicate_subsets.append([key for key in df if int(key[-2:])==ii+1])

    # no replicates case
    else:
        replicate_subsets = [df]

    for ii in range(len(replicate_subsets)):
        subset = dict([(key,df[key]) for key in replicate_subsets[ii]])

        # prepare data to slice cases by migration methods
        try:
            # try and check if there are velocity modifiers present
            vel_mods = dict([(key, subset[key]['case_properties']['velocity_modifiers']) for key in subset])
            v_vel = dict([(key, vel_mods[key][0]['mean']) for key in subset])
            
            # get case keys for each migration method
            drifting = [key for key in vel_mods if 
                (vel_mods[key][0]['class_name']=='oceantracker.particle_velocity.terminal_velocity.AddTerminalVelocity') 
                and (vel_mods[key][0]['mean']==0)]
            monotonic  = [key for key in vel_mods if 
                (vel_mods[key][0]['class_name']=='oceantracker.particle_velocity.terminal_velocity.AddTerminalVelocity')
                and (vel_mods[key][0]['mean']!=0)]
            diel  = [key for key in vel_mods if
                (vel_mods[key][0]['class_name']=='oceantracker.particle_velocity.terminal_velocity.AddDielVelocity') 
                and (vel_mods[key][0]['mean']!=0)]
            
        except IndexError:
            # if not just set them all to zero and add them to drifting
            v_vel = dict([(key, 0) for key in subset])
            drifting = [key for key in vel_mods]
            monotonic = []
            diel = []
            
        try:
            # if splitting is present
            splitting = dict([(key, subset[key]['case_properties']['trajectory_modifiers'][0]['probability_of_splitting']) for key in subset])
        except:
            # if not set it to zero
            splitting = dict([(key, 0) for key in subset])

        
        # plot each migration pattern seperatly
        if len(drifting) != 0:
            _draw_sa_figure(subset, ii, drifting, 'drifting', splitting, v_vel, poly_range, output_path)
        if len(monotonic) != 0:
            _draw_sa_figure(subset, ii, monotonic, 'monotonic', splitting, v_vel, poly_range, output_path)
        if len(diel) != 0:
            _draw_sa_figure(subset, ii, diel, 'diel', splitting, v_vel, poly_range, output_path)



def _draw_sa_figure(rep_subset, ii_rep_subset, migration_pattern, migration_name, splitting, v_vel, poly_range, fig_path):
    cases = [rep_subset[key] for key in migration_pattern]
    
    # then check for lenght of dimensions
    v_vel_of_migration_pattern = [v_vel[key] for key in migration_pattern]
    splitting_of_migration_pattern = [splitting[key] for key in migration_pattern]
    
    # reshape it accordingly
    # plot_shape = (max(1,len(set(v_vel_of_migration_pattern))),max(1,len(set(splitting_of_migration_pattern))))
    plot_shape = (6,4)
    
    # draw figure
    fig,ax = plt.subplots(plot_shape[0],plot_shape[1],sharex=True,figsize=(24,12))
    if plot_shape[0] == 1: ax = [ax]
    if plot_shape[1] == 1: ax = [[item] for item in ax]

    # draw individual subfigure frames for each case
    for kk in range(plot_shape[0]):
        for jj in range(plot_shape[1]):
            
            case = cases[plot_shape[1]*kk+jj]
            total = case['m'] + case['s'] + case['b']
            
            # draws the indivual polygon stats
            for ll in np.arange(poly_range[0],poly_range[1]):
                ax[kk][jj].plot(case['time'].astype('datetime64[s]'),
                    total[:,0,ll],label=case['polygon_names'][ll])
            
            # draws total particles
            ax[kk][jj].plot(case['time'].astype('datetime64[s]'),
                np.sum(total[:,0,:],axis=1),label='total')

            # set title  
            velo = v_vel[case['name']]
            ratio = splitting[case['name']]
            ax[kk][jj].set_title(f'run: {case["name"]} velo: {round(velo,3)} ratio: {round(ratio,5)}')

            #for label in ax.xaxis.get_ticklabels():
            #        label.set_rotation(45)
            ax[kk][jj].xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    
    handles, labels = ax[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', ncol=1, bbox_to_anchor=(1.0001, 0.9))
    # plt.tight_layout()

    if fig_path is not None: 
        filename = os.path.join(fig_path,os.path.split(fig_path)[-1])
        filename = f'{filename}_{migration_name}_R{ii_rep_subset}'
        plt.savefig(filename)
        print(f'Saved overview figure (polygon counts) to {filename}')

# %%
def animate_cases(cases):
    
    for case in cases:
        case_info = load_output_files.read_case_info_file(case)
        track = load_output_files.load_particle_track_vars(case)
        # run_name = case_info['shared_params']['output_file_base']
        case_name_long = case_info['output_files']['output_file_base']
        case_name_short = case_name_long.split('_')[-1]
        output_dir = case_info['output_files']['run_output_dir']
        
        try:
            v_vel = case_info['full_params']['velocity_modifiers'][0]['mean']
        except: 
            v_vel = 0
            
        try: 
            split = case_info['full_params']['trajectory_modifiers'][0]['probability_of_splitting']
        except:
            split = 0

        filename = f'{os.path.join(output_dir,case_name_long)}.mp4'
        plot_tracks.animate_particles(track,
            movie_file=filename,fps=30,dpi=300,
            title=f'{case_name_short} vel:{round(v_vel,3)} split:{round(split,6)}')
        print(f'Saved particle animation to {filename}')
