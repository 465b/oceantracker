import multiprocessing
import os
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from oceantracker.post_processing.plotting import plot_tracks
from oceantracker.post_processing.read_output_files import load_output_files
from oceantracker.util import json_util

sns.set_context("talk")


# %%

class retention_data:

    def __init__(self,run_dir_path,drop_dublicates=True,compress_verticle_modes=True):
        # keep it try'y
        run_case_info = load_output_files.get_run_info_files_from_dir(run_dir_path)
        
        metadata = {}

        # read runInfo to get run-meta-info
        metadata["n_rep"] = run_case_info['performance']['replicates']
        
        metadata["n_indiv_cases"] = len(run_case_info['user_params']['case_list'])
        metadata["n_total_cases"] = metadata["n_indiv_cases"] * metadata["n_rep"]
       
        metadata["stat_poly_names"] = \
            [item['user_polygon_name'] for item in run_case_info['user_params']['base_case_params']['particle_statistics'][0]['polygon_list']]

        coulmn_header_type = [
            'case_info','case_info','case_info','case_info',
            'statistical_polygon_data','statistical_polygon_data',
            'statistical_polygon_data','statistical_polygon_data',
            'statistical_polygon_data',
            'track_derived_data','track_derived_data',
            'track_derived_data','track_derived_data',
            'track_derived_data','track_derived_data'
            ]
        column_header = [
            'name','vert_mode',"vert_vel","split_frac",
            'model_time','n_particle_moving','n_particle_stranded','n_particle_bottom','n_particle_total',
            'depth_below_free_surface_long_living', 'depth_below_free_surface_short_living',
            'distance_traveled_long_living', 'distance_traveled_short_living',
            'ratio_stranded_long_living', 'ratio_stranded_short_living']
        
        data = pd.DataFrame(
            columns=column_header,
            # columns=pd.MultiIndex.from_tuples(list(zip(coulmn_header_type,column_header))),
            index=range(metadata["n_total_cases"])
            )

        case_infos = load_output_files.get_case_info_files_from_dir(run_dir_path)
       
        for ii,case in enumerate(case_infos):
            print(f"Current file: {ii}/{len(case_infos)}\nCurrent size of data: {round(sys.getsizeof(data)/1e9,2)}GB", end="\r")
            case_info = load_output_files.read_case_info_file(case)

            # name of the current case
            data.iloc[ii]["name"] = case_info['output_files']['output_file_base'].split('_')[-1]

            try:
                data.iloc[ii]["vert_vel"] = case_info['full_params']['velocity_modifiers'][0]['mean']
            except KeyError:
                data.iloc[ii]["vert_vel"] = 0
            
            try:
                data.iloc[ii]["split_frac"] =  case_info['full_params']['trajectory_modifiers'][0]['probability_of_splitting']
            except KeyError:
                data.iloc[ii]["split_frac"] = 0


            try:
                data.iloc[ii]["vert_mode"] = case_info['full_params']['velocity_modifiers'][0]['class_name']
            except KeyError:
                data.iloc[ii]["vert_mode"] = ['oceantracker.velocity_modifiers.terminal_velocity.TerminalVelocity']

            diel = 'oceantracker.velocity_modifiers.terminal_velocity.DielVelocity'
            mono = 'oceantracker.velocity_modifiers.terminal_velocity.TerminalVelocity'

            if (data.iloc[ii]["vert_mode"] == diel) and (data.iloc[ii]["vert_vel"] != 0):
                data.iloc[ii]["vert_mode"] = 'diel'
            elif (data.iloc[ii]["vert_mode"] == mono) and (data.iloc[ii]["vert_vel"] != 0):
                data.iloc[ii]["vert_mode"] = 'mono'
            elif (data.iloc[ii]["vert_mode"] == mono) and (data.iloc[ii]["vert_vel"] == 0):
                data.iloc[ii]["vert_mode"] = 'drift'

            try:
                df = self.load_stat_poly_info(case)
                
                data.iloc[ii]["model_time"] = df['time']
                data.iloc[ii]["n_particle_moving"] = df['m']
                data.iloc[ii]["n_particle_stranded"] = df['s']
                data.iloc[ii]["n_particle_bottom"] = df['b']
                data.iloc[ii]["n_particle_total"] = df['m']+df['s']+df['b']
            
            except KeyError:
                data.iloc[ii]["model_time"] = np.nan
                data.iloc[ii]["n_particle_moving"] = np.nan
                data.iloc[ii]["n_particle_stranded"] = np.nan
                data.iloc[ii]["n_particle_bottom"] = np.nan
                data.iloc[ii]["n_particle_total"] = np.nan
            
            try:
                df = self.load_track_derived_metrics(case)
                
                data.iloc[ii]["depth_below_free_surface_long_living"] = df['bfs_long']
                data.iloc[ii]["depth_below_free_surface_short_living"] = df['bfs_short']
                data.iloc[ii]["distance_traveled_long_living"] = df['dt_long']
                data.iloc[ii]["distance_traveled_short_living"] = df['dt_short']
                data.iloc[ii]["ratio_stranded_long_living"] = df['str_long_ratio']
                data.iloc[ii]["ratio_stranded_short_living"] = df['str_short_ratio']

            except KeyError:
                data.iloc[ii]["depth_below_free_surface_long_living"] = np.nan
                data.iloc[ii]["depth_below_free_surface_short_living"] = np.nan
                data.iloc[ii]["distance_traveled_long_living"] = np.nan
                data.iloc[ii]["distance_traveled_short_living"] = np.nan
                data.iloc[ii]["ratio_stranded_long_living"] = np.nan
                data.iloc[ii]["ratio_stranded_short_iving"] = np.nan


        # reshaping the dataframe into cases and repetitions 
        multi_index = np.empty((2,metadata["n_total_cases"]),dtype=int)

        if metadata["n_rep"] == 1:
            case = [int(item.split('C')[1]) for item in data["name"]]
            multi_index[0,:] = case
            multi_index[1,:] = 1

        else:            
            for ii in range(metadata["n_rep"]):
                current_rep = 'R'+str(ii+1).zfill(2) 
                rep_subset = np.where([current_rep in item for item in data["name"]])[0]
                case = [int(item.split('R')[0][1:]) for item in data["name"][rep_subset]]
                
                multi_index[0,rep_subset] = case
                multi_index[1,rep_subset] = ii+1

        index = pd.MultiIndex.from_tuples(list(zip(*multi_index)),names=['case','rep'])
        data = data.set_index(index)

        # removing unneccessary repreats
        cond = (np.array(data["vert_mode"] != 'mono') * np.array(data["vert_mode"] != 'diel') * np.array(data["vert_mode"] != 'drift'))
        if cond.any():
            print('Warning: unknown verticle modes found or unneccessary drift repeats found')
            data = data[~cond]
            metadata['n_total_cases'] = int(data.shape[0])
            metadata['n_indiv_cases'] = int(data.shape[0]/metadata['n_rep'])
        
        metadata["vert_mode"] = sorted(list(set(data["vert_mode"])))
        metadata["vert_vel"] = sorted(list(set(data["vert_vel"])))
        metadata["vert_vel"].remove(0)
        metadata['split_frac'] = sorted(list(set(data['split_frac'])),reverse=True)

        if drop_dublicates == True:
            # tmp add rep to columns for dublication selection 
            # (otherwise the reps are labeled as well)

            data['rep'] = data.index.get_level_values(1)
            data = data.drop_duplicates(subset=['rep','vert_mode','vert_vel','split_frac'])
            data = data.drop(labels=["rep"],axis=1)


        if compress_verticle_modes == True:
            drift = data[data["vert_mode"] == "drift"]

            data = data.drop(drift.index)

            for vert_mode in set(data["vert_mode"]):
                tmp = drift.copy()
                tmp["vert_mode"] = vert_mode
                data = pd.concat([data,tmp])

            data = data.sort_values(by=["vert_mode","vert_vel","split_frac"],
                                    ascending=[True,True,False])

            metadata["n_total_cases"] = len(data)
            metadata["n_indiv_cases"] = int(len(data)/metadata["n_rep"])
            metadata["vert_mode"] = sorted(list(set(data["vert_mode"])))
            metadata["vert_vel"] = sorted(set(data["vert_vel"]))

        self.metadata = metadata
        self.data = data

            
    def plot_retention_sa_polycount_overview(self,title='',poly_range=(0,13),mode=None,fig_path=None):

        for ii in range(self.metadata['n_rep']):
            for mode in set(self.metadata['vert_mode']):
                self._draw_retention_polycounts_sa_figure(ii_replicate=ii+1, migration_mode=mode, fig_path=fig_path)


    def plot_retention_sa_sucess_overview(self,title='',poly_range=(0,13),mode=None,fig_path=None):

        for ii in range(self.metadata['n_rep']):
            for mode in set(self.metadata['vert_mode']):
                self._draw_retention_success_sa_figure(ii_replicate=ii+1, migration_mode=mode, fig_path=fig_path)
    
    
    def plot_retention_box_plots(self,fig_path=None, average_repetitions=True):

        if average_repetitions:

            fig,ax = plt.subplots()
            ax = self._draw_avg_dbf_box_plot(ax)

            if fig_path is not None:
                filename = filename = os.path.split(fig_path)[-1]
                filename = os.path.join(fig_path,filename+'_dbf_avg.png')
                plt.savefig(filename)
                print(f'Saved depth below free surface averaged figure to {filename}')
            
            fig,ax = plt.subplots()
            ax = self._draw_avg_dt_box_plot(ax)

            if fig_path is not None:
                filename = filename = os.path.split(fig_path)[-1]
                filename = os.path.join(fig_path,filename+'_dt_avg.png')
                plt.savefig(filename)
                print(f'Saved distance traveled averaged figure to {filename}')

            fig,ax = plt.subplots()
            ax = self._draw_avg_rs_box_plot(ax)

            if fig_path is not None:
                filename = filename = os.path.split(fig_path)[-1]
                filename = os.path.join(fig_path,filename+'_rs_avg.png')
                plt.savefig(filename)
                print(f'Saved distance traveled averaged figure to {filename}')

    def _draw_retention_box_plot(self, ii_replicate=1, migration_mode='drift', fig_path=None):
        
        df = self.data.xs(ii_replicate,level=1)
        df = df[df["vert_mode"]==migration_mode]
        
        if migration_mode == 'drift':
            plot_shape = (1 ,max(1,len(self.metadata['split_frac'])))
        else:
            plot_shape = (max(1,len(set(self.metadata['vert_vel']))) ,max(1,len(self.metadata['split_frac'])))


    def _draw_avg_dbf_box_plot(self,ax):

        short = np.concatenate([item for item in self.data['depth_below_free_surface_short_living']])
        long = np.concatenate([item for item in self.data['depth_below_free_surface_long_living']])

        ax.set_title('Plankton depth below surface')
        ax.boxplot(short,positions=[1],labels=['short\nliving'],
                   #whis=(0,100)
                   showfliers=False
                   )
        ax.boxplot(long,positions=[2],labels=['long\nliving'],
                   #whis=(0,100)
                   showfliers=False
                   )
        ax.set_ylabel('depth (m)')

        return ax


    def _draw_avg_dt_box_plot(self,ax):

        short = np.concatenate([item for item in self.data['distance_traveled_short_living']])
        long = np.concatenate([item for item in self.data['distance_traveled_long_living']])

        ax.set_title('Plankton distance traveled')
        # ax.set_yscale('log')
        ax.boxplot(short,positions=[1],labels=['short\nliving'],
                   #whis=(0,100),
                   showfliers=False
                   )
        ax.boxplot(long,positions=[2],labels=['long\nliving'],
                   #whis=(0,100),
                   showfliers=False
                   )
        ax.set_ylabel('depth (m)')
        # try:
        #     ax.set_ylim(0.9*min(np.min(short),np.min(long))),1.1*max(np.max(short),np.max(long))
        # except ValueError:
        #     pass

        return ax


    def _draw_avg_rs_box_plot(self,ax):

        short = np.concatenate([[item] for item in self.data['ratio_stranded_short_living']])
        long = np.concatenate([[item] for item in self.data['ratio_stranded_long_living']])

        ax.set_title('Ratio of plankton stranded')
        ax.boxplot(short,positions=[1],labels=['short\nliving'],whis=(0,100))
        ax.boxplot(long,positions=[2],labels=['long\nliving'],whis=(0,100))
        ax.set_ylabel('ratio of plankton stranded')

        return ax


    


    def _draw_retention_success_sa_figure(self, ii_replicate=1, migration_mode='drift', poly_range=(0,13), fig_path=None):

        # slicing for 
        df = self.data.xs(ii_replicate,level=1)
        df = df[df["vert_mode"]==migration_mode]
        
        if migration_mode == 'drift':
            plot_shape = (1 ,max(1,len(self.metadata['split_frac'])))
        else:
            plot_shape = (max(1,len(self.metadata['vert_vel'])) ,
                          max(1,len(self.metadata['split_frac'])))

        
        # filter particle counts for the surviving particles
        total = [np.sum(case[:,0,:],axis=1) for case in df['n_particle_total']]
        threshold = np.average([case[0] for case in total])

        surviving = [case[-1] for case in total]
        if len(surviving) > plot_shape[0]*plot_shape[1]:
            print("Warning. Dataset surpisingly large. Cutting down to expected size")
            surviving = surviving[:plot_shape[0]*plot_shape[1]]
        surviving = np.reshape(surviving, plot_shape)

        # draw figure
        fig,ax = plt.subplots(1,1,figsize=(24,12))

        plt.title(f'Retention success - migration mode: {migration_mode}, replicate: {ii_replicate}')

        if migration_mode == 'drift':
            sns.heatmap(surviving, 
                cmap='RdYlGn',
                center=threshold,
                xticklabels=['{:,.1e}'.format(item) for item in  self.metadata['split_frac']],
                yticklabels=[0])
        else:
            sns.heatmap(surviving,
                cmap='RdYlGn',
                center=threshold,
                xticklabels=['{:,.1e}'.format(item) for item in  self.metadata['split_frac']],
                yticklabels=[round(item,4) for item in self.metadata["vert_vel"]]
                )

        plt.xlabel('growth rate as a splitting fraction')
        plt.ylabel('vertical velocity (m/s)')

        
        if fig_path is not None: 
            filename = os.path.join(fig_path,os.path.split(fig_path)[-1])
            filename = f'{filename}_retention_success_{migration_mode}_R{ii_replicate}'
            plt.savefig(filename)
            print(f'Saved overview figure (retenion success) to {filename}')


    def _draw_retention_polycounts_sa_figure(self, ii_replicate=1, migration_mode='drift', poly_range=(0,13), fig_path=None):
        

        # slicing for 
        df = self.data.xs(ii_replicate,level=1)
        df = df[df["vert_mode"]==migration_mode]
        
        if migration_mode == 'drift':
            plot_shape = (1 ,max(1,len(self.metadata['split_frac'])))
        else:
            plot_shape = (max(1,len(self.metadata['vert_vel'])), max(1,len(self.metadata['split_frac'])))

        # draw figure
        fig,ax = plt.subplots(plot_shape[0],plot_shape[1],sharex=True,sharey=True,figsize=(24,12))

        fig.suptitle('Polycounts - mode: '+migration_mode+' replicate: '+str(ii_replicate))

        if plot_shape[0] == 1: ax = [ax]
        if plot_shape[1] == 1: ax = [[item] for item in ax]

        # draw individual subfigure frames for each case
        for kk in range(plot_shape[0]):
            for jj in range(plot_shape[1]):
                ii = plot_shape[1]*kk+jj                

                # might need to sort first or select not for index but for entry in col
                case = df.iloc[ii]
                ax[kk][jj] = self.plot_indiv_polycount(ax[kk][jj], case) 
        
        plt.setp(ax[-1][:], xlabel='Particle Counts')
        [plt.setp(axis[0],  ylabel='Time') for axis in ax]

        handles, labels = ax[0][0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper right', ncol=1, bbox_to_anchor=(1.0001, 0.9))
        plt.tight_layout()

        if fig_path is not None: 
            filename = os.path.join(fig_path,os.path.split(fig_path)[-1])
            filename = f'{filename}_{migration_mode}_R{ii_replicate}'
            plt.savefig(filename)
            print(f'Saved overview figure (polygon counts) to {filename}')


    def plot_indiv_polycount(self,ax,case,poly_range=(0,13),labels=False,legend=False,savefig=False,output_path=''):
            
        total = case['n_particle_total']
        title =  f"run: {case['name']} velo: {round(case['vert_vel'],3)} ratio: {round(case['split_frac'],5)}"

        for ll in np.arange(poly_range[0],poly_range[1]):
            ax.plot(case['model_time'].astype('datetime64[s]'),
                total[:,0,ll],label=self.metadata['stat_poly_names'][ll])
        
        # draws total particles
        ax.plot(case['model_time'].astype('datetime64[s]'),
            np.sum(total[:,0,:],axis=1),label='total')

        try:
            # set title  
            velo = v_vel[case['name']]
            ratio = splitting[case['name']]
            ax.set_title(f'run: {case["name"]} velo: {round(velo,3)} ratio: {round(ratio,5)}')
        except NameError:
            print('Warning: No case title set')


        if labels == True:
            ax.set_xlabel('Time')
            ax.set_ylabel('Particle Count')
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))

        if legend == True:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels)#, loc='upper right', ncol=1, bbox_to_anchor=(1.0001, 0.9))

        if savefig==True:
            filename = case['name']
            filename = os.path.join(output_path,filename+'.png')
            plt.savefig(filename)
            print(f'Saved individual polycount figure to {filename}')

        return ax





    @staticmethod
    def load_stat_poly_info(case,replicate_set=0):

        moving =   load_output_files.load_stats_file(case,nsequence=1)
        stranded = load_output_files.load_stats_file(case,nsequence=2)
        bottom =   load_output_files.load_stats_file(case,nsequence=3)

        df = {
            'time': moving['time'],
            'm': moving['count'],
            's': stranded['count'],
            'b': bottom['count']
        }

        return df


    @staticmethod
    def load_track_derived_metrics(case, 
        var_list=['x', 'time','status', 'age', 'tide', 'water_depth','distance_travelled'],
        min_status=0,age_threshold = 3*60*60*24*28):
        
        
        case_info = load_output_files.read_case_info_file(case) 
        tracks = load_output_files.load_particle_track_vars(case,var_list)

        # drop non-living particles
        sel = tracks['status'][:, :] < min_status
        tracks['x'][sel] = np.nan

        # seperate all particles into long- and short-lived
        long_lived =  tracks['age'][-1] > age_threshold
        short_lived = ~long_lived

        metrics = {}

        # 1st statistical metric    - below_free_surface
        if 'tide' in var_list:
            try:
                below_free_surface = tracks['x'][:,:,2] - tracks['tide']
                metrics['bfs_long'] = flatclean(below_free_surface[:,long_lived])
                metrics['bfs_short'] = flatclean(below_free_surface[:,short_lived])
            except KeyError:
                pass
        
        # 2nd                       - distance traveled
        if 'distance_travelled' in var_list:
            try:
                distance_travelled = tracks['distance_travelled']
                metrics['dt_long'] = flatclean(distance_travelled[:,long_lived])
                metrics['dt_short'] = flatclean(distance_travelled[:,short_lived])
            except KeyError:
                pass

        # 3rd                       - stranded
        try:
            status_long = tracks['status'][:,long_lived]
            status_short = tracks['status'][:,short_lived]

            str_long = flatclean(status_long[status_long == 3])
            nstr_long = flatclean(status_long[(status_long != 3)&(status_long > 0)])
            str_short = flatclean(status_short[status_short == 3])
            nstr_short = flatclean(status_short[(status_short != 3)&(status_short > 0)])
            
            if len(str_long) != 0:
                metrics['str_long_ratio'] = len(str_long)/(len(str_long) + len(nstr_long))
            else:
                metrics['str_long_ratio'] = np.nan

            if len(str_short) != 0:
                metrics['str_short_ratio'] = len(str_short)/(len(str_short) + len(nstr_short))
            else:
                metrics['str_short_ratio'] = np.nan

        except KeyError:
            pass


        return metrics


def flatclean(x):
    x = x.flatten()
    x = x[~np.isnan(x)]
    return x







def plot_all_polycounts_individually(cases,title='',poly_range=(0,13),savefig=True,output_path=''):
    df = load_stat_poly_info(cases) 

    for run in df:
        case = df[run]

        fig,ax = plt.subplots(1,1)
        plot_indiv_polycount(ax, case, labels=True, legend=True, savefig=savefig, output_path=output_path)


# %%
def animate_cases(cases):
    with multiprocessing.Pool(processes=16) as pool:
        pool.map(animate1case,cases)


def animate1case(case):
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

    filename = f'{os.path.join(output_dir,case_name_long)}.gif'
    plot_tracks.animate_particles(track,
        movie_file=filename,fps=30,dpi=300,
        title=f'{case_name_short} vel:{round(v_vel,3)} split:{round(split,6)}')
    print(f'Saved particle animation to {filename}')
