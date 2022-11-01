# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%

# %%
import numpy as np

import oceantracker.user_post_processing.loadOutputFiles as loadOutputFiles
from oceantracker.user_post_processing import particlePlot

import matplotlib.pyplot as plt
import seaborn as sns 
sns.set_palette('husl',8)


# %%
def load_multicase_msb_stats(path_runInfo,replicate_set=0):

    runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path)

    num_cases = runCaseInfo['performance']['num_cases']
    replicates = runCaseInfo['performance']['replicates']

    multicase_df = {}
    
    for kk in range(0,int(num_cases/replicates),replicates):
        ii = kk+replicate_set


        #moving
        m = loadOutputFiles.load_stats_file(path_runInfo,ncase=ii,nsequence=0)
        #stranded
        s = loadOutputFiles.load_stats_file(path_runInfo,ncase=ii,nsequence=1)
        #bottom
        b = loadOutputFiles.load_stats_file(path_runInfo,ncase=ii,nsequence=2)
        
        file_name = m['file_name']
        run_name = file_name[-28:-22]
        print(kk/replicates,run_name,file_name)

        polygon_names = []
        for item in m['info']['polygon_list']:
            polygon_names.append(item['__comment'])
        
        singlecase_df = {
            'name': run_name,
            'time': m['time'],
            'm': m['count'],
            's': s['count'],
            'b': b['count'],
            'polygon_names': polygon_names,
            'case_properties': runCaseInfo['caseInfo'][ii]['case_params']}
        
        
        multicase_df[run_name] = singlecase_df

    return multicase_df



# %%
def plot_sa_total_polycount(df,sa_shape,title='',poly_range=(0,12)):


    cases = [df[item] for item in df]
    fig,ax = plt.subplots(sa_shape[0],sa_shape[1],sharex=True,figsize=(12,6))


    if sa_shape[0] == 1: ax = [ax]
    if sa_shape[1] == 1: ax = [[item] for item in ax]

    kk = 0

    for kk in range(sa_shape[0]):
        for jj in range(sa_shape[1]):
            print(kk*sa_shape[1]+jj,kk,jj)
            case = cases[sa_shape[1]*kk+jj]
            velo = case['case_properties']['user_velocity_modifiers'][0]['mean']
            ratio = case['case_properties']['user_trajectory_modifiers'][1]['fraction_to_split']
            total = case['m'] + case['s'] + case['b']
            for ii in np.arange(kk*sa_shape[1]+jj+1):
                ax[kk][jj].plot(case['time'].astype('datetime64[s]'),total[:,0,ii],label=case['polygon_names'][ii])
            #ax[kk][jj].plot(case['time'].astype('datetime64[s]'),np.sum(total[:,0,:],axis=1),label='total',color='black')
            ax[kk][jj].set_title(f'run: {case["name"]} velo: {velo} ratio: {ratio}')

    #for label in ax.xaxis.get_ticklabels():
    #        label.set_rotation(45)
    handles, labels = ax[kk][jj].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', ncol=3, bbox_to_anchor=(.75, 0.98))

    plt.tight_layout()
    plt.savefig('sa_output.svg')
    plt.show()


# %%
path = '/scratch/local1/output/22_02_09_horizontal_migration_experiment_v02/22_02_09_horizontal_migration_experiment_v02_runInfo.json'

# %%
df = load_multicase_msb_stats(path) 

# %%
plot_sa_total_polycount(df,[6,4])

# %%
#%%time
runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path)

for ii,item in enumerate(runCaseInfo['output_files']['caseInfo']):
    particlePlot.animate_particles(runCaseInfo,movie_file=item+'.mp4',fps=10,ncase=ii,title='test',dpi=300)



# %%



