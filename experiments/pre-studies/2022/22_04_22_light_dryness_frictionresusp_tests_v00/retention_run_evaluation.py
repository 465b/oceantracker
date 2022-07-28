# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
import matplotlib
import numpy as np

import matplotlib.pyplot as plt

import oceantracker.user_post_processing.loadOutputFiles as loadOutputFiles
from oceantracker.user_post_processing import particlePlot
from oceantracker.user_post_processing import loadOutputFiles

import seaborn as sns 
sns.set_palette('husl',12)


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

        #polygon_names = []
        #for item in m['info']['polygon_list']:
        #    polygon_names.append(item['__comment'])

        polygon_names = ['geesthacht',
                         'neuengamme',
                         'kirchwerder',
                         'wilhelmsburg',
                         'harbor',
                         'schulau',
                         'stade',
                         'glückstadt',
                         'freiburg',
                         'brunsbüttel',
                         'ottendorf',
                         'cuxhafen',
                         'north_sea']

        
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
            ratio = case['case_properties']['user_trajectory_modifiers'][1]['probability_of_splitting']
            total = case['m'] + case['s'] + case['b']
            for ii in np.arange(poly_range[0],poly_range[1]):
                ax[kk][jj].plot(case['time'].astype('datetime64[s]'),total[:,0,ii],label=case['polygon_names'][ii])
            ax[kk][jj].plot(case['time'].astype('datetime64[s]'),np.sum(total[:,0,:],axis=1),label='total')
            ax[kk][jj].set_title(f'run: {case["name"]} velo: {velo} ratio: {ratio}')

    #for label in ax.xaxis.get_ticklabels():
    #        label.set_rotation(45)
    handles, labels = ax[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', ncol=3, bbox_to_anchor=(.75, 0.98))

    plt.tight_layout()
    plt.savefig('sa_output.svg')
    #plt.show()

# %%
def animate_cases(runCaseInfo):
        
    for ii,item in enumerate(runCaseInfo['output_files']['caseInfo']):

        case_name = item.split('.')[0][:-9]
        velo = runCaseInfo['caseInfo'][ii]['case_params']['user_velocity_modifiers'][0]['mean']
        try:
            split = runCaseInfo['caseInfo'][ii]['case_params']['user_trajectory_modifiers'][1]['fraction_to_split']
        except:
            split = 0
        title = f'run: {case_name[-6:]} velo: {velo} ratio: {split}'
        particlePlot.animate_particles(
            runCaseInfo,movie_file=case_name+'.mp4',fps=30,ncase=ii,title=title,dpi=300)

# %%
path = '/scratch/local1/output/22_04_22_light_dryness_frictionresusp_tests_v00/22_04_22_light_dryness_frictionresusp_tests_v00_runInfo.json'


runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path)
# %%
df = load_multicase_msb_stats(path) 

# %%
%matplotlib widget
plot_sa_total_polycount(df,[3,3])

#%%time
animate_cases(runCaseInfo)




# %%
#particlePlot.plot_relative_height(runCaseInfo,plot_file_name='depth_test.png')   


# %%
#for ii in np.arange(100,120):
#    particlePlot.plot_path_in_vertical_section(runCaseInfo,plot_file_name='depth_test.png',ncase=22,particleID=ii)   

# %%
