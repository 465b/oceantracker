from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 
import matplotlib.pyplot as plt

import os
import numpy as np

def draw_average_first_polygone_occurance(cases):
    
    df = stats_plot.load_multicase_msb_stats(cases) 
   
    # if there are replicates
    if len(list(df)[0]) != 4:
        n_cases = int(max([key[1:4] for key in df]))
        n_replicates = int(max([key[-2:] for key in df]))

        cases_subsets = []
        for ii in range(n_cases):
            cases_subsets.append([key for key in df if int(key[1:4])==ii+1])

    else:
        cases_subsets = [df]

    active_particles_in_polys = {}
    for ii in range(len(cases_subsets)):
        subset = dict([(key,df[key]) for key in cases_subsets[ii]])


        ii_case = cases_subsets[ii][0][:4]
        active_particles_in_polys[ii_case] = {}
        active_particles_in_polys[ii_case]['counts'] = np.zeros((n_replicates,)+df[list(subset)[0]]['m'].shape)

        for jj,key in enumerate(subset):
            active_particles_in_polys[ii_case]['counts'][jj] = subset[key]['m'] + subset[key]['b'] + subset[key]['s']
            active_particles_in_polys[ii_case]['time'] = subset[key]['time']
            active_particles_in_polys[ii_case]['n_sub_step'] = subset[key]['case_properties']['solver']['n_sub_steps']

    time_of_arrival = {}
    for key in active_particles_in_polys:
        buffer = np.zeros((n_replicates))
        for ii in range(n_replicates):
            # idx_arrival = np.where(active_particles_in_polys[key]['counts'][ii,:,0,9] != 0)[0][0]
            idx_arrival = np.where(np.sum(active_particles_in_polys[key]['counts'][ii,:,0,8:],axis=1) > 250)[0][0]
            buffer[ii] = active_particles_in_polys[key]['time'][idx_arrival]
        time_of_arrival[key] = buffer #- min(active_particles_in_polys[key]['time'])

        
    avg_time_of_arrival = [np.average(time_of_arrival[item]) for item in  time_of_arrival]
    std_time_of_arrival = [np.std(time_of_arrival[item]) for item in  time_of_arrival]
    n_sub_steps = [active_particles_in_polys[item]['n_sub_step'] for item in active_particles_in_polys]


    plt.errorbar(n_sub_steps[:int(len(n_sub_steps)/2)],avg_time_of_arrival[:int(len(n_sub_steps)/2)],yerr=std_time_of_arrival[:int(len(n_sub_steps)/2)],label='RK2')
    plt.errorbar(n_sub_steps[int(len(n_sub_steps)/2):],avg_time_of_arrival[int(len(n_sub_steps)/2):],yerr=std_time_of_arrival[int(len(n_sub_steps)/2):],label='RK4')
    plt.legend()
    plt.savefig('integration_scheme_analysis.png')
    print('saved fig')



    

path_to_dir = '/scratch/local1/output/22_07_19_test_integration_methods_v01'
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)
draw_average_first_polygone_occurance(cases)