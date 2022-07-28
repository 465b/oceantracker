from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 

import os
import numpy as np

path_to_dir = '/scratch/local1/output/22_07_19_test_integration_methods_v01'
cases = load_output_files.get_case_info_files_from_dir(path_to_dir)

for item in cases:
    case_info = load_output_files.read_case_info_file(item)
    name = case_info['output_files']['output_file_base'][-4:]
    RK_order = case_info['full_params']['solver']['RK_order']
    n_sub_steps = case_info['full_params']['solver']['n_sub_steps']
    dispersion = int(case_info['full_params']['dispersion']['A_H'])

    print(name,RK_order,str(n_sub_steps).zfill(3),dispersion)
cases