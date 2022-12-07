#%%
from oceantracker.post_processing.read_output_files import load_output_files 
from oceantracker.post_processing.plotting import stats_plot 

import os
import numpy as np

import importlib
importlib.reload(stats_plot)


#%% load and draw statistical overview
path_to_dir = '/scratch/local1/output/22_11_28_mini_test_set_w_rep'

retention_mini = stats_plot.retention_data(path_to_dir)

# retention_mini.plot_retention_sa_polycount_overview(fig_path='/scratch/local1/output/22_11_28_mini_test_set_w_rep')
retention_mini.plot_retention_sa_sucess_overview()
retention_mini.plot_retention_box_plots()