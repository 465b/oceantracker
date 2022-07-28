# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from oceantracker.user_post_processing import loadOutputFiles
from oceantracker.user_post_processing import statsPlot
from oceantracker.user_post_processing import particlePlot
import seaborn as sns 
sns.set_palette('husl',12)

# %%
name = '22_03_29_resuspsension_test_v01'

root = '/scratch/local1/output/'
file = '_runInfo.json'

path = os.path.join(root,name,'') + name + file
runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path)
# %%
df = statsPlot.load_multicase_msb_stats(path) 
# %%
statsPlot.plot_sa_total_polycount(df,[4,4])
# %%
statsPlot.animate_cases(runCaseInfo)
# %%
particlePlot.plot_relative_height(runCaseInfo,plot_file_name='depth_test.png')   
# %%
for ii in np.arange(100,120):
    particlePlot.plot_path_in_vertical_section(runCaseInfo,plot_file_name='depth_test.png',ncase=22,particleID=ii)   
