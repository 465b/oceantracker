# %%
from oceantracker.user_post_processing import loadOutputFiles

# %%
path = '/scratch/local1/output/22_04_11_light_dryness_tests_v01/22_04_11_light_dryness_tests_v01_runInfo.json'                                                                                                                  

# %%
ncase=0 
runCaseInfo = loadOutputFiles.load_Run_and_CaseInfo(path) 
tracks = loadOutputFiles.load_particle_track_vars(runCaseInfo, ['x', 'status', 'time','light_limitation'], ncase=ncase)


# %%
