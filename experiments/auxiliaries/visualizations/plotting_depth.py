# %%
from oceantracker.user_post_processing import loadOutputFiles
from oceantracker.user_post_processing.particlePlot import draw_grid
path = '/scratch/local1/output/22_04_11_light_dryness_tests_v01/22_04_11_light_dryness_tests_v01_runInfo.json'
runCaseInfo =  loadOutputFiles.load_Run_and_CaseInfo(path, ncase = 0)
grid = loadOutputFiles.load_grid(runCaseInfo)


# %%
from netCDF4 import Dataset
df = Dataset('/scratch/local1/hzg2/schout_1.nc')
# %%
topo = df.variables['depth']
water_level = df.variables['elev']
water_depth = water_level[0] - topo
water_depth[water_depth > 0] = 0
water_depth[water_depth < -20] = -20

#%%
import matplotlib.pyplot as plt
fig,ax = plt.subplots(figsize=(15,10))
#draw_grid(grid)
tc = ax.tripcolor(grid['x'][:,0,],grid['x'][:,1],water_depth,cmap='magma')
plt.colorbar(tc, ax=ax)
plt.xlim(560e3,580e3)
plt.ylim(5.92e6,5.94e6)
plt.show()

# %%
%matplotlib widget
fig,ax = plt.subplots(figsize=(15,10))
#tc = ax.tripcolor(grid['x'][:,0,],grid['x'][:,1],water_depth,triangles=grid['triangles'],cmap='plasma')
tc = ax.tripcolor(grid['x'][:,0,],grid['x'][:,1],water_depth,cmap='plasma')
plt.colorbar(tc, ax=ax)
plt.show()

