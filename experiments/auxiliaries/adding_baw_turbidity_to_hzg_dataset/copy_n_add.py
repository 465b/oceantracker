# %%
## 

# %%
from netCDF4 import Dataset

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import matplotlib.tri as tri



# %% [markdown]
# ## Merging Turbidity data from BAWs Untrim data set with Johannes' SCHISM grid
# 
# 

# %% [markdown]
# ### BAW data

# %%
baw = Dataset('/scratch/local1/baw/f.AZHel_FT_REF.2D.cut.nc', 'r')

# %%
# _ = [print(item) for item in baw.variables.keys()]

# %%
baw_face_x = baw.variables['Mesh2_face_x'][:]
baw_face_y = baw.variables['Mesh2_face_y'][:]

baw_spm = baw.variables['Mesh2_face_Schwebstoffgehalt_2d']
baw_spm = np.average(baw_spm[:, 0, 0, :], axis=0)

# %%
fig,ax = plt.subplots(figsize=(10,10),dpi=300)
ax.scatter(baw_face_x, baw_face_y, c=baw_spm, cmap='viridis')
ax.set_aspect('equal')

# %% [markdown]
# ### HZG data

# %%
hzg = Dataset('/scratch/local1/output/22_11_01_depth_losses_v10/22_11_01_depth_losses_v10_grid.nc', 'r')

# %%
# _ = [print(item) for item in hzg.variables.keys()]

# %%

hzg_node_x = hzg.variables['x'][:,0]
hzg_node_y = hzg.variables['x'][:,1]

hzg_tri = hzg.variables['triangles'][:]

# %%
fig,ax = plt.subplots(figsize=(10,10),dpi=300)
ax.triplot(hzg_node_x, hzg_node_y, hzg_tri)
ax.set_aspect('equal')

# %% [markdown]
# ### Interpolating and merging BAW data onto HZG grid

# %% [markdown]
# #### delaunay triangulation
# 

# %% [markdown]
# #### griddata interpolation

# %%


# Create source points from baw data
source_points = np.array([baw_face_x, baw_face_y]).T

# Create target grid points from hzg data
target_points = np.array([hzg_node_x, hzg_node_y]).T

# Interpolate baw data onto hzg grid
interpolated_spm = griddata(source_points, baw_spm, target_points, method='linear')
# Note: The 'method' argument can be 'linear', 'nearest', or 'cubic'.


# %%
# set 0- where nan
interpolated_spm[np.isnan(interpolated_spm)] = 0


# %%


# Create source points from baw data
source_points = np.array([baw_face_x, baw_face_y]).T

# Create target grid points from hzg data
target_points = np.array([hzg_node_x, hzg_node_y]).T

triangulation = tri.Triangulation(baw_face_x, baw_face_y)

# set masked values to 0
baw_spm_unmasked = np.ma.filled(baw_spm, 0)
interpolator = tri.LinearTriInterpolator(triangulation, baw_spm_unmasked)

# plot the tripcolor with baw spm
# %matplotlib widget
# fig, ax = plt.subplots()
# tcf = ax.tripcolor(triangulation, baw_spm_unmasked, shading='gouraud', cmap='viridis')
# ax.set_aspect('equal')

# interpolate baw spm to hzg grid
hzg_spm = interpolator(hzg_node_x, hzg_node_y)






# %% [markdown]
# ## Creating a copy of the HZG dataset with the interpolated SPM values

# %%
# for safety make folder and all files in folder read only for all
import os, stat
path_to_dir = '/scratch/local1/hzg'
for root, dirs, files in os.walk(path_to_dir):
    for momo in dirs:
        os.chmod(os.path.join(root, momo), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    for momo in files:
        os.chmod(os.path.join(root, momo), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


# %%
# copy files from hzg to hzg2 and 

import shutil

path_to_dir = '/scratch/local1/hzg'
path_to_dir2 = '/scratch/local1/hzg2'


for root, dirs, files in os.walk(path_to_dir):
    n_files = len(files)
    for momo in dirs:
        os.makedirs(os.path.join(path_to_dir2, momo))
    for ii,momo in enumerate(files):
        print(f'( {ii+1} / {n_files} )')
        shutil.copy(os.path.join(root, momo), path_to_dir2)



# %%
# make files in hzg2 writeable

path_to_dir = '/scratch/local1/hzg2'
for root, dirs, files in os.walk(path_to_dir):
    for momo in dirs:
        os.chmod(os.path.join(root, momo), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    for momo in files:
        os.chmod(os.path.join(root, momo), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)

# %%
import netCDF4 as nc

# Define your turbidity data - for this example, I'm just creating dummy data.
# In reality, you should replace this with your actual turbidity data.
turbidity_data = interpolated_spm

path_to_dir2 = '/scratch/local1/hzg2'

for root, dirs, files in os.walk(path_to_dir2):
    n_files = len(files)
    for ii,momo in enumerate(files):
        if momo.endswith('.nc'):
            print(momo, f'( {ii+1} / {n_files} )')
            with nc.Dataset(os.path.join(root, momo), 'a') as ds:
                # Ensure the dimension exists
                if 'nSCHISM_hgrid_node' not in ds.dimensions:
                    raise ValueError("The dimension 'nSCHISM_hgrid_node' does not exist in the file.")

                # Check if turbidity variable already exists, if not, create it
                if 'turbidity' not in ds.variables:
                    turbidity_var = ds.createVariable('turbidity', 'f4', ('nSCHISM_hgrid_node',))

                # Assign the turbidity data to the variable
                turbidity_var[:] = turbidity_data

                # If needed, set attributes for the turbidity variable
                turbidity_var.units = 'kg/l'  # Adjust this as necessary
                turbidity_var.long_name = 'Water Turbidity'
                # add description
                turbidity_var.description = 'Turbidity based on BAW data from Arne Hammrich. Average of 3D data over time and depth'





