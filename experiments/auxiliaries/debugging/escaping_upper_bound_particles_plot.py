
import os
import numpy as np

import netCDF4 as cdf
import json

import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib

from matplotlib.animation import FuncAnimation
from IPython.display import HTML

import oceantracker.post_processing.particlePlot as otPlot
import oceantracker.util.rUtil  as rUtil
import oceantracker.post_processing.loadOutputFiles as loadOutputFiles

ncase = 10

runCaseInfo = loadOutputFiles.load_runCaseInfo('/scratch/local1/output/21_10_09_depth_accuracy_test_v01/21_10_09_depth_accuracy_test_v01_runInfo.json',ncase=ncase)
tracks = loadOutputFiles.load_particle_track_vars(runCaseInfo, ['x', 'status', 'time','water_depth'])

ii = 1
fig,ax0 = plt.subplots()
ax0.plot(tracks['time'].astype('datetime64[s]'),tracks['status'][:,ii],c='green',label='particle status',zorder=0)
ax0.plot(tracks['time'].astype('datetime64[s]'),tracks['x'][:,ii,2],label='vertical position',zorder=1)
ax0.plot(tracks['time'].astype('datetime64[s]'),-tracks['water_depth'][:,ii],c='red',label='water depth')


for label in ax0.xaxis.get_ticklabels():
        label.set_rotation(45)

plt.ylim(-20,15)
plt.legend()
plt.tight_layout()
plt.show()