from re import A
from matplotlib.ft2font import HORIZONTAL
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

cdf = Dataset('/scratch/local1/hzg2/schout_1.nc')
diff = cdf.variables['diffusivity']


for ii in range(diff.shape[2]):
    print(ii)

    layer = diff[:,:,ii]
    layer = layer.flatten()

    plt.figure()
    plt.hist(layer,bins=100)
    plt.title(f'diffusivity at vlayer {ii}')
    plt.savefig('/home/zmaw/u301513/Documents/scr/phd/bicest/oceantracker/experiments/auxiliaries/diffusivity_hist_'+str(ii)+'.png')
    plt.xlabel('diffusivity bin')

for ii in np.linspace(0,diff.shape[1],20,dtype=int):
    print(ii)
    plt.figure()
    for jj in range(diff.shape[0]):
        print(jj)
        x = diff[jj,ii,]
        y = np.arange(len(x))
        plt.plot(x,y)
    plt.xlabel('diffusivity')
    plt.ylabel('verticle layer')
    plt.title(f'diffusivity at node {ii} every hour for 24h')
    plt.savefig('/home/zmaw/u301513/Documents/scr/phd/bicest/oceantracker/experiments/auxiliaries/diffusivity_vert_'+str(ii)+'.png')

# diff = diff[diff != 1e-6]