import numpy as np
import matplotlib.pyplot as plt
import oceantracker.post_processing.plotting.plot_utilities as plot_utilities

# time series of near bottom or near surface
def plot_relative_height(tracks_data,  particleID =0, ax = plt.gca(), title='', bottom=True,ncase= 0,plot_file_name=None,credit=None):

    if bottom:
        zb=tracks_data['z'] +tracks_data['water_depth']
        status = tracks_data['status'][:, particleID] -10
        statusLab='Status-10'
        ylab= 'Above bottom, m'
    else:
        zb =  tracks_data['z']-tracks_data['tide']
        ylab = 'Below free surface, m'
        status = tracks_data['status'][:, particleID]
        statusLab = 'Status'

    t = tracks_data['time'] / 24. / 3600.
    t = t - t[0]

    ax.plot(t, status, label=statusLab)
    ax.plot(t, zb, linewidth=0.5,alpha =0.5)
    ax.plot(t, zb[:,particleID],label= 'particle', color = 'g')
    ax.set( xlabel='Time, days', ylabel=ylab, title= title)
    ax.legend()

    plot_utilities.show_output(plot_file_name=plot_file_name)

def plot_path_in_vertical_section(tracks_data,  particleID =0,title='', ncase= 0, plot_file_name=None,credit=None):

    threshold = 1

    t = tracks_data['time'] / 24. / 3600.
    t = t - t[0]
    ax = plt.gca()
    ax.plot(t, tracks_data['status'][:, particleID], label='Status')
    ax.plot(t, tracks_data['tide'][:, particleID], label='Tide, m', color='k', linewidth=.5)
    ax.plot(t, tracks_data['tide'][:, particleID]-threshold, label='Tide, m', color='r', linewidth=.5)
    ax.plot(t, -tracks_data['water_depth'][:, particleID], label='Water depth, m', color='k')

    z = tracks_data['z'][:,particleID]
    tide = tracks_data['tide'][:,particleID]
    above = z[z >= tide-threshold]
    below = z[z < tide-threshold]
        
    ax.plot(t[z >= tide-threshold], above,label='Particle z, m', color = 'g')
    ax.plot(t[z < tide-threshold], below,label='Particle z, m', color = 'r')
    # ax.plot(t, tracks_data['z'][:,particleID],label='Particle z, m', color = 'g')

    ax.set(xlabel='Time, days',ylabel='z, m', title=title)
    ax.legend()
    plot_utilities.show_output(plot_file_name=plot_file_name)