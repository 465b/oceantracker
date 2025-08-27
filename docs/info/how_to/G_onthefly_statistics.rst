On-the-fly statistics
=====================

[This note-book is in oceantracker/tutorials_how_to/]

Scaling up particle numbers to millions will create large volumes of
particle track data. Storing and analyzing these tracks is slow and
rapidly becomes overwhelming. For example, building a heat map from a
terabyte of particle tracks after a run has completed. Ocean tracker can
build some particle statistics on the fly, without recording any
particle tracks. This results in more manageable data volumes and
analysis.

On-the-fly statistics record particle counts separately for each release
group. It is also possible to subset the counts, ie only count particles
which are stranded by the tide by designating a range of particle status
values to count. Or, only count particles in a given vertical “z” range.
Users can add multiple statistics, all calculated in from the same
particles during the run. Eg. could add a particle statistic for each
status type, for different depth ranges.

Statistics can be read, plotted or animated with OceanTrackers
post-processing code, see below

The available “particle_statistics” classes with their individual
settings are at …. add link

Currently there are two main classes of 2D particle statistics “gridded”
which counts particles inside cells of a regular grid, and “polygon”
which counts particles in a given list of polygons.

The user can add many particle statistics classes, all based on the same
particles. For both types it is possible to only count a subset of these
particles, by setting a min. and/or max status to count, or setting a
min. and/or max. “z”, the vertical location. So could add several
statistics classes, each counting particles in different layers, or
classes to separately count those moving and those on the bottom hoping
to be re-suspended.

Gridded statistics
------------------

These are heat maps of counts binned into cells of a regular grid. Along
with heat maps of particle counts, users can optionally build a heat
maps of named particle properties, eg. the value decaying particle
property. To ensure the heat map grids are not too large or too coarse,
by default grids are centred on each release group, thus there are
different grid locations for each release group.

Polygon statistics
------------------

These particle counts can be used to calculate the connectivity between
each release group and a user given list of “statistics” polygons. Also,
used to estimate the influence of each release group on a particle
property with each given statistics polygon. Polygon statistics count
the particles from each point or polygon release within each statistics
polygons. The statistics polygons are are completely independent of the
polygons that might be used in any polygon release (they can be the same
if the user gives both the same point coordinates). A special case of a
polygon statistic, is the “residence_time” class, which can be used to
calculate the fraction of particles from each release group remaining
within each statistics polygon at each ‘update_interval’ as one way to
estimate particle residence time for each release group.

Particle property statistics
----------------------------

Both types of statistics can also record sums of user designated
particle properties within each given grid cell or statistics polygon,
which originate from each release group. These sums enabling mean values
of designated particle properties within each grid cell or polygon to be
calculated. They can also be used to estimate the relative influence of
each release group on the value of a particle property within each given
grid cell or polygon.

A future version with allow estimating the variance of the designated
property values and particle counts in each grid cell or polygon, for
each release group.

Gridded/Heat map example
------------------------

The below uses the helper class method to extends the minimal_example to
add

-  Decaying particle property, eg. breakdown of a pollutant
-  Gridded time series of particle statistics as heat maps, which also
   builds a heat map of the pollutant
-  Plot the particle counts and pollutant as animated heatmap.

.. code:: ipython3

    # Gridded Statistics example.py using class helper method
    #------------------------------------------------
    from oceantracker.main import OceanTracker
    
    # make instance of oceantracker to use to set parameters using code, then run
    ot = OceanTracker()
    
    # ot.settings method use to set basic settings
    ot.settings(output_file_base='heat_map_example', # name used as base for output files
                root_output_dir='output',             #  output is put in dir   'root_output_dir'\\'output_file_base'
                time_step= 600., #  10 min time step as seconds
                write_tracks = False # particle tracks not needed for on fly 
                )
    # ot.set_class, sets parameters for a named class
    ot.add_class('reader',input_dir=  './demo_hindcast/schsim3D',  # folder to search for hindcast files, sub-dirs will, by default, also be searched
                          file_mask=  'demo_hindcast_schisim3D*.nc')  # hindcast file mask
    
    # add one release locations 
    ot.add_class('release_groups', 
                   name='my_release_point', # optional name used to refere to group in plotting
                    points= [ [1599000, 5486200]],       # ust be 1 by N list pairs of release locations
                    release_interval= 900,           # seconds between releasing particles
                    pulse_size= 1000,                   # number of particles released each release_interval
                )
    # add a decaying particle property
    # add and Age decay particle property, with exponential decay based on age, with time scale 1 hour                             
    ot.add_class('particle_properties', # add a new property to particle_properties role
                name ='a_pollutant', # must have a user given name to 
                class_name='oceantracker.particle_properties.age_decay.AgeDecay', #  class_role is resuspension
                initial_value= 1000,
                decay_time_scale = 3600.) # time scale of age decay ie decays initial_value* exp(-age/decay_time_scale)
    
    # add a gridded particle statistic 
    ot.add_class('particle_statistics', 
                    name = 'my_heatmap',
                    class_name= 'oceantracker.particle_statistics.gridded_statistics2D.GriddedStats2D_timeBased',
                    # the below settings are optional
                    update_interval = 900, # time interval in sec, between doing particle statists counts 
                    particle_property_list = ['a_pollutant'], # request a heat map for the decaying part. prop. added above
                    status_list =['moving'], # only count the particles which are moving 
                    z_min =-2.,  # only count particles at locations above z=-2m
                    grid_size= [220, 221],  # number of east and north cells in the heat map
                    grid_span = [20000, 10000],
                    release_group_centered_grids= True
                    )
    
    
    # run oceantracker
    case_info_file_name = ot.run()


.. parsed-literal::

    prelim:     Starting package set up
    helper: ----------------------------------------------------------------------
    helper: Starting OceanTracker helper class,  version 0.50.03.0100-2025-08-13 
    helper:      Python version: 3.11.10 | packaged by conda-forge | (main, Oct 16 2024, 01:17:14) [MSC v.1941 64 bit (AMD64)]
    helper: >>> Warning: Oceantracker is compatible with Python 3.11,  however not all external imported packages have been updated to be compatible with 3.11
    helper:     hint: Down grade to python 3.10 if unexplained issues in external packages
    helper: ----------------------------------------------------------------------
    helper: OceanTracker version 0.50.03.0100-2025-08-13  starting setup helper "main.py":
    helper: Output is in dir "f:\H_Local_drive\ParticleTracking\oceantracker\tutorials_how_to\output\heat_map_example"
    helper:     hint: see for copies of screen output and user supplied parameters, plus all other output
    helper:     >>> Note: to help with debugging, parameters as given by user  are in "heat_map_example_raw_user_params.json"
    helper: ----------------------------------------------------------------------
    helper: Numba setup: applied settings, max threads = 32, physical cores = 32
    helper:     hint:  cache code = False, fastmath= False
    helper: ----------------------------------------------------------------------
    loading oceantracker read files
    helper:       - Built OceanTracker package tree,	  1.234 sec
    helper:       - Built OceanTracker sort name map,	  0.000 sec
    helper:   - Done package set up to setup ClassImporter,	  1.234 sec
    setup: ----------------------------------------------------------------------
    setup:  OceanTracker version 0.50.03.0100-2025-08-13 
    setup:     Starting user param. runner: "heat_map_example" at  2025-08-26T14:10:48.164809
    setup: ----------------------------------------------------------------------
    setup:   - Start  field group manager and readers setup
    setup:   - Found input dir "./demo_hindcast/schsim3D"
    setup:   - Detected reader class_name = "oceantracker.reader.SCHISM_reader.SCHISMreader"
    setup:     Hydro-model is "3D", type "SCHISMreader"
    setup:         hint: Files found in dir and sub-dirs of "./demo_hindcast/schsim3D"
    setup:         Geographic coords = "False" 
    setup:         Hindcast start: 2017-01-01T00:30:00  end:  2017-01-01T23:30:00
    setup:           time step = 0 days 1 hrs 0 min 0 sec, number of time steps= 24 
    setup:           grid bounding box = [1589789.000 5479437.000] to [1603398.000 5501640.000]
    setup:           has:  A_Z profile=True  bottom stress=False
    setup: ----------------------------------------------------------------------
    setup:       - Starting grid setup
    setup:       - built node to triangles map,	  1.067 sec
    setup:       - built triangle adjacency matrix,	  0.194 sec
    setup:       - found boundary triangles,	  0.000 sec
    setup:       - built domain and island outlines,	  1.150 sec
    setup:       - calculated triangle areas,	  0.000 sec
    setup:       - Finished grid setup
    setup:       - built barycentric-transform matrix,	  0.274 sec
    setup:   - Loading reader fields ['water_velocity', 'water_depth', 'tide']
    setup:   - Finished field group manager and readers setup,	  4.949 sec
    setup: ----------------------------------------------------------------------
    setup:   - Added 1 release group(s) and found run start and end times,	  1.887 sec
    setup:   - Done initial setup of all classes,	  0.364 sec
    setup: ----------------------------------------------------------------------
    setup:   - Starting "heat_map_example",  duration: 0 days 23 hrs 0 min 0 sec
    setup:       From 2017-01-01T00:30:00 to  2017-01-01T23:30:00
    setup:       Time step 600.0 sec
    setup:         using: A_Z_profile = False bottom_stress = False
    setup: ----------------------------------------------------------------------
    setup:   -  Reading 24 time steps,  for hindcast time steps 00:23 into ring buffer offsets 000:023 ,  for run "heat_map_example"
    setup:       -  read  24 time steps in  2.9 sec, from ./demo_hindcast/schsim3D 
    setup:   - Starting time stepping: 2017-01-01T00:30:00 to 2017-01-01T23:30:00
    setup:     duration  0 days 23 hrs 0 min 0 sec, time step=  0 days 0 hrs 10 min 0 sec 
    S: 0000: 00%:H0000b00-01 Day +00 00:00 2017-01-01 00:30:00: Rel:1,000: Active:1,000  Move:1,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 1%  step time = 13321.8 ms
    S: 0001: 01%:H0000b00-01 Day +00 00:10 2017-01-01 00:40:00: Rel:1,000: Active:1,000  Move:1,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 1%  step time = 5304.7 ms
    S: 0007: 05%:H0001b01-02 Day +00 01:10 2017-01-01 01:40:00: Rel:4,000: Active:4,000  Move:4,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 5%  step time =  3.7 ms
    S: 0013: 09%:H0002b02-03 Day +00 02:10 2017-01-01 02:40:00: Rel:7,000: Active:7,000  Move:7,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:10%  step time =  4.4 ms
    S: 0019: 14%:H0003b03-04 Day +00 03:10 2017-01-01 03:40:00: Rel:10,000: Active:10,000 Move:10,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:14%  step time =  4.6 ms
    S: 0025: 18%:H0004b04-05 Day +00 04:10 2017-01-01 04:40:00: Rel:13,000: Active:13,000 Move:13,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:18%  step time =  5.0 ms
    S: 0031: 22%:H0005b05-06 Day +00 05:10 2017-01-01 05:40:00: Rel:16,000: Active:16,000 Move:16,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:22%  step time =  6.2 ms
    S: 0037: 27%:H0006b06-07 Day +00 06:10 2017-01-01 06:40:00: Rel:19,000: Active:19,000 Move:19,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:27%  step time =  6.3 ms
    S: 0043: 31%:H0007b07-08 Day +00 07:10 2017-01-01 07:40:00: Rel:22,000: Active:22,000 Move:22,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:31%  step time =  6.4 ms
    S: 0049: 36%:H0008b08-09 Day +00 08:10 2017-01-01 08:40:00: Rel:25,000: Active:25,000 Move:25,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:35%  step time =  6.7 ms
    S: 0055: 40%:H0009b09-10 Day +00 09:10 2017-01-01 09:40:00: Rel:28,000: Active:28,000 Move:28,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:40%  step time =  7.1 ms
    S: 0061: 44%:H0010b10-11 Day +00 10:10 2017-01-01 10:40:00: Rel:31,000: Active:31,000 Move:31,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:44%  step time =  6.8 ms
    S: 0067: 49%:H0011b11-12 Day +00 11:10 2017-01-01 11:40:00: Rel:34,000: Active:34,000 Move:34,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:48%  step time =  7.4 ms
    S: 0073: 53%:H0012b12-13 Day +00 12:10 2017-01-01 12:40:00: Rel:37,000: Active:37,000 Move:37,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:52%  step time =  6.9 ms
    S: 0079: 57%:H0013b13-14 Day +00 13:10 2017-01-01 13:40:00: Rel:40,000: Active:40,000 Move:40,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:57%  step time =  7.6 ms
    S: 0085: 62%:H0014b14-15 Day +00 14:10 2017-01-01 14:40:00: Rel:43,000: Active:43,000 Move:42,972 Bottom:0     Strand:28     Dead:0     Out:   0 Buffer:61%  step time =  7.6 ms
    S: 0091: 66%:H0015b15-16 Day +00 15:10 2017-01-01 15:40:00: Rel:46,000: Active:46,000 Move:45,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:65%  step time =  8.5 ms
    S: 0097: 70%:H0016b16-17 Day +00 16:10 2017-01-01 16:40:00: Rel:49,000: Active:49,000 Move:48,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:70%  step time =  7.6 ms
    S: 0103: 75%:H0017b17-18 Day +00 17:10 2017-01-01 17:40:00: Rel:52,000: Active:52,000 Move:51,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:74%  step time =  8.3 ms
    S: 0109: 79%:H0018b18-19 Day +00 18:10 2017-01-01 18:40:00: Rel:55,000: Active:55,000 Move:54,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:78%  step time =  8.3 ms
    S: 0115: 83%:H0019b19-20 Day +00 19:10 2017-01-01 19:40:00: Rel:58,000: Active:58,000 Move:57,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:82%  step time =  8.4 ms
    S: 0121: 88%:H0020b20-21 Day +00 20:10 2017-01-01 20:40:00: Rel:61,000: Active:61,000 Move:60,827 Bottom:0     Strand:173    Dead:0     Out:   0 Buffer:87%  step time =  9.0 ms
    S: 0127: 92%:H0021b21-22 Day +00 21:10 2017-01-01 21:40:00: Rel:64,000: Active:64,000 Move:63,982 Bottom:0     Strand:18     Dead:0     Out:   0 Buffer:91%  step time =  8.4 ms
    S: 0133: 96%:H0022b22-23 Day +00 22:10 2017-01-01 22:40:00: Rel:67,000: Active:67,000 Move:67,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:95%  step time =  8.6 ms
    S: --- Closing all classes ----------------------------------------------
    end: ----------------------------------------------------------------------
    end: Finished "heat_map_example",  started: 2025-08-26 14:10:46.647960, ended: 2025-08-26 14:11:23.064098
    end:       Computational time =0:00:36.416138
    end:     Timings: total =  36.4 sec
    end:         Setup                        1.51 s	  4.2%
    end:         Reading hindcast             2.93 s	  8.0%
    end:         Initial cell guess           1.54 s	  4.2%
    end:         RK integration               2.37 s	  6.5%
    end:         Find horizontal cell         1.63 s	  4.5%
    end:         Find vertical cell           1.10 s	  3.0%
    end:         Interpolate fields           1.56 s	  4.3%
    end:         Update statistics            0.96 s	  2.6%
    end:         Update custom particle prop. 0.28 s	  0.8%
    end:         resuspension                 0.36 s	  1.0%
    end:         dispersion                   0.48 s	  1.3%
    end:     0 errors,   1 warnings,   1 notes
    end: --- Finished: output in "f:\H_Local_drive\ParticleTracking\oceantracker\tutorials_how_to\output\heat_map_example" 
    


.. image:: G_onthefly_statistics_files%5CG_onthefly_statistics_2_1.png


Read and plot heat maps
~~~~~~~~~~~~~~~~~~~~~~~

The statistics output from the above run is in file
output:raw-latex:`\heat`\_map_example:raw-latex:`\heat`\_map_example_stats_gridded_time_my_heatmap.nc

This netcdf file can be read and organized as a python dictionary by
directly with read_ncdf_output_files.read_stats_file.

To plot use, load_output_files.load_stats_data, which also loads grid
etc for plotting

.. code:: ipython3

    # read stats files
    from oceantracker.read_output.python import read_ncdf_output_files, load_output_files
    from oceantracker.plot_output import plot_statistics
    from IPython.display import HTML
    from os import path
    from oceantracker.util import json_util
    
    # basic read of net cdf, first get file name from case_info.json
    case_info = json_util.read_JSON(case_info_file_name)
    stats_file = path.join(case_info['output_files']['run_output_dir'], case_info['output_files']['particle_statistics']['my_heatmap'])
    raw_stats = read_ncdf_output_files.read_stats_file(stats_file)
    
    print('raw_stats', raw_stats.keys())
    
    # better,  load netcdf plus grid and other data useful in plotting 
    # uses case_info name returned from run above
    stats_data = load_output_files.load_stats_data(case_info_file_name,'my_heatmap')
    print('stats',stats_data.keys())
    
    # use stats_data variable to plot heat map at last time step, by default plots var= "count"
    ax= [1591000, 1601500, 5478500, 5491000] 
    anim= plot_statistics.animate_heat_map(stats_data, release_group_name='my_release_point', axis_lims=ax,
                        heading='Particle count heatmap built on the fly, no tracks recorded', fps=1)
    HTML(anim.to_html5_video())# this is slow to build!
    
    # animate the pollutant
    anim= plot_statistics.animate_heat_map(stats_data, var='a_pollutant',release_group_name= 'my_release_point', axis_lims=ax,
                        heading='Decaying particle property , a_pollutant built on the fly, no tracks recorded', fps=1)
    
    # this line only used in note books, in python scripts use show = True above and set moive file name
    # this is slow to build! 
    HTML(anim.to_html5_video())# this is slow to build!
    
    
    # static heat map
    plot_statistics.plot_heat_map(stats_data, var='a_pollutant',release_group_name= 'my_release_point', axis_lims=ax,  heading='a_pollutant at last time step  depth built on the fly, no tracks recorded')


.. parsed-literal::

    raw_stats dict_keys(['global_attributes', 'dimensions', 'limits', 'variable_attributes', 'grid_spacings', 'y', 'time', 'y_grid', 'num_released', 'cell_area', 'num_released_total', 'x', 'count_all_alive_particles', 'number_released_each_release_group', 'sum_a_pollutant', 'count', 'x_grid', 'time_var', 'date', 'stats_type', 'connectivity_matrix', 'a_pollutant'])
    stats dict_keys(['global_attributes', 'dimensions', 'limits', 'variable_attributes', 'grid_spacings', 'y', 'time', 'y_grid', 'num_released', 'cell_area', 'num_released_total', 'x', 'count_all_alive_particles', 'number_released_each_release_group', 'sum_a_pollutant', 'count', 'x_grid', 'time_var', 'date', 'stats_type', 'connectivity_matrix', 'a_pollutant', 'particle_status_flags', 'particle_release_groups', 'grid'])
    animate_heat_map> colour axis limits [np.int64(0), np.int64(1057)] [np.int64(0), np.int64(1057)]
    

::


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    Cell In[2], line 22
         20 # use stats_data variable to plot heat map at last time step, by default plots var= "count"
         21 ax= [1591000, 1601500, 5478500, 5491000] 
    ---> 22 anim= plot_statistics.animate_heat_map(stats_data, release_group_name='my_release_point', axis_lims=ax,
         23                     heading='Particle count heatmap built on the fly, no tracks recorded', fps=1)
         24 HTML(anim.to_html5_video())# this is slow to build!
         26 # animate the pollutant
    

    File f:\h_local_drive\particletracking\oceantracker\oceantracker\plot_output\plot_statistics.py:49, in animate_heat_map(stats_data, release_group_name, var, axis_lims, credit, interval, heading, vmin, vmax, show_grid, title, logscale, caxis, cmap, movie_file, fps, dpi, back_ground_depth, back_ground_color_map)
         44 plot_utilities.draw_base_map(stats_data['grid'], ax=ax, axis_lims=axis_lims, show_grid=show_grid, title=title, credit=credit,
         45                              back_ground_depth=back_ground_depth, back_ground_color_map=back_ground_color_map)
         47 plot_utilities.plot_release_points_and_polygons(stats_data, ax= ax, release_group_name=release_group_name)
    ---> 49 plot_utilities.show_particleNumbers(stats_data['total_num_particles_released'])
         50 plot_utilities.add_heading(heading)
         52 time_text = plt.text(.05, .05, '', transform=ax.transAxes,c='k', zorder=5)
    

    KeyError: 'total_num_particles_released'



.. image:: G_onthefly_statistics_files%5CG_onthefly_statistics_4_2.png


Polygon example
---------------

add polygon stats example with plotting
=======================================

.. code:: ipython3

    # Polygon Statistics example.py run using dictionary of parameters
    #------------------------------------------------
    # make instance of oceantracker to use to set parameters using code, then run
    from oceantracker.main import OceanTracker
    ot = OceanTracker()
    
    # ot.settings method use to set basic settings
    ot.settings(output_file_base='heat_map_example', # name used as base for output files
                root_output_dir='output',             #  output is put in dir   'root_output_dir'\\'output_file_base'
                time_step= 600., #  10 min time step as seconds
                write_tracks = False # particle tracks not needed for on fly 
                )
    # ot.set_class, sets parameters for a named class
    ot.add_class('reader',input_dir=  './demo_hindcast/schsim3D',  # folder to search for hindcast files, sub-dirs will, by default, also be searched
                          file_mask=  'demo_hindcast_schisim3D*.nc')  # hindcast file mask
    
    # add one release locations 
    ot.add_class('release_groups', 
                   name='my_release_point', # optional name used to refere to group in plotting
                    points= [ [1599000, 5486200]],       # ust be 1 by N list pairs of release locations
                    release_interval= 900,           # seconds between releasing particles
                    pulse_size= 1000,                   # number of particles released each release_interval
                )
    # add a decaying particle property
    # add and Age decay particle property, with exponential decay based on age, with time scale 1 hour                             
    ot.add_class('particle_properties', # add a new property to particle_properties role
                name ='a_pollutant', # must have a user given name to 
                class_name='oceantracker.particle_properties.age_decay.AgeDecay', #  class_role is resuspension
                initial_value= 1000,
                decay_time_scale = 3600.) # time scale of age decay ie decays initial_value* exp(-age/decay_time_scale)
    
    # add a gridded particle statistic 
    ot.add_class('particle_statistics',
                   name='my_polygon',
                   class_name= 'PolygonStats2D_timeBased',
                   polygon_list = [
                            dict(points= [   [1597682.1237, 5489972.7479],# list of one or more polygons
                                     [1598604.1667, 5490275.5488],
                                     [1598886.4247, 5489464.0424],
                                     [1597917.3387, 5489000],
                                     [1597300, 5489000], [1597682.1237, 5489972.7479]]
                            )] 
                   )
    
    # run oceantracker
    poly_case_info_file_name = ot.run()
    


.. parsed-literal::

    helper: ----------------------------------------------------------------------
    helper: Starting OceanTracker helper class,  version 0.50.03.0100-2025-08-13 
    helper:      Python version: 3.11.10 | packaged by conda-forge | (main, Oct 16 2024, 01:17:14) [MSC v.1941 64 bit (AMD64)]
    helper: >>> Warning: Oceantracker is compatible with Python 3.11,  however not all external imported packages have been updated to be compatible with 3.11
    helper:     hint: Down grade to python 3.10 if unexplained issues in external packages
    helper: ----------------------------------------------------------------------
    helper: OceanTracker version 0.50.03.0100-2025-08-13  starting setup helper "main.py":
    helper: Output is in dir "f:\H_Local_drive\ParticleTracking\oceantracker\tutorials_how_to\output\heat_map_example"
    helper:     hint: see for copies of screen output and user supplied parameters, plus all other output
    helper:     >>> Note: to help with debugging, parameters as given by user  are in "heat_map_example_raw_user_params.json"
    helper: >>> Warning: Numba has already been imported, some numba options may not be used (ignore SVML warning)
    helper:     hint: Ensure any code using Numba is imported after Oceantracker is run, eg Oceantrackers "load_output_files.py" and "read_ncdf_output_files.py"
    helper: ----------------------------------------------------------------------
    helper: Numba setup: applied settings, max threads = 32, physical cores = 32
    helper:     hint:  cache code = False, fastmath= False
    helper: ----------------------------------------------------------------------
    helper:       - Built OceanTracker package tree,	  0.019 sec
    helper:       - Built OceanTracker sort name map,	  0.000 sec
    helper:   - Done package set up to setup ClassImporter,	  0.019 sec
    setup: ----------------------------------------------------------------------
    setup:  OceanTracker version 0.50.03.0100-2025-08-13 
    setup:     Starting user param. runner: "heat_map_example" at  2025-08-26T14:12:44.475145
    setup: ----------------------------------------------------------------------
    setup:   - Start  field group manager and readers setup
    setup:   - Found input dir "./demo_hindcast/schsim3D"
    setup:   - Detected reader class_name = "oceantracker.reader.SCHISM_reader.SCHISMreader"
    setup:     Hydro-model is "3D", type "SCHISMreader"
    setup:         hint: Files found in dir and sub-dirs of "./demo_hindcast/schsim3D"
    setup:         Geographic coords = "False" 
    setup:         Hindcast start: 2017-01-01T00:30:00  end:  2017-01-01T23:30:00
    setup:           time step = 0 days 1 hrs 0 min 0 sec, number of time steps= 24 
    setup:           grid bounding box = [1589789.000 5479437.000] to [1603398.000 5501640.000]
    setup:           has:  A_Z profile=True  bottom stress=False
    setup: ----------------------------------------------------------------------
    setup:       - Starting grid setup
    setup:       - built node to triangles map,	  0.000 sec
    setup:       - built triangle adjacency matrix,	  0.000 sec
    setup:       - found boundary triangles,	  0.000 sec
    setup:       - built domain and island outlines,	  0.733 sec
    setup:       - calculated triangle areas,	  0.000 sec
    setup:       - Finished grid setup
    setup:       - built barycentric-transform matrix,	  0.000 sec
    setup:   - Loading reader fields ['water_velocity', 'water_depth', 'tide']
    setup:   - Finished field group manager and readers setup,	  0.814 sec
    setup: ----------------------------------------------------------------------
    setup:   - Added 1 release group(s) and found run start and end times,	  0.001 sec
    setup:   - Done initial setup of all classes,	  1.229 sec
    setup: ----------------------------------------------------------------------
    setup:   - Starting "heat_map_example",  duration: 0 days 23 hrs 0 min 0 sec
    setup:       From 2017-01-01T00:30:00 to  2017-01-01T23:30:00
    setup:       Time step 600.0 sec
    setup:         using: A_Z_profile = False bottom_stress = False
    setup: ----------------------------------------------------------------------
    setup:   -  Reading 24 time steps,  for hindcast time steps 00:23 into ring buffer offsets 000:023 ,  for run "heat_map_example"
    setup:       -  read  24 time steps in  0.0 sec, from ./demo_hindcast/schsim3D 
    setup:   - Starting time stepping: 2017-01-01T00:30:00 to 2017-01-01T23:30:00
    setup:     duration  0 days 23 hrs 0 min 0 sec, time step=  0 days 0 hrs 10 min 0 sec 
    S: 0000: 00%:H0000b00-01 Day +00 00:00 2017-01-01 00:30:00: Rel:1,000: Active:1,000  Move:1,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 1%  step time = 4206.3 ms
    S: 0001: 01%:H0000b00-01 Day +00 00:10 2017-01-01 00:40:00: Rel:1,000: Active:1,000  Move:1,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 1%  step time =  2.6 ms
    S: 0007: 05%:H0001b01-02 Day +00 01:10 2017-01-01 01:40:00: Rel:4,000: Active:4,000  Move:4,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 5%  step time =  2.9 ms
    S: 0013: 09%:H0002b02-03 Day +00 02:10 2017-01-01 02:40:00: Rel:7,000: Active:7,000  Move:7,000  Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:10%  step time =  4.0 ms
    S: 0019: 14%:H0003b03-04 Day +00 03:10 2017-01-01 03:40:00: Rel:10,000: Active:10,000 Move:10,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:14%  step time =  4.8 ms
    S: 0025: 18%:H0004b04-05 Day +00 04:10 2017-01-01 04:40:00: Rel:13,000: Active:13,000 Move:13,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:18%  step time =  5.4 ms
    S: 0031: 22%:H0005b05-06 Day +00 05:10 2017-01-01 05:40:00: Rel:16,000: Active:16,000 Move:16,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:22%  step time =  6.1 ms
    S: 0037: 27%:H0006b06-07 Day +00 06:10 2017-01-01 06:40:00: Rel:19,000: Active:19,000 Move:19,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:27%  step time =  6.0 ms
    S: 0043: 31%:H0007b07-08 Day +00 07:10 2017-01-01 07:40:00: Rel:22,000: Active:22,000 Move:22,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:31%  step time =  6.4 ms
    S: 0049: 36%:H0008b08-09 Day +00 08:10 2017-01-01 08:40:00: Rel:25,000: Active:25,000 Move:25,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:35%  step time =  6.2 ms
    S: 0055: 40%:H0009b09-10 Day +00 09:10 2017-01-01 09:40:00: Rel:28,000: Active:28,000 Move:28,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:40%  step time =  6.4 ms
    S: 0061: 44%:H0010b10-11 Day +00 10:10 2017-01-01 10:40:00: Rel:31,000: Active:31,000 Move:31,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:44%  step time =  6.5 ms
    S: 0067: 49%:H0011b11-12 Day +00 11:10 2017-01-01 11:40:00: Rel:34,000: Active:34,000 Move:34,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:48%  step time =  6.5 ms
    S: 0073: 53%:H0012b12-13 Day +00 12:10 2017-01-01 12:40:00: Rel:37,000: Active:37,000 Move:37,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:52%  step time =  7.0 ms
    S: 0079: 57%:H0013b13-14 Day +00 13:10 2017-01-01 13:40:00: Rel:40,000: Active:40,000 Move:40,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:57%  step time =  6.9 ms
    S: 0085: 62%:H0014b14-15 Day +00 14:10 2017-01-01 14:40:00: Rel:43,000: Active:43,000 Move:42,970 Bottom:0     Strand:30     Dead:0     Out:   0 Buffer:61%  step time =  7.4 ms
    S: 0091: 66%:H0015b15-16 Day +00 15:10 2017-01-01 15:40:00: Rel:46,000: Active:46,000 Move:45,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:65%  step time =  7.5 ms
    S: 0097: 70%:H0016b16-17 Day +00 16:10 2017-01-01 16:40:00: Rel:49,000: Active:49,000 Move:48,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:70%  step time =  7.6 ms
    S: 0103: 75%:H0017b17-18 Day +00 17:10 2017-01-01 17:40:00: Rel:52,000: Active:52,000 Move:51,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:74%  step time =  7.5 ms
    S: 0109: 79%:H0018b18-19 Day +00 18:10 2017-01-01 18:40:00: Rel:55,000: Active:55,000 Move:54,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:78%  step time =  7.6 ms
    S: 0115: 83%:H0019b19-20 Day +00 19:10 2017-01-01 19:40:00: Rel:58,000: Active:58,000 Move:57,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:82%  step time =  8.3 ms
    S: 0121: 88%:H0020b20-21 Day +00 20:10 2017-01-01 20:40:00: Rel:61,000: Active:61,000 Move:60,806 Bottom:0     Strand:194    Dead:0     Out:   0 Buffer:87%  step time =  8.4 ms
    S: 0127: 92%:H0021b21-22 Day +00 21:10 2017-01-01 21:40:00: Rel:64,000: Active:64,000 Move:63,983 Bottom:0     Strand:17     Dead:0     Out:   0 Buffer:91%  step time =  7.9 ms
    S: 0133: 96%:H0022b22-23 Day +00 22:10 2017-01-01 22:40:00: Rel:67,000: Active:67,000 Move:67,000 Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:95%  step time =  8.0 ms
    S: --- Closing all classes ----------------------------------------------
    end: ----------------------------------------------------------------------
    end: Finished "heat_map_example",  started: 2025-08-26 14:12:44.431103, ended: 2025-08-26 14:12:56.533813
    end:       Computational time =0:00:12.102710
    end:     Timings: total =  12.1 sec
    end:         Setup                        1.56 s	 12.9%
    end:         Reading hindcast             2.96 s	 24.5%
    end:         Initial cell guess           1.54 s	 12.7%
    end:         RK integration               2.92 s	 24.1%
    end:         Find horizontal cell         1.78 s	 14.7%
    end:         Find vertical cell           1.22 s	 10.1%
    end:         Interpolate fields           1.63 s	 13.5%
    end:         Update statistics            5.20 s	 43.0%
    end:         Update custom particle prop. 0.29 s	  2.4%
    end:         resuspension                 0.02 s	  0.1%
    end:         dispersion                   0.04 s	  0.3%
    end:     0 errors,   2 warnings,   1 notes
    end: --- Finished: output in "f:\H_Local_drive\ParticleTracking\oceantracker\tutorials_how_to\output\heat_map_example" 
    

Read polygon/connectivity statistics
------------------------------------

.. code:: ipython3

    #Read polygon stats and calculate connectivity matrix 
    from oceantracker.read_output.python import load_output_files
    
    poly_stats_data = load_output_files.load_stats_data(poly_case_info_file_name,'my_polygon')
    print('stats',poly_stats_data.keys())
    
    import matplotlib.pyplot as plt
    plt.plot(poly_stats_data['date'], poly_stats_data['connectivity_matrix'][:,0,0])
    plt.title('Connectivity time series between release point and polygon')
    
    #print(poly_stats_data['date'])


.. parsed-literal::

    stats dict_keys(['global_attributes', 'dimensions', 'limits', 'variable_attributes', 'time', 'num_released', 'count_all_alive_particles', 'num_released_total', 'number_released_each_release_group', 'count', 'time_var', 'date', 'stats_type', 'polygon_list', 'connectivity_matrix', 'particle_status_flags', 'particle_release_groups', 'grid'])
    



.. parsed-literal::

    Text(0.5, 1.0, 'Connectivity time series between release point and polygon')




.. image:: G_onthefly_statistics_files%5CG_onthefly_statistics_8_2.png


Time verses Age statistics
--------------------------

Both gridded and polygon statistics come in two types, “time” and “age”.

-  “time” statistics are time series, or snapshots, of particle numbers
   and particle properties at a time interval given by
   “calculation_interval” parameter. Eg. gridded stats showing how the
   heat map of a source’s plume evolves over time.

-  “age” statistics are particle counts and properties binned by
   particle age. The result are age based histograms of counts or
   particle proprieties. This is useful to give numbers in each age band
   arriving at a given grid cell or polygon, from each release group.
   Eg. counting how many larvae are old enough to settle in a polygon or
   grid cell from each potential source location.
