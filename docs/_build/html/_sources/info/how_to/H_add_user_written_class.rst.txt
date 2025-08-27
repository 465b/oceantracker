Add your own class
==================

[This note-book is in oceantracker/tutorials_how_to/]

First create your own class create as a python file in the same diretory
as your run script. Your class must in import the base case for the role
it performs or one of its children in that role.

.. container::

.. code:: ipython3

    # minimal_example.py using class helper method
    #------------------------------------------------
    from oceantracker.main import OceanTracker
    
    # make instance of oceantracker to use to set parameters using code, then run
    ot = OceanTracker()
    
    # ot.settings method use to set basic settings
    ot.settings(output_file_base='add_user_written_class', # name used as base for output files
                root_output_dir='output',             #  output is put in dir   'root_output_dir'\\'output_file_base'
                time_step= 120. #  2 min time step as seconds
                )
    # ot.set_class, sets parameters for a named class
    ot.add_class('reader',input_dir= './demo_hindcast/schsim3D',  # folder to search for hindcast files, sub-dirs will, by default, also be searched
                          file_mask=  'demo_hindcast_schisim3D*.nc')  # hindcast file mask
    # add  release locations from two points,
    #               (ie locations where particles are released at the same times and locations)
    # note : can add multiple release groups
    ot.add_class('release_groups', 
                        points= [[1595000, 5482600],        #[x,y] pairs of release locations
                                 [1599000, 5486200]],      # must be an N by 2 or 3 or list, convertible to a numpy array
                        release_interval= 3600,           # seconds between releasing particles
                        pulse_size= 10,                   # number of particles released each release_interval
                )
    
    # add user written particle property in same dir as notebook
    
    ot.add_class('particle_properties',name='on_bottom_time',
                                            class_name='my_part_prop.TimeAtStatus')
    
    # run oceantracker
    case_info_file_name = ot.run()
    
    # output now in folder output/minimal_example
    # case_info_file_name the name a json file with useful info for post processing, eg output file names
    print(case_info_file_name)
    


.. parsed-literal::

    helper: ----------------------------------------------------------------------
    helper: Starting OceanTracker helper class,  version 0.50.03.0100-2025-08-13 
    helper:      Python version: 3.10.18 | packaged by conda-forge | (main, Jun  4 2025, 14:45:41) [GCC 13.3.0]
    helper: ----------------------------------------------------------------------
    helper: OceanTracker version 0.50.03.0100-2025-08-13  starting setup helper "main.py":
    helper: Output is in dir "/mnt/c/Users/Laurin.Steidle/OneDrive - Cawthron/Documents/projects/oceantracker_dev/oceantracker/tutorials_how_to/output/minimal_example"
    helper:     hint: see for copies of screen output and user supplied parameters, plus all other output
    helper:     >>> Note: to help with debugging, parameters as given by user  are in "minimal_example_raw_user_params.json"
    helper: >>> Warning: Numba has already been imported, some numba options may not be used (ignore SVML warning)
    helper:     hint: Ensure any code using Numba is imported after Oceantracker is run, eg Oceantrackers "load_output_files.py" and "read_ncdf_output_files.py"
    helper: ----------------------------------------------------------------------
    helper: Numba setup: applied settings, max threads = 8, physical cores = 8
    helper:     hint:  cache code = False, fastmath= False
    helper: ----------------------------------------------------------------------
    helper:       - Built OceanTracker package tree,	  0.368 sec
    helper:       - Built OceanTracker sort name map,	  0.000 sec
    helper:   - Done package set up to setup ClassImporter,	  0.368 sec
    setup: ----------------------------------------------------------------------
    setup:  OceanTracker version 0.50.03.0100-2025-08-13 
    setup:     Starting user param. runner: "minimal_example" at  2025-08-27T10:02:25.798314
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
    setup:       - built domain and island outlines,	  0.455 sec
    setup:       - calculated triangle areas,	  0.000 sec
    setup:       - Finished grid setup
    setup:       - built barycentric-transform matrix,	  0.000 sec
    setup:   - Loading reader fields ['water_depth', 'water_velocity', 'tide']
    setup:   - Finished field group manager and readers setup,	  1.007 sec
    setup: ----------------------------------------------------------------------
    setup:   - Added 1 release group(s) and found run start and end times,	  0.006 sec
    setup:   - Done initial setup of all classes,	  0.014 sec
    setup: ----------------------------------------------------------------------
    setup:   - Starting "minimal_example",  duration: 0 days 23 hrs 0 min 0 sec
    setup:       From 2017-01-01T00:30:00 to  2017-01-01T23:30:00
    setup:       Time step 120.0 sec
    setup:         using: A_Z_profile = False bottom_stress = False
    setup: ----------------------------------------------------------------------
    setup:   -  Reading 24 time steps,  for hindcast time steps 00:23 into ring buffer offsets 000:023 ,  for run "minimal_example"
    setup:       -  read  24 time steps in  0.2 sec, from ./demo_hindcast/schsim3D 
    setup:   - Starting time stepping: 2017-01-01T00:30:00 to 2017-01-01T23:30:00
    setup:     duration  0 days 23 hrs 0 min 0 sec, time step=  0 days 0 hrs 2 min 0 sec 
    S:   - Opened tracks output and done written first time step in: "minimal_example_tracks_compact_000.nc",	  0.093 sec
    S: 0000: 00%:H0000b00-01 Day +00 00:00 2017-01-01 00:30:00: Rel:20   : Active:20     Move:20     Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 4%  step time = 686.9 ms
    S: 0001: 00%:H0000b00-01 Day +00 00:02 2017-01-01 00:32:00: Rel:20   : Active:20     Move:20     Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 4%  step time = 2484.7 ms
    S: 0031: 04%:H0001b01-02 Day +00 01:02 2017-01-01 01:32:00: Rel:40   : Active:40     Move:40     Bottom:0     Strand:0      Dead:0     Out:   0 Buffer: 8%  step time = 12.7 ms
    S: 0061: 09%:H0002b02-03 Day +00 02:02 2017-01-01 02:32:00: Rel:60   : Active:60     Move:60     Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:12%  step time =  2.3 ms
    S: 0091: 13%:H0003b03-04 Day +00 03:02 2017-01-01 03:32:00: Rel:70   : Active:70     Move:60     Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:14%  step time =  3.2 ms
    S: 0121: 18%:H0004b04-05 Day +00 04:02 2017-01-01 04:32:00: Rel:80   : Active:80     Move:70     Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:16%  step time =  2.4 ms
    S: 0151: 22%:H0005b05-06 Day +00 05:02 2017-01-01 05:32:00: Rel:90   : Active:90     Move:80     Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:18%  step time =  4.9 ms
    S: 0181: 26%:H0006b06-07 Day +00 06:02 2017-01-01 06:32:00: Rel:100  : Active:100    Move:90     Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:20%  step time =  3.4 ms
    S: 0211: 31%:H0007b07-08 Day +00 07:02 2017-01-01 07:32:00: Rel:110  : Active:110    Move:100    Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:22%  step time =  3.4 ms
    S: 0241: 35%:H0008b08-09 Day +00 08:02 2017-01-01 08:32:00: Rel:120  : Active:120    Move:110    Bottom:0     Strand:10     Dead:0     Out:   0 Buffer:25%  step time =  4.8 ms
    S: 0271: 39%:H0009b09-10 Day +00 09:02 2017-01-01 09:32:00: Rel:140  : Active:140    Move:140    Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:29%  step time =  1.8 ms
    S: 0301: 44%:H0010b10-11 Day +00 10:02 2017-01-01 10:32:00: Rel:160  : Active:160    Move:160    Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:33%  step time =  2.5 ms
    S: 0331: 48%:H0011b11-12 Day +00 11:02 2017-01-01 11:32:00: Rel:180  : Active:180    Move:180    Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:37%  step time =  4.3 ms
    S: 0361: 52%:H0012b12-13 Day +00 12:02 2017-01-01 12:32:00: Rel:200  : Active:200    Move:200    Bottom:0     Strand:0      Dead:0     Out:   0 Buffer:41%  step time =  2.4 ms
    S: 0391: 57%:H0013b13-14 Day +00 13:02 2017-01-01 13:32:00: Rel:220  : Active:220    Move:206    Bottom:0     Strand:14     Dead:0     Out:   0 Buffer:45%  step time =  3.2 ms
    S: 0421: 61%:H0014b14-15 Day +00 14:02 2017-01-01 14:32:00: Rel:240  : Active:240    Move:219    Bottom:0     Strand:21     Dead:0     Out:   0 Buffer:50%  step time =  3.6 ms
    S: 0451: 65%:H0015b15-16 Day +00 15:02 2017-01-01 15:32:00: Rel:250  : Active:250    Move:192    Bottom:0     Strand:58     Dead:0     Out:   0 Buffer:52%  step time =  1.9 ms
    S: 0481: 70%:H0016b16-17 Day +00 16:02 2017-01-01 16:32:00: Rel:260  : Active:260    Move:199    Bottom:0     Strand:61     Dead:0     Out:   0 Buffer:54%  step time =  2.5 ms
    S: 0511: 74%:H0017b17-18 Day +00 17:02 2017-01-01 17:32:00: Rel:270  : Active:270    Move:201    Bottom:0     Strand:69     Dead:0     Out:   0 Buffer:56%  step time =  4.1 ms
    S: 0541: 78%:H0018b18-19 Day +00 18:02 2017-01-01 18:32:00: Rel:280  : Active:280    Move:211    Bottom:0     Strand:69     Dead:0     Out:   0 Buffer:58%  step time =  2.6 ms
    S: 0571: 83%:H0019b19-20 Day +00 19:02 2017-01-01 19:32:00: Rel:290  : Active:290    Move:223    Bottom:0     Strand:67     Dead:0     Out:   0 Buffer:60%  step time =  1.5 ms
    S: 0601: 87%:H0020b20-21 Day +00 20:02 2017-01-01 20:32:00: Rel:300  : Active:300    Move:236    Bottom:0     Strand:64     Dead:0     Out:   0 Buffer:62%  step time =  2.6 ms
    S: 0631: 91%:H0021b21-22 Day +00 21:02 2017-01-01 21:32:00: Rel:320  : Active:320    Move:299    Bottom:0     Strand:21     Dead:0     Out:   0 Buffer:66%  step time =  4.0 ms
    S: 0661: 96%:H0022b22-23 Day +00 22:02 2017-01-01 22:32:00: Rel:340  : Active:340    Move:326    Bottom:0     Strand:14     Dead:0     Out:   0 Buffer:70%  step time =  5.8 ms
    S: --- Closing all classes ----------------------------------------------
    S: Converting compact track files to rectangular format (to disable set reader param convert=False)
    S:     - Finished "minimal_example_tracks_rectangular_000.nc",	  0.449 sec
    S:     - Conversion complete,	  0.449 sec
    end: ----------------------------------------------------------------------
    end: Finished "minimal_example",  started: 2025-08-27 10:02:25.368104, ended: 2025-08-27 10:02:35.613188
    end:       Computational time =0:00:10.245094
    end:     Timings: total =  10.2 sec
    end:         Setup                        3.88 s	 37.8%
    end:         Reading hindcast             1.40 s	 13.7%
    end:         Initial cell guess           1.13 s	 11.0%
    end:         RK integration               2.05 s	 20.1%
    end:         Find horizontal cell         2.19 s	 21.4%
    end:         Find vertical cell           0.88 s	  8.6%
    end:         Interpolate fields           1.77 s	 17.3%
    end:         Update statistics            0.00 s	  0.0%
    end:         Update custom particle prop. 0.27 s	  2.6%
    end:         resuspension                 0.33 s	  3.2%
    end:         dispersion                   0.41 s	  4.0%
    end:         tracks_writer                0.08 s	  0.7%
    end:     0 errors,   1 warnings,   1 notes
    end: --- Finished: output in "/mnt/c/Users/Laurin.Steidle/OneDrive - Cawthron/Documents/projects/oceantracker_dev/oceantracker/tutorials_how_to/output/minimal_example" 
    /mnt/c/Users/Laurin.Steidle/OneDrive - Cawthron/Documents/projects/oceantracker_dev/oceantracker/tutorials_how_to/output/minimal_example/minimal_example_caseInfo.json
    

Read and plot output
--------------------

A first basic plot of particle tracks

.. code:: ipython3

    from oceantracker.read_output.python import load_output_files
    
    # read particle track data into a dictionary using case_info_file_name
    tracks = load_output_files.load_track_data(case_info_file_name)
    print(tracks.keys()) # show what is in tracks dictionary holds
    
    from oceantracker.plot_output import plot_tracks
    
    ax= [1591000, 1601500, 5479500, 5491000]  # area to plot 
    plot_tracks.plot_tracks(tracks, axis_lims=ax, show_grid=True)


.. parsed-literal::

    Merging rectangular track files
    	 Reading rectangular track file "minimal_example_tracks_rectangular_000.nc"
    dict_keys(['particles_written_per_time_step', 'time_step_range', 'time', 'num_part_released_so_far', 'x', 'x0', 'status', 'status_last_good', 'age', 'ID', 'IDrelease_group', 'user_release_groupID', 'IDpulse', 'hydro_model_gridID', 'time_released', 'on_bottom_time', 'water_depth', 'tide', 'dry_cell_index', 'grid', 'particle_status_flags', 'particle_release_groups', 'axis_lim'])
    


.. image:: H_add_user_written_class_files%5CH_add_user_written_class_3_1.png

