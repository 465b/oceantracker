#%%
from oceantracker.post_processing.transect_projection.transect import Transect
from oceantracker.post_processing.plotting import plot_transects

import numpy as np

#%%
path_to_dir = '/scratch/local1/output/22_08_19_testing_stranding_no_resus_v01'
#%%

test_polygone = [{'points': np.array([[561516,5933791],
                    [566084,5933766],
                    [565935,5932215],
                    [561500,5932557]])}]

test_plane = [{'points': np.array([
                                [570000,5933500],
                                [561000,5933000]
                             ])}]

polygone = [
    {'points':
        np.array([
            [569979,5930284],
            [570239,5930483],
            [569147,5931534],
            [568979,5931251]
        ])
    },
    {'points':
        np.array([
            [568988,5931254],
            [569147,5931537],
            [567622,5932312],
            [567502,5932043]
        ])
    },
    {'points':
        np.array([
            [567584,5931996],
            [567704,5932265],
            [565924,5932656],
            [565900,5932289]
        ])
    },
    {'points':
        np.array([
            [565900,5932282],
            [565943,5932663],
            [564086,5933401],
            [563971,5933043]
        ])
    },
    {'points':
        np.array([
            [564000,5933064],
            [564096,5933417],
            [562586,5933178],
            [562643,5932757]
        ])
    },
    {'points':
        np.array([
            [562638,5932757],
            [562614,5933141],
            [561359,5933128],
            [561359,5932690]
        ])
    }
]

planes = [
    {'points':
        np.array([
            [570181,5930237],
            [568926,5931456]
        ])
    },
    {'points':
        np.array([
            [569123,5931325],
            [567444,5932218]
        ])
    },
    {'points':
        np.array([
            [567771,5932080],
            [565804,5932477]
        ])
    },
    {'points':
        np.array([
            [565977,5932423],
            [563961,5933259]
        ])
    },
    {'points':
        np.array([
            [564207,5933266],
            [562446,5932925]
        ])
    },
    {'points':
        np.array([
            [562704,5932919],
            [561174,5932947]
        ])
    }
]

transect = [{'polygon_vertices': poly['points'],'plane_vertices': plane['points']} for poly,plane in zip(polygone,planes)]

#%%
# 

#%%
norderelbe_transect = Transect(transect)
norderelbe_transect.project_track_data(path_to_dir,var_list=['water_depth','age'],n_case=0)

#%%

nt = 30
axis_lims = [560500,572500,5.923e6,5.935e6]
plot_transects.plot_transect_map(norderelbe_transect,nt=nt,axis_lims=axis_lims,plot_file_name='transect_map.png')
plot_transects.plot_projected_horizontal_tracks(norderelbe_transect,nt=nt)
plot_transects.plot_projected_verticle_tracks(norderelbe_transect,nt=nt)

# %%
