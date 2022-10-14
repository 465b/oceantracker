#%%
from oceantracker.post_processing.read_output_files.load_output_files import load_particle_track_vars
from oceantracker.post_processing.read_output_files.load_output_files import get_case_info_files_from_dir
from oceantracker.post_processing.plotting import plot_utilities
from oceantracker.util import time_util

from copy import deepcopy

from oceantracker.util.polygon_util import InsidePolygon

import numpy as np
import matplotlib.pyplot as plt




#%%
color_palette={'land': (np.asarray([146, 179, 140])/256).tolist(), 'land_edge': [.5, .5, .5]}
def draw_base_map(grid, ax=plt.gca(), axis_lims=None, back_ground_depth=True,
                  show_grid=False, back_ground_color_map='Blues', title=None, text1=None, credit=None):

    # get grid bounds to fill a recgtangle
    bounds= [np.min(grid['x'][:, 0]), np.max(grid['x'][:, 0]), np.min(grid['x'][:, 1]), np.max(grid['x'][:, 1])]
    dx,dy = bounds[1]- bounds[0], bounds[3]- bounds[2]
    f= 0.05
    bounds =np.asarray([ [bounds[0]-f*dx, bounds[1]+f*dx], [bounds[2]-f*dy,  bounds[3]+f*dy]]) # l
    b = np.asarray([bounds[0,:], [bounds[1,0], bounds[0,1] ], bounds[1, :], [bounds[0,0],bounds[1,1] ], bounds[0,:] ] )

    # fill background land retangle
    ax.fill(b[:,0] , b[:, 1],  facecolor=color_palette['land'],  zorder=0)

    if axis_lims is None: axis_lims= bounds.flatten().tolist()
    ax.set_xlim(axis_lims[:2])
    ax.set_ylim(axis_lims[2:])

    # fill domain as white
    ax.fill(grid['grid_outline']['domain']['points'][:,0], grid['grid_outline']['domain']['points'][:,1],
            edgecolor= None, facecolor=(1., 1., 1.), linewidth=.5, zorder=0)

    # plot islands from outline
    for g in grid['grid_outline']['islands']:
            ax.fill(g['points'][:, 0], g['points'][:, 1], edgecolor=color_palette['land_edge'],
                    facecolor=color_palette['land'], linewidth= 0.5, zorder= 3)
    ax.plot(grid['grid_outline']['domain']['points'][:, 0], grid['grid_outline']['domain']['points'][:, 1], c=color_palette['land_edge'], linewidth=0.5, zorder=3)

    if  back_ground_depth:
        plot_utilities.plot_coloured_depth(grid, ax=ax,
                                           color_map= back_ground_color_map,
                                           zorder=1)

    if show_grid:
        ax.triplot(grid['x'][:, 0], grid['x'][:, 1], grid['triangles'],
                   color=(0.8, 0.8, 0.8), linewidth=.5, zorder=1)
    for o in grid['grid_outline']['open_boundary_nodes']:
        plt.plot(grid['x'][o, 0], grid['x'][o, 1], '--','red')
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])
    # ax.tick_params(axis="both", direction="in", right=True, top=True)

    if title is not None:  ax.set_title(title)
    if text1 is not None:  text_norm(.4, .1, text1, fontsize=8)
    # add_credit(credit)
    plot_utilities.add_map_scale_bar(axis_lims, ax=ax)
    return grid

#%%

def plot_polygon_map(track_data, fraction_to_plot=None, show_grid=False, credit=None,
                     heading =None,title=None, axis_lims=None, back_ground_depth=True,
                     back_ground_color_map= None,plot_file_name=None, 
                     polygon_list_to_plot = None):

    x = track_data['x'][nt, :, :2].copy() # copy so as not to change original data
    sel = track_data['status'][nt, :] < min_status # get rid of dead particles
    x[sel,:] = np.nan
    
    fig = plt.figure()
    ax = plt.gca()

    fig.tight_layout()

    draw_base_map(track_data['grid'],ax=ax, axis_lims=axis_lims, show_grid= show_grid, title=title, credit=credit,
                  back_ground_depth=back_ground_depth,back_ground_color_map= back_ground_color_map)

    ax.scatter(x[:,0],x[:,1],s=10)

    # plot_utilities.plot_release_points_and_polygons(track_data, ax=ax) # these are nominal starts
    plot_utilities.draw_polygon_list(polygon_list_to_plot,ax=ax)

#%%
# PROJECTION

class Line():
    # https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines
    def __init__(self,X1,X2):
        # X1 = (x1,y1)
        # (x1,y1) -> (x2,y2)

        self.X1 = X1
        self.X2 = X2    
        self.A = (X1[1] - X2[1])
        self.B = (X2[0] - X1[0])
        self.C = (X1[0]*X2[1] - X2[0]*X1[1])

    @staticmethod
    def intersection(L1, L2):
        # (L1, L2) = ((X1,Y1))

        #intersection interval
        interval = [max( min(L1.X1[0],L1.X2[0]), min(L2.X1[0],L2.X2[0]) ),
                    min(max(L1.X1[0],L1.X2[0]), max(L2.X1[0],L2.X2[0]) )]

        D  = L1.A * L2.B - L1.B * L2.A
        Dx = L1.C * L2.B - L1.B * L2.C
        Dy = L1.A * L2.C - L1.C * L2.A
        if D != 0:
            x = - Dx / D
            y = - Dy / D

            # fig, ax = plt.subplots(1, 1)
            # plot_utilities.draw_polygon_list(polygone,ax=ax)
            # A = np.stack([L1.X1,L1.X2])
            # plt.plot(A[:,0],A[:,1],linewidth=2.0)
            # B = np.stack([L2.X1,L2.X2])
            # plt.plot(B[:,0],B[:,1],linewidth=3.0)
            # plt.scatter(x=[x],y=[y])
            # plt.savefig('tmp.png')

            if  (x < interval[0]) or (x > interval[1]):
                return
            else:
                return x,y
        else:
            return 


class Plane():
    """ Creates a 2D plane in a 3D space """
    def __init__(self,vertices):
        """
        Takes two points on a plane to calculate an origin and a normal
        """
        if vertices.shape[1] == 2:
            buffer = np.zeros((2,3))
            buffer[:,:2] = vertices
            vertices = buffer


        point_1 = vertices[0,:]
        point_2 = vertices[1,:]

        span_1 = point_1-point_2
        span_1[2] = 0
        span_2 = span_1
        span_2[2] = 1

        self.vertices = vertices
        self.normal = np.cross(span_1,span_2)
        self.origin = point_1


class Transect():

    def __init__(self,transect):
        """
        Input:
        ------
        
        transect: list
            [
                #section 1
                {
                    polygon_vertices: (n,2)-array,
                    plane_vertices: (2,2)-array
                }
                #section 2
                ...
            ]
        """
        self.transect = transect

        self.build_polygons()
        self.build_planes()
        self.fix_planes()
        self.length_of_section()
        self.calc_distance_downstream()
        # create the rotation matrix to transform the coordinate space into
        # x1 pointing downstream and x2 pointing to the right hand side of it
        self.build_transformation()
        print(self)


    def project_track_data(self,track_data,nt_range=None,min_status=0):

        # This is supposed to take the track_data['x'] and transform 

        #   track_data['x'][time_step,particle_index,(lon,lat,elevation)]
        x = track_data['x'][:, :, :]
        transect_track_data = np.zeros_like(x)


        # np.nan all the dead/deselected particles
        sel = track_data['status'][:, :] < min_status
        x[sel] = np.nan
        track_data['x'] = x 

        for ii,section in enumerate(self.transect):

            for nt in range(len(x)):
            #   each point in time
                particles_in_poly = section['polygon'].is_inside(x[nt])
            
                # slice track_data down to those inside poly and save their indicies
                sub_track_data = x[nt,particles_in_poly]
                # particle_index = 'placeholder'

                local_tracks = section['transform'](sub_track_data)

                # at this point the local tracks give the position 
                #   x1 - downstream
                #   x2 - orthogonal to downstream
                #   x3 - elevation
                
                # downstream-position of particle in global transect 
                if ii == 0:
                    global_tracks = local_tracks
                else:
                    global_tracks = local_tracks
                    global_tracks[:,0] += section['distance_downstream']

                transect_track_data[nt,particles_in_poly] = global_tracks

        self.track_data = transect_track_data


    def calc_distance_downstream(self):
        for nn,section in enumerate(self.transect):
            if nn == 0:
                section['distance_downstream'] = 0
            else:
                distance_downstream = 0
                for ii in range(nn-1):
                    distance_downstream += self.transect[ii]['length']
                section['distance_downstream'] = distance_downstream


    def length_of_section(self):

        for section in self.transect:
            section['length'] = np.linalg.norm(
                section['plane_vertices'][1]-section['plane_vertices'][0]
                )


    def build_transformation(self):
       
       # we use a transformation of the form
       #    x'_i =  A_ij (X_j + B_i)
        for section in self.transect:
            
            # fetch 2d plane vertices
            A = section['plane_vertices'][0]
            B = section['plane_vertices'][1]

            # make them 3d
            A = np.append(A,0)
            B = np.append(B,0)


            AB = B-A


            standard_base = [[1,0,0],[0,1,0],[0,0,1]]

            # here we start to assume that we only rotate along the verticle axis
            # to generalize angle_of_vectors needs to be adjusted and the transform
            #  needs to be a sum of rotations instead of only a single one
            theta = self.angle_of_vectors(standard_base[0],AB)
            # print(f'theta {theta}')
            

            rotation_matrix = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta),  0],
                [0            , 0,              1]
            ])
            

            standard_base_rotated = np.dot(rotation_matrix,standard_base[0])
            # print(f'standard base {standard_base[0]}')
            # print(f'downstream base {AB/np.linalg.norm(AB)}')
            # print(f'rotation_matrix {rotation_matrix}')
            # print(f'standard base rotated {standard_base_rotated}')

            def section_transform(x):

                # we assume an input in the shape of (3,n)
                # HERE - this is ugly. reshape.
                if x.shape == (3,):
                    x_prime = np.dot(rotation_matrix,x - A)
                    # x_prime = np.dot(rotation_matrix,x) + A
                    return [x_prime]
                
                elif x.shape == (2,):
                    x = self.two_to_three_dim(x)
                    x_prime = np.dot(rotation_matrix,x - A)
                    # x_prime = np.dot(rotation_matrix,x) + A
                    return [x_prime]

                elif x.shape[1] == 3:
                    x_prime = np.zeros((x.shape[0],3))
                    for ii in np.arange(x.shape[0]):
                        x_prime[ii] = np.dot(rotation_matrix,x[ii] - A)
                        # x_prime[ii] = np.dot(rotation_matrix,x[ii]) + A
                   
                    return x_prime

                elif x.shape[1] == 2:
                    x = self.two_to_three_dim(x)

                    x_prime = np.zeros((x.shape[0],3))
                    for ii in np.arange(x.shape[0]):
                        x_prime[ii] = np.dot(rotation_matrix,x[ii] - A)
                        # x_prime[ii] = np.dot(rotation_matrix,x[ii]) + A
                   
                    return x_prime

                else:
                    raise ValueError


            section['transform'] = section_transform


    def build_polygons(self):
        
        for section in self.transect:
            section['polygon'] = InsidePolygon(section['polygon_vertices'])


    def build_planes(self):
        for section in self.transect:
            section['plane'] = Plane(section['plane_vertices'])


    def fix_planes(self):

        # find intersection of the projection plane with the selection polygon
        # to use intersection as new start&stop of projection plane 
        self.plane_poly_intersection()
        # orient the plane intersection such that the first one is the upstream
        # one based on the user input (who is supposed to define upstream first)
        self.oriente_plane()

    @staticmethod
    def two_to_three_dim(x):
        if x.shape == (2,):
            return np.append(x,0)
        elif x.shape[1] == 2:
            padding = np.zeros((x.shape[0],1))
            return np.concatenate((x,padding),axis=1)


    @staticmethod
    def angle_of_vectors(a,b):
        inner = np.inner(a, b)
        norms = np.linalg.norm(a) * np.linalg.norm(b)

        cos = inner / norms
        rad = np.arccos(np.clip(cos, -1.0, 1.0))
        # deg = np.rad2deg(rad)

        return rad


    #@staticmethod
    def plane_poly_intersection(self):

        for section in self.transect:
            intersections = []        
            # for plane (there is only one)
            plane_line_segment = Line(section['plane'].vertices[0,:2], 
                                      section['plane'].vertices[1,:2])

            for n in range(len(section['polygon'].points)-1):
                points = section['polygon'].points[ [n,  (n+1) % len(section['polygon'].points) ],:]
                polygon_line_segment = Line(points[0],points[1])

                # 
                intersection = Line.intersection(plane_line_segment, polygon_line_segment)
                if intersection is not None: intersections.append(intersection)
            
            section['intersections'] = intersections

    def oriente_plane(self):

        for section in self.transect:

            section['user_plane_vertices'] = deepcopy(section['plane_vertices'])
            x1 = section['intersections'][0]
            x2 = section['intersections'][1]

            A = section['plane_vertices'][0]
            B = section['plane_vertices'][1]

            # alpha as a measure how far "downstream" from A
            alpha_1 = ((x1 - A)/(B-A))[0]
            alpha_2 = ((x2 - A)/(B-A))[0]

            if alpha_1 < alpha_2:
                section['plane_vertices'][0] = x1
                section['plane_vertices'][1] = x2
            else:
                section['plane_vertices'][0] = x2
                section['plane_vertices'][1] = x1
                

    def basis_transform_on_plane(self):
        """ doc """
        for section in self.transect:
            # find right most point of plane-polygon-intersections
            # this point will be the point of origin to measure distance downstream
            origin = 'placeholder'

            x_sub = section['sub_track_data']

            x_transformed = section['plane'].transform_basis(origin,x_sub)
            section['sub_track_data']['x_plane'] = x_transformed


    def transform_on_transect(self,track_data):
        # each section already has x in coordinates relative to the rightmost
        # aka upstream original
        # now we need to "concatinate" the section and sum up the x_transect
        # coordinate to receive a continues transect

        for section in self.transect:
            pass

            section['sub_track_data']['x_transect'] = 'placeholder'

        # put all the individual transects together
        # since they are all disjoint we should be able to simply loop over 
        # them and write them into a new transect_track_data

        transect_track_data = track_data.copy()
        for section in self.transect:
            sel = section['index']
            transect_track_data['x'][:,sel,:] = section['x_transect']


    def transform_plane_basis(self):
        for section in self.transect:
            x = section['sub_track_data']
            pass

#%%
path_to_dir = '/scratch/local1/output/22_08_19_testing_stranding_no_resus_v01'
cases = get_case_info_files_from_dir(path_to_dir)
track_data = load_particle_track_vars(cases[0])
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

#%%
# 
nt = 30
min_status = 0

x = track_data['x'][nt, :, :2].copy() # copy so as not to change original data
sel = track_data['status'][nt, :] < min_status # get rid of dead particles
x[sel,:] = np.nan

track_data['x'][nt, :, :2] = x
plot_polygon_map(track_data,polygon_list_to_plot=polygone+planes)
                             


print('pre projection')
#%%

transect = [{'polygon_vertices': poly['points'],'plane_vertices': plane['points']} for poly,plane in zip(polygone,planes)]
test_transect = Transect(transect)

#%%
test_transect.project_track_data(track_data)


#%%
# poly = InsidePolygon(polygone[0]['points'])

# nt = 30
# min_status = 0

# track_data = load_particle_track_vars(cases[0])

# sel = np.where(~np.isnan(track_data['x'][nt,:,0]) == True)[0]
# sel = np.random.choice(sel)
# sample_point = deepcopy(track_data['x'][nt,sel,:])



# fig, ax = plt.subplots(1, 1)
# plot_utilities.draw_polygon_list(polygone,ax=ax)

# plot_utilities.draw_polygon_list(plane,ax=ax)
# plt.scatter(x=sample_point[0],y=sample_point[1],s=100,color='black')

# A = test_section.transect[0]['plane_vertices'][0]
# B = test_section.transect[0]['plane_vertices'][1]
# AB = B-A

# plt.scatter(A[0],A[1])
# plt.quiver(A[0],A[1],AB[0],AB[1])
# plt.plot([A[0],B[0]],[A[1],B[1]])
# plt.savefig('tmp.png')

#%%

# local_sample_point = test_section.transect[0]['transform'](sample_point)[0]
# local_polygone = [{'points': test_section.transect[0]['transform'](polygone[0]['points']) }]
# local_plane = [{'points': test_section.transect[0]['transform'](plane[0]['points']) }]


# fig, ax = plt.subplots(1, 1)
# plot_utilities.draw_polygon_list(local_polygone,ax=ax)

# plot_utilities.draw_polygon_list(local_plane,ax=ax)
# plt.scatter(x=local_sample_point[0],y=local_sample_point[1],s=100)
# A = local_plane[0]['points'][0]
# B = local_plane[0]['points'][1]
# AB = B-A

# plt.scatter(A[0],A[1])
# plt.quiver(A[0],A[1],AB[0],AB[1])

#%%
# projected_sample_point = deepcopy(local_sample_point)
# projected_sample_point[1] = 0
# print(projected_sample_point)

# fig, ax = plt.subplots(1, 1)
# plot_utilities.draw_polygon_list(local_polygone,ax=ax)
# A = local_plane[0]['points'][0]
# B = local_plane[0]['points'][1]
# AB = B-A
# plt.scatter(A[0],A[1])
# plt.quiver(A[0],A[1],AB[0],AB[1])

# plot_utilities.draw_polygon_list(local_plane,ax=ax)
# plt.scatter(x=projected_sample_point[0],y=projected_sample_point[1],s=100)


#%%

# testing inside_polygon

# nt = 30
# min_status = 0

# x = track_data['x'][nt, :, :2].copy() # copy so as not to change original data
# sel = track_data['status'][nt, :] < min_status # get rid of dead particles
# x[sel,:] = np.nan

# x_inside = x.copy()
# x_inside[poly.is_inside(x)] = x[poly.is_inside(x)]
# x_inside[~poly.is_inside(x)] = np.nan
# x_outside = x.copy()
# x_outside[~poly.is_inside(x)] = x[~poly.is_inside(x)]
# x_outside[poly.is_inside(x)] = np.nan

# track_data['x'][nt, :, :2] = x
# plot_polygon_map(track_data,polygon_list_to_plot=polygone)
# track_data['x'][nt, :, :2] = x_inside
# plot_polygon_map(track_data,polygon_list_to_plot=polygone)
# track_data['x'][nt, :, :2] = x_outside
# plot_polygon_map(track_data,polygon_list_to_plot=polygone)
