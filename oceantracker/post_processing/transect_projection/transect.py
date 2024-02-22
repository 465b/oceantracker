import numpy as np
from oceantracker.post_processing.read_output_files.load_output_files import load_track_data
from oceantracker.post_processing.read_output_files.load_output_files import get_case_info_files_from_dir
from oceantracker.util.polygon_util import InsidePolygon

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

        self.vertices = vertices


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
        self.fix_transect_line()
        self.length_of_section()
        self.distance_downstream()
        # create the rotation matrix to transform the coordinate space into
        # x1 pointing downstream and x2 pointing to the right hand side of it
        self.build_basis_transform_matrix()


    def project_track_data(self,case_info_dir,n_case,
                           var_list=[],nt_range=None,min_status=0):

        case = get_case_info_files_from_dir(case_info_dir)[n_case]

        default_vars = ['x', 'time','status', 'IDrelease_group', 'IDpulse', 'x0','x_last_good']
        var_list = default_vars + var_list
        
        track_data = load_track_data(case,var_list)

        # This is supposed to take the track_data['x'] and transform 

        #   track_data['x'][time_step,particle_index,(lon,lat,elevation)]
        x = track_data['x'][:, :, :]
        
        # np.nan all the dead/deselected particles
        sel = track_data['status'][:, :] < min_status
        x[sel] = np.nan

        transect_track_data = np.zeros_like(x)*np.nan

        for ii,section in enumerate(self.transect):

            for nt in range(len(x)):
            #   each point in time
                particles_in_poly = section['polygon'].is_inside(x[nt])
            
                # slice track_data down to those inside poly and save their indicies
                sub_track_data = x[nt,particles_in_poly]

                local_tracks = self.transform_basis(
                    sub_track_data,section['rotation_matrix'],section['plane'].vertices[0])

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

        track_data['x'] = transect_track_data

        self.track_data = track_data


    def distance_downstream(self):
        for nn,section in enumerate(self.transect):
            if nn == 0:
                section['distance_downstream'] = 0
            else:
                distance_downstream = 0
                for ii in range(nn):
                    distance_downstream += self.transect[ii]['length']
                section['distance_downstream'] = distance_downstream


    def length_of_section(self):

        for section in self.transect:
            section['length'] = np.linalg.norm(
                section['plane'].vertices[1]-section['plane'].vertices[0]
                )


    def build_basis_transform_matrix(self):
       
       # we use a transformation of the form
       #    x'_i =  A_ij (X_j + B_i)
        for section in self.transect:
            
            # fetch 2d plane vertices
            A = section['plane'].vertices[0]
            B = section['plane'].vertices[1]

            AB = B-A


            standard_base = [[1,0,0],[0,1,0],[0,0,1]]

            # here we start to assume that we only rotate along the verticle axis
            # to generalize angle_of_vectors needs to be adjusted and the transform
            #  needs to be a sum of rotations instead of only a single one
            theta = -self.angle_of_vectors(standard_base[0],AB)
            # print(f'theta {theta}')
            

            rotation_matrix = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta),  0],
                [0            , 0,              1]
            ])
            
            downstream_base_rotated = np.dot(rotation_matrix,AB/np.linalg.norm(AB))

            section['rotation_matrix'] = rotation_matrix

    @staticmethod
    def transform_basis(x,rotation_matrix,base_origin):

        if x.shape[1] == 3:
            x_prime = np.zeros((x.shape[0],3))
            for ii in np.arange(x.shape[0]):
                x_prime[ii] = np.dot(rotation_matrix,x[ii] - base_origin)
                # x_prime[ii] = np.dot(rotation_matrix,x[ii]) + A
            
            return x_prime

        else:
            raise ValueError


    def build_polygons(self):
        
        for section in self.transect:
            section['polygon'] = InsidePolygon(section['polygon_vertices'])


    def build_planes(self):
        for section in self.transect:
            section['plane'] = Plane(section['plane_vertices'])


    def fix_transect_line(self):

        # plane_poly_intersection
        # find intersection of the projection plane with the selection polygon
        # to use intersection as new start&stop of projection plane 

        for section in self.transect:
            intersections = []        
            # for plane (there is only one)
            plane_line_segment = Line(section['plane'].vertices[0,:2], 
                                      section['plane'].vertices[1,:2])

            for n in range(len(section['polygon'].points)-1):
                points = section['polygon'].points[
                    [n,  (n+1) % len(section['polygon'].points) ],:]
                polygon_line_segment = Line(points[0],points[1])

                intersection = Line.intersection(plane_line_segment, polygon_line_segment)
                if intersection is not None: intersections.append(intersection)
            
            section['plane'].intersections = intersections

        # oriant_plane_downstream
        # orient the plane intersection such that the first one is the upstream
        # one based on the user input (who is supposed to define upstream first)
        
        for section in self.transect:

            x1 = section['plane'].intersections[0]
            x2 = section['plane'].intersections[1]

            A = section['plane'].vertices[0,:2]
            B = section['plane'].vertices[1,:2]
    

            # alpha as a measure how far "downstream" from A
            alpha_1 = ((x1 - A)/(B-A))[0]
            alpha_2 = ((x2 - A)/(B-A))[0]

            if alpha_1 < alpha_2:
                section['plane'].vertices[0,:2] = x1
                section['plane'].vertices[1,:2] = x2
            else:
                section['plane'].vertices[0,:2] = x2
                section['plane'].vertices[1,:2] = x1


    
    @staticmethod
    def two_to_three_dim(x):
        if x.shape == (2,):
            return np.append(x,0)
        elif x.shape[1] == 2:
            padding = np.zeros((x.shape[0],1))
            return np.concatenate((x,padding),axis=1)


    @staticmethod
    def angle_of_vectors(a,b):

        a = a[0:2]
        b = b[0:2]

        inner = np.inner(a, b)
        norms = np.linalg.norm(a) * np.linalg.norm(b)

        cos = inner / norms
        rad = np.arccos(np.clip(cos, -1.0, 1.0))

        sign = np.sign(np.linalg.det(np.swapaxes(np.stack((a,b)),0,1)))

        return sign*rad
