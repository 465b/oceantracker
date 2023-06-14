import numpy as np
from numba import njit,prange, types as nbt, typeof, from_dtype
from oceantracker.util.profiling_util import function_profiler
from oceantracker.common_info_default_param_dict_templates import particle_info
# record varaible to hold walk info/couts
# to reduce number of args required in numba functions and be morr readable

# globals
# todo make numpy stucture
status_moving = int(particle_info['status_flags']['moving'])
status_on_bottom = int(particle_info['status_flags']['on_bottom'])
status_stranded_by_tide = int(particle_info['status_flags']['stranded_by_tide'])

status_outside_open_boundary = int(particle_info['status_flags']['outside_open_boundary'])
status_dead = int(particle_info['status_flags']['dead'])
status_bad_cord = int(particle_info['status_flags']['bad_cord'])
status_cell_search_failed = int(particle_info['status_flags']['cell_search_failed'])

#below is called by another numba function which will work out signature on first call
@njit()
def _get_single_BC_cord_numba(x, BCtransform, bc):
    # get BC cord of x for one triangle from DT transform matrix inverse, see scipy.spatial.Delaunay
    # also return index smallest BC for walk and largest
    # returns n_min the index of smallest bc used to choose next triangle
    # bc is (3,) pre-allocated working scale, used to return BC's

    # do (2x2) matrix multiplication of  bc[:2]=BCtransform[:2,:2]*(x-transform[:,2]
    # for i in range(2): bc[i] = 0.
    for i in range(2):
        # for j in range(2):
        #  bc[i] +=  BCtransform[i,j]*(x[j]-BCtransform[2,j])
        # replace loop with faster explicit adds, as no need to zero bc[:] above
        bc[i] = BCtransform[i, 0] * (x[0] - BCtransform[2, 0]) + BCtransform[i, 1] * (x[1] - BCtransform[2, 1])

    bc[2] = 1.0 - bc[0] - bc[1]

    return np.argmin(bc), np.argmax(bc)

# ________ Barycentric triangle walk________
@njit()
def BCwalk_with_move_backs(xq, grid,
                           x_last_good, n_cell,status,bc_cords,
                           active, step_info):
    # Barycentric walk across triangles to find cells

    bc = np.zeros((3,), dtype=np.float64) # working space
    # shortcuts needed to use prange


    # loop over active particles in place
    for nn in prange(active.size):
        n= active[nn]

        if xq[n, 0] == np.nan or xq[n, 1] == np.nan:
            # if any is nan copy all and move on
            _move_back(xq[n,:], x_last_good[n, :])
            step_info['nans_encountered_triangle_walk'] += 1  # count nans
            return

        n_tri = n_cell[n]  # starting triangle
        # do BC walk
        n_steps = 0
        move_back = False

        while n_steps < step_info['max_triangle_walk_steps']:
            # update barcentric cords of xq
            n_min, n_max = _get_single_BC_cord_numba(xq[n, :], grid['bc_transform'][n_tri, :, :], bc)

            if bc[n_min] > -step_info['bc_walk_tol'] and bc[n_max] < 1. + step_info['bc_walk_tol']:
                # are now inside triangle, leave particle status as is
                break  # with current n_tri as found cell

            n_steps += 1
            # move to neighbour triangle at face with smallest bc then test bc cord again
            next_tri = grid['adjacency'][n_tri, n_min]  # n_min is the face num in  tri to move across

            if next_tri < 0:
                # if no new adjacent triangle, then are trying to exit domain at a boundary triangle,
                # keep n_cell, bc  unchanged
                if step_info['open_boundary_type'] > 0 and next_tri == -2:  # outside domain
                    # leave x, bc, cell, location  unchanged as outside
                    status[n] = status_outside_open_boundary
                    break
                else:  # n_tri == -1 outside domain and any future
                    # solid boundary, so just move back
                    move_back = True
                    break

            # check for dry cell
            if step_info['block_dry_cells']:  # is faster split into 2 ifs, not sure why
                if grid['dry_cell_index'][next_tri] > 128:
                    # treats dry cell like a lateral boundary,  move back and keep triangle the same
                    move_back = True
                    break

            n_tri = next_tri

        # not found in given number of search steps
        if n_steps >= step_info['max_triangle_walk_steps']:  # dont update cell
            status[n] = status_cell_search_failed
            # move_back = True# todo shoul it just move back, not retyr?do move back externally

        if move_back:
            # move back dont update
            _move_back(xq[n,:], x_last_good[n, :])
        else:
            # update cell anc BC for new triangle
            n_cell[n] = n_tri
            for i in range(3): bc_cords[n, i] = bc[i]

        step_info['particles_located_by_walking'] += 1  # particles walked
        step_info['number_of_triangles_walked'] += n_steps  # steps taken
        step_info['longest_triangle_walk'] = max(n_steps, step_info['longest_triangle_walk'])  # longest walk

@njit()
def _move_back(x, x_old):
    for i in range(x.shape[0]): x[i] = x_old[i]

@njit
def get_BC_cords_numba(x, n_cells, BCtransform, bc):
    # get BC cords of set of points x inside given cells and return in bc

    for n in range(x.shape[0]):
        _get_single_BC_cord_numba(x[n, :], BCtransform[n_cells[n], :, :], bc[n, :])

@njit
def check_if_point_inside_triangle_connected_to_node(x, node, node_to_tri_map,tri_per_node, BCtransform, bc_walk_tol):
    # get BC cords of set of points x inside given cells and return in bc
    bc = np.zeros((3,), dtype=np.float64)  # working space
    n_cell = np.full((x.shape[0],),-1, np.int32)
    for n in range(x.shape[0]):
        nn = node[n]
        # loop over tri attached to node
        for m in range(tri_per_node[nn]):
            n_tri = node_to_tri_map[nn, m]
            n_min,n_max= _get_single_BC_cord_numba(x[n, :2], BCtransform[n_tri, :, :], bc)
            if    bc[n_min] > -bc_walk_tol and bc[n_max] < 1. + bc_walk_tol:
                # found triangle
                n_cell[n] = n_tri
                continue
    return n_cell

@njit
def get_BC_transform_matrix(points, simplices):
    # pre-build barycectric tranforms for 2D triangles based in scipy spatial qhull as used by scipy.Delauny

    """ from scipy ............
    Compute barycentric affine coordinate transformations for given simplices.
    Returns
    -------
    Tinvs : array, shape (nsimplex, ndim+1, ndim)
        Barycentric transforms for each simplex.
        Tinvs[i,:ndim,:ndim] contains inverse of the matrix ``T``,
        and Tinvs[i,ndim,:] contains the vector ``r_n`` (see below).
    Notes
    -----
    Barycentric transform from ``x`` to ``c`` is defined by::
        T c = x - r_n
    where the ``r_1, ..., r_n`` are the vertices of the simplex.
    The matrix ``T`` is defined by the condition::
        T e_j = r_j - r_n
    where ``e_j`` is the unit axis vector, e.g, ``e_2 = [0,1,0,0,...]``
    This implies that ``T_ij = (r_j - r_n)_i``.
    For the barycentric transforms, we need to compute the inverse
    matrix ``T^-1`` and store the vectors ``r_n`` for each vertex.
    These are stacked into the `Tinvs` returned.
    """

    ndim = 2  # only works on 2D triangles
    nsimplex = simplices.shape[0]

    T = np.empty((ndim, ndim), dtype=np.double)
    Tinvs = np.zeros((nsimplex, ndim + 1, ndim), dtype=np.double)

    for isimplex in range(nsimplex):
        for i in range(ndim):
            Tinvs[isimplex, ndim, i] = points[simplices[isimplex, ndim], i]  # puts cords of last point as extra column, ie r_n vector
            for j in range(ndim):
                T[i, j] = (points[simplices[isimplex, j], i] - Tinvs[isimplex, ndim, i])
            Tinvs[isimplex, i, i] = np.nan

        # form inverse of 2 by 2, https://mathworld.wolfram.com/MatrixInverse.html
        # compute matrix determinate of 2 by 2
        det = T[0, 0] * T[1, 1] - T[0, 1] * T[1, 0]

        # inverse  matrix term by term
        Tinvs[isimplex, 0, 0] = T[1, 1] / det
        Tinvs[isimplex, 1, 1] = T[0, 0] / det
        Tinvs[isimplex, 0, 1] = -T[0, 1] / det
        Tinvs[isimplex, 1, 0] = -T[1, 0] / det

    return Tinvs

#no signature needed as called by a numba function
@njit()
def _eval_z_at_nz_cell( tf,nb, nz_cell, z_level_at_nodes,  nz_bottom_nodes, BCcord, nodes):
    # eval zlevel at particle location and depth cell, return z and nodes required for evaluation
    z = 0.
    for m in range(3):
        nz = max(nz_cell, nz_bottom_nodes[m]) # move up to bottom, so not out of range
        z += z_level_at_nodes[nb[0], nodes[m], nz] * BCcord[m] * tf[1] \
             + z_level_at_nodes[nb[1], nodes[m], nz] * BCcord[m] * tf[0]
    return z

@njit()
def set_hindcast_buffer_steps(time_sec, step_info):
    # get next two buffer time steps around the given time in reader ring buffer
    # plus global time step locations and time ftactions od timre step
    # put results in interpolators step info numpy structure
    hindcast_fraction = (time_sec - step_info['hindcast_first_time']) / (step_info['hindcast_last_time'] - step_info['hindcast_first_time'])
    step_info['current_hydro_model_step'] = (step_info['n_time_steps_in_hindcast'] - 1) * hindcast_fraction # global hindcast time step

    # ring buffer locations of surounding steps
    step_info['nb'][0] = step_info['current_hydro_model_step'] % step_info['time_buffer_size']
    step_info['nb'][1] = (step_info['current_hydro_model_step']  + step_info['model_direction']) % step_info['time_buffer_size']

@njit()
def set_time_fractions(time_sec, time_hindcast, step_info):
    # sets the fraction of time step that current time is between
    # surrounding hindcast time steps
    # abs makes it work when backtracking
    s  = abs(time_sec - time_hindcast) / step_info['hydro_model_time_step']
    step_info['fractional_time_steps'][0] = s
    step_info['fractional_time_steps'][1] = 1.0 - s

@njit()
def get_depth_cell_time_varying_Slayer_or_LSCgrid(xq, grid,
                                                  n_cell, status, bc_cords,nz_cell,z_fraction,z_fraction_bottom_layer,
                                                  active, step_info):
    # find the zlayer for each node of cell containing each particle and at two time slices of hindcast  between nz_bottom and number of z levels
    # LSC grid means must track vertical nodes for each particle
    # nz_with_bottom is lowest cell in grid, is 0 for slayer vertical grids, but may be > 0 for LSC grids
    # nz_with_bottom must be time independent

    for nn in prange(active.size): # loop over active particles
        n = active[nn]

        # vertical walk to search for a particle's layer in the grid, nz_cell
        nb = step_info['nb']
        top_nz_cell = grid['zlevel'].shape[2] - 2
        nodes = grid['triangles'][n_cell[n], :]  # nodes for the particle's cell
        bottom_nz_nodes = grid['bottom_cell_index'][nodes]
        bottom_nz_cell = np.min(bottom_nz_nodes)  # cell at bottom is smallest of those in triangle

        # preserve status if stranded by tide
        if status[n] == status_stranded_by_tide:
            nz_cell[n] = bottom_nz_cell
            # update nodes above and below
            z_below = _eval_z_at_nz_cell(step_info['fractional_time_steps'], step_info['nb'], bottom_nz_cell, grid['zlevel'], grid['bottom_cell_index'], bc_cords[n, :], nodes)
            xq[n, 2] = z_below
            z_fraction[n] = 0.0
            return

        n_vertical_steps = 0

        zq = xq[n, 2]

        # make any already on bottom active, may be flagged on bottom if found on bottom, below
        if status[n] == status_on_bottom: status[n] = status_moving

        # find zlevel above and below  current vertical cell
        nz = nz_cell[n]
        z_below = _eval_z_at_nz_cell(step_info['fractional_time_steps'], nb, nz, grid['zlevel'],grid['bottom_cell_index'],bc_cords[n, :], nodes)
        #print('zz1 ',n_cell[n],status[n],    nz, z_below, zq, bottom_nz_cell, BCcord[n,:])

        if zq >= z_below:
            # search upwards, do nothing if z_above > zq[n] > z_below, ie current nodes are correct
            z_above = _eval_z_at_nz_cell(step_info['fractional_time_steps'], nb, nz + 1,  grid['zlevel'], grid['bottom_cell_index'], bc_cords[n, :], nodes)
            while zq > z_above:
                #print('up ',  nz, z_below, zq, z_above)
                if nz >= top_nz_cell:
                    if zq > z_above:
                        zq = z_above   # clip to free surface height
                    break  # stop if in top cell
                nz += 1
                z_below = z_above  # retain for dz calc
                z_above = _eval_z_at_nz_cell(step_info['fractional_time_steps'], nb, nz, grid['zlevel'], grid['bottom_cell_index'], bc_cords[n, :], nodes)
                n_vertical_steps += 1
        else:
            # search downwards
            z_above  = z_below
            z_below = _eval_z_at_nz_cell(step_info['fractional_time_steps'], nb, nz - 1, grid['zlevel'], bottom_nz_nodes, bc_cords[n, :], nodes)
            while zq < z_below:
                #print('down ', nz, z_below, zq, z_above)
                if nz <= bottom_nz_cell:
                    if zq < z_below:
                        zq = z_below  # clip to bottom depth
                    break  # found cell
                nz -= 1
                z_above = z_below  # retain for dz calc.
                z_below = _eval_z_at_nz_cell(step_info['fractional_time_steps'], nb, nz, grid['zlevel'], grid['bottom_cell_index'], bc_cords[n, :], nodes)
                n_vertical_steps += 1

        # nz now holds required cell
        dz = z_above - z_below
        # get z linear z_fraction
        if dz < grid['z0']:
            z_fraction[n] = 0.0
        else:
            z_fraction[n] = (zq - z_below) / dz

        # extra work if in bottom cell
        z_fraction_bottom_layer[n] = -999.  # flag as not in bottom layer, will become >= 0 if in layer

        if nz == bottom_nz_cell:
            z_bot = z_below
            # set status if on the bottom set status
            if zq < z_bot + grid['z0']:
                status[n] = status_on_bottom
                zq = z_bot

            # get z_fraction for log layer
            if dz < grid['z0']:
                z_fraction_bottom_layer[n] = 0.0
            else:
                # adjust z fraction so that linear interp acts like log layer
                z0p = grid['z0'] / dz
                z_fraction_bottom_layer[n] = (np.log(z_fraction[n] + z0p) - np.log(z0p)) / (np.log(1. + z0p) - np.log(z0p))

        # record new depth cell
        nz_cell[n] = nz
        xq[n, 2] = zq
        #print('zz2 ',status[n],  nz,z_below,zq , z_above, z_fraction[n], z_fraction_bottom_layer[n] ,dz)
        # record number of vertical search steps made for this particle
        # step count stats, tidal stranded particles are not counted
        step_info['total_vertical_steps_walked'] += n_vertical_steps
        step_info['longest_vertical_walk'] = max(step_info['longest_vertical_walk'], n_vertical_steps) # record max number of steps

def _dev_combined_cell_and_depth_walk(xq, grid, part_prop, active, step_info):
    # faster as done as two kernals
    bc = np.zeros((3,), dtype=np.float64) # working space for BC cords

    for nn in prange(active.size):  # loop over active particles
        n = active[nn]
        _kernal_BCwalk_with_move_backs(xq[n, :], grid, part_prop, n, step_info, bc)
        _kernal_get_depth_cell_time_varying_Slayer(xq[n, :], grid, part_prop, step_info, n)


# Below is numpy version of numba BC cord code, now only used as check
#________________________________________________
    def get_cell_cords_check(self,x,n_cell):
        # barycentric cords, only for use with non-improved scipy and KDtree for al time steps
        # numba code does this faster
        si = self.shared_info
        grid = si.classes['reader'].grid

        TT = np.take(grid['bc_transform'], n_cell, axis=0,)
        b = np.full((x.shape[0],3), np.nan, order='C')
        b[:,:2] = np.einsum('ijk,ik->ij', TT[:, :2], x[:, :2] - TT[:, 2], order='C')  # Einstein summation
        b[:,2] = 1. - b[:,:2].sum(axis=1)
        return b

#________ old versions

