import numpy as np
from oceantracker.user_trajectory_modifiers._trajectory_modifers_base import TrajectoryModifiersBase
from oceantracker.util.polygonUtil import  InsidePolygon
from oceantracker.util.parameterObject import ParamDictValueChecker as PVC

class cull_in_polygon(TrajectoryModifiersBase):
    # fallows particles to freeze if inside a polygon
    def __init__(self):
        # set up info/attributes
        super().__init__()  # required in children to get parent defaults
        self.add_default_params({'name': PVC('cull_in_polygon', str),
                                 'polygon': {'points': PVC(None,'vector', is_required=True)},
                                 'fraction_to_cull': PVC(0.,float)
                                 })
        self.polygons = []

    def initialize(self, **kwargs):

        super().initialize()

        # set up polygons to test if particles inside
        if 'points' not in self.params['polygon']:
            self.exit('Polygon settlement, each polygon must be a dictionary with at least a "points"  key a list of corrdinates')

        a = np.asarray(self.params['polygon']['points'])
        if a.shape[1] != 2:
            self.throw_error('parameter polygon_ points must be a list of lists of cordinate pairs, eg [[p1_x1,p1_y1],[p1_x2,p1_y2],..]]')

        self.polygon = InsidePolygon(verticies=a)# do set up to speed inside using pre-calculation

        # add particle prop to track which are inside polygon, which will be automatically written to output
        particle= self.shared_info.pointers['particleGroupManager']
        particle.create_particle_property(dict(name='is_culled_in_polygon', type='manual_update', dtype=np.int8))

        # ad a parameter to record when last released
        particle.create_particle_property(dict(name='time_of_release', type='manual_update', initial_value=0.))

    def can_be_added(self):
        # required method to decide if the operation can be add to simulation
        return True

    # all particles checked to see if they need status changing
    def update(self, buffer_index, time, active):
        si = self.shared_info
        p = si.particles

        # find those inside and freeze, only if they haven't been recently frozen
        those_inside = self.polygon.inside_indices(p.prop_dict['x'].dataInBufferPtr(), out= self.get_particle_index_buffer(), active=active)
        if those_inside.shape[0] > 0:
            last_released = p.prop_dict['time_of_release'].get_values(those_inside)

            # selecting particles to cull based on culling ratio
            particles_to_cull = np.random.choice(those_inside, int(those_inside.shape[0] * self.params['fraction_to_cull']), replace=False)

            # record those frozen
            p.prop_dict['status'].set_values(p.status_flags['dead'], particles_to_cull)
            p.prop_dict['is_culled_in_polygon'].set_values(1, particles_to_cull)
            p.prop_dict['time_of_release'].set_values(0., particles_to_cull)