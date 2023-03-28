from oceantracker.particle_concentrations._base_user_triangle_properties import _BaseTriangleProperties
from oceantracker.util.parameter_checking import ParamDictValueChecker as PVC
from numba import njit
import numpy as np

class  ParticleConcentrations2D(_BaseTriangleProperties):
    def __init__(self):
        super().__init__()
        # set up info/attributes


    def check_requirements(self):
        self.check_class_required_fields_prop_etc(required_grid_var_list=['triangle_area', 'x'],
                                                            required_props_list=['total_water_depth'])

    def set_up_data_buffers(self):
        si = self.shared_info
        grid = si.classes['reader'].grid
        # set up data buffer
        s = (grid['triangles'].shape[0],)

        self.particle_count = np.full(s, 0, dtype=np.int32)
        self.particle_concentration = np.full(s, 0.)
        self.data_buffers={} # for other particle prop to get concentrations of

    def set_up_output_file(self):
        super().set_up_output_file()
        # add 2D variables
        nc= self.nc
        nc.create_a_variable('particle_count', ['time_dim','triangle_dim'], dtype=self.particle_count.dtype)
        nc.create_a_variable('particle_concentration', ['time_dim', 'triangle_dim'], dtype=self.particle_concentration.dtype)

    def update(self,n_buffer, time):
        si=self.shared_info
        grid = si.classes['reader'].grid

        sel = self.select_particles_to_count()
        self.calcuate_concentration2D(si.classes['particle_properties']['n_cell'].data,
                                      si.classes['particle_properties']['total_water_depth'].data,
                                      grid['triangle_area'],
                                       self.particle_count,
                                        self.particle_concentration,  sel)
        self.write(n_buffer, time)
        self.record_time_stats_last_recorded(time)

    @staticmethod
    @njit()
    def calcuate_concentration2D(n_cell, total_water_depth, triangle_area,particle_count, particle_concentration, sel_to_count):
        particle_count[:] = 0
        particle_concentration[:] = 0.
        for n in sel_to_count:
            c = n_cell[n]
            twd= total_water_depth[n]
            vol = triangle_area[c] * twd
            if vol > .1:  # only calculate  if at least .1 cubic meter
                particle_count[c] += 1
                particle_concentration[c] += 1.0 / vol


class ParticleConcentrationsDepthLayer(ParticleConcentrations2D):

    # def __init__(self):
    #     super().__init__()
    #     self.add_default_params({'layer_depth' : PVC(1, float, min=0)})

    # def check_requirements(self):
    #     super().__init__()
    #     self.check_class_required_fields_prop_etc(required_props_list=['tide'])

    def select_particles_to_count(self):
        si= self.shared_info
        part_prop =  si.classes['particle_properties']
        status_sel = part_prop['status'].compare_all_to_a_value('eq', si.particle_status_flags['moving'], out=self.get_particle_index_buffer())
        distance_from_surface = part_prop['tide'].data - part_prop['x'].data[:,2]
        layer_sel = np.where(distance_from_surface < 1)[0]
        return np.intersect1d(status_sel, layer_sel)