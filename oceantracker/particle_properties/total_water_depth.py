from oceantracker.particle_properties._base_particle_properties import ParticleProperty
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from numba import njit
from oceantracker.util.numba_util import njitOT
import numpy as np

class TotalWaterDepth(ParticleProperty):
    def __init__(self):
        super().__init__()
        self.add_default_params({'time_varying': PVC(True,bool),
                                 'is3D': PVC(False,bool)})

    def check_requirements(self):
        si = self.shared_info
        self.check_class_required_fields_prop_etc(required_props_list=['tide', 'water_depth'])

    def update(self,active):
        si = self.shared_info
        self.get_time_dependent_total_water_depth_from_tide_and_water_depth(
            si.classes['particle_properties']['tide'].data,
            si.classes['particle_properties']['water_depth'].data,
            self.data,
            active)


    @staticmethod
    @njitOT
    def get_time_dependent_total_water_depth_from_tide_and_water_depth(tide, water_depth, total_water_depth, active):
        # get total time dependent water depth as 4D field  from top and bottom cell of LSC grid zlevels
        for n in active:
            total_water_depth[n] = tide[n] + water_depth[n]