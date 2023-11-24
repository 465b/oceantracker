from oceantracker.fields._base_field import CustomFieldBase
from oceantracker.util.parameter_checking import ParamValueChecker as PVC

from numba import njit
import numpy as np

class FrictionVelocity(CustomFieldBase):
    def __init__(self):
        super().__init__()
        self.add_default_params({'is_time_varying': PVC(True,bool),
                                 'num_components': PVC(1, int),
                                 'is3D': PVC(False,bool)})


    def check_requirements(self):
        si = self.shared_info

        self.check_class_required_fields_prop_etc(
            required_fields_list=['water_velocity'],
            requires3D=True)

    def update(self, buffer_index):
        si = self.shared_info
        grid = si.classes['reader'].grid
        fields = si.classes['fields']
        if 'sigma' in grid:
            # sigma model
            self.calc_friction_velocity_from_sigma_levels(buffer_index,
                                                          grid['sigma'],
                                                          fields['total_water_depth'].data,
                                                          fields['water_velocity'].data,
                                                          si.z0, self.data)
        else:
            # native vertical grid
            self.calc_friction_velocity_from_native_zlevels(buffer_index, grid['zlevel'], grid['bottom_cell_index'], si.z0, fields['water_velocity'].data, self.data)
        pass

    @staticmethod
    @njit()
    def calc_friction_velocity_from_sigma_levels(buffer_index, sigma, total_water_depth,
                                                 water_velocity, z0, out):
        # get friction velocity from bottom cell, if velocity is zero at base of bottom cell
        # based on log layer  u= u_* log(z/z0)/kappa
        for nt in buffer_index:
            for n in np.arange(out.shape[1]):  # loop over nodes
                # size of bottom cell from its fraction of the water depth
                dz = (sigma[1] - sigma[0]) * total_water_depth[nt, n, 0, 0]
                speed = np.sqrt(water_velocity[nt, n, 1, 0] ** 2 + water_velocity[nt, n, 1, 1] ** 2)
                out[nt, n, 0, 0] = 0.4 * speed / np.log((dz + z0) / z0)
                # will give np.inf for very thin lower layers, ie small total water depth

    @staticmethod
    @njit()
    def calc_friction_velocity_from_native_zlevels(buffer_index, zlevel, bottom_cell_index, z0, water_velocity, out):
        # get friction velocity from bottom cell, if velocity is zero at base of bottom cell
        # based on log layer  u= u_* log(z/z0)/kappa
        for nt in buffer_index:
            for n in np.arange(zlevel.shape[1]): # loop over nodes
                nz1=bottom_cell_index[n]+1
                dz =  zlevel[nt, n, nz1] - zlevel[nt, n, bottom_cell_index[n]] # size of bottom cell

                if dz >  0.2:
                    speed = np.sqrt(water_velocity[nt, n, nz1, 0]**2 + water_velocity[nt, n, nz1, 1]**2)
                    out[nt, n, 0, 0] = 0.4*speed/np.log((dz+z0)/z0)
                else:
                    out[nt, n, 0, 0] = 0.

