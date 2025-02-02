from oceantracker.fields._base_field import CustomFieldBase
import numpy as np
from numba import njit
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from oceantracker.util.numba_util import njitOT

class VerticalGradient(CustomFieldBase):

    def __init__(self):
        super().__init__()
        self.add_default_params({'name_of_field': PVC(None, str, is_required=True),
                                 # below are not required as acquired from named field
                                 'time_varying': PVC(True, bool, is_required=False),
                                 'is3D': PVC(True, bool, is_required=False),
                                 })
        self.class_doc(description='Calculated a vertical gradient field with name  "name_of_field" param, as a field named "name_of_field_vertical_grad"')

    def initial_setup(self, grid, fields):
        si = self.shared_info
        # get fields prop from named field
        params= self.params
        field_name = params['name_of_field']
        params['time_varying'] =fields[field_name].is_time_varying() # match base field

        if field_name not in fields:
            si.msg_logger.msg(f'Field vertical gradient >> can not find field {field_name} to setup its vertical gradient class', fatal_error=True, exit_now=True)

        super().initial_setup(grid,fields)  # set up self.data with above params
        pass

    def check_requirements(self):
        self.check_class_required_fields_prop_etc(requires3D=True,)
    def update(self,fields,grid,active):
        si = self.shared_info
        if 'sigma' in grid:
            _calc_field_vert_grad_from_sigma_levels(fields[self.params['name_of_field']].data, grid['sigma'],
                                               fields['tide'].data,fields['water_depth'].data,
                                               grid['bottom_cell_index'], si.z0, fields[self.info['name']].data)
        else:
            # z levels
            _calc_field_vert_grad_from_zlevels(fields[self.params['name_of_field']].data,grid['zlevel'],
                                    grid['bottom_cell_index'], si.z0, fields[self.info['name']].data)

@njitOT
def _calc_field_vert_grad_from_zlevels(field4D,zlevel,bottom_cell_index,z0,gradient_field):

    for nt in range(field4D.shape[0]):
        for node  in  range(field4D.shape[1]):
            for nz in  range(bottom_cell_index[node],field4D.shape[2]-1):
                dz = zlevel[nt,node,nz+1] - zlevel[nt,node,nz]
                if dz > z0:
                    for ncomp in range(field4D.shape[3]):
                        gradient_field[nt, node, nz, ncomp] = (field4D[nt, node, nz+1, ncomp] - field4D[nt, node, nz, ncomp])/dz
                else:
                    gradient_field[nt, node, nz, :] = 0.

                # top cell, assume gradient same as cell below
                gradient_field[nt, node, -1, :] = gradient_field[nt, node, -2, :]

@njitOT
def _calc_field_vert_grad_from_sigma_levels(field4D,sigma, tide, water_depth,bottom_cell_index,z0,gradient_field):

    for nt in range(field4D.shape[0]):
        for node  in  range(field4D.shape[1]):
            twd = abs(tide[nt,node,0,0] +water_depth[0,node,0,0])

            for nz in  range(bottom_cell_index[node],field4D.shape[2]-1):
                dz = (sigma[nz+1] - sigma[nz]) * twd
                if dz < 0.1: dz = 0.1
                if dz > z0:
                    for ncomp in range(field4D.shape[3]):
                        gradient_field[nt, node, nz, ncomp] = (field4D[nt, node, nz+1, ncomp] - field4D[nt, node, nz, ncomp])/dz
                else:
                    gradient_field[nt, node, nz, :] = 0.

                # top cell, assume gradient same as cell below
                gradient_field[nt, node, -1, :] = gradient_field[nt, node, -2, :]
