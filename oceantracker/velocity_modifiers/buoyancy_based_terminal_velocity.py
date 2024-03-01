from oceantracker.velocity_modifiers._base_velocity_modifer import VelocityModiferBase
from oceantracker.velocity_modifiers.terminal_velocity import TerminalVelocity
from oceantracker.util.parameter_checking import ParamValueChecker as PVC
from numba import njit
import numpy as np


class BuoyancyBasedTerminalVelocity(TerminalVelocity):
    
    def __init__(self):
        super().__init__()
        # self.add_default_params({})
        

    def check_requirements(self):
        super().check_requirements()
        self.check_class_required_fields_prop_etc(requires3D=True, required_props_list=['buoyancy'])


    def update(self, time_sec, active):
        # modify vertical velocity, if backwards, make negative
        si = self.shared_info
        part_prop = si.classes['particle_properties']
        velocity_modifier = part_prop['velocity_modifier']
        
        self._add_individual_vertical_vel(velocity_modifier.data, part_prop['buoyancy'].data, si.model_direction, active)
        
