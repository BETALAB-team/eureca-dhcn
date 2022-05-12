__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math
import logging
from eureca_dhcs.exceptions import DuplicateBranch


class Branch:

    """
    Class for the network Branch
    """

    _idx_list = []
    _counter = 0
    default_starting_temperature = 80  # Celsius

    def __init__(
        self,
        idx: str,
        supply_node: str,
        demand_node: str,
        pipe_diameter: float,  # [m]
        pipe_len=None,  # [m]
        roughness=None,  # [-]
        starting_temp=None,  # [°C]
        nodes_objects_dict=None,
    ):
        self._idx = idx
        self._supply_node_idx = supply_node
        self._demand_node_idx = demand_node
        self._pipe_diameter = pipe_diameter

        # Check if pipe len is passed
        if pipe_len == None:
            try:
                supply_node_obj = nodes_objects_dict[supply_node]
                demand_node_obj = nodes_objects_dict[demand_node]
                self._pipe_len = supply_node_obj.distance_from_node(demand_node_obj)
            except TypeError:
                raise NameError(
                    f"Branch {idx}. If pipe lenght not passed, the dictionary of the nodes objects must be passed to calculate the pipe lenght"
                )
        else:
            self._pipe_len = pipe_len
        if roughness == None:
            # Default roughness
            self._roughness = 1.0
        else:
            self._roughness = roughness
        if starting_temp == None:
            # Default starting temperature
            self._previous_temp = Branch.default_starting_temperature
        else:
            self._previous_temp = float(starting_temp)
        # set a unique integer for the incidence matrix
        self._unique_matrix_idx = Branch._counter
        # Other useful properties
        self._perimeter = self._pipe_diameter * math.pi
        self._external_area = self._perimeter * self._pipe_len

    @property
    def _idx(self) -> str:
        return self.__idx

    @_idx.setter
    def _idx(self, value: str):
        try:
            value = str(value)
        except ValueError:
            raise TypeError(f"Branch class, idx must be a str: {value}")
        if value in Branch._idx_list:
            raise DuplicateBranch("Duplicate branch id: {value}")
        self.__idx = value

    @property
    def _pipe_diameter(self) -> float:
        return self.__pipe_diameter

    @_pipe_diameter.setter
    def _pipe_diameter(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe diameter must be a float: {value}"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, pipe diameter very high: {value} [m]")
        self.__pipe_diameter = value

    @property
    def _pipe_diameter(self) -> float:
        return self.__pipe_diameter

    @_pipe_diameter.setter
    def _pipe_diameter(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe diameter must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, pipe diameter negative pipe diameter: {value} [m]"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, pipe diameter very high: {value} [m]")
        self.__pipe_diameter = value

    @property
    def _pipe_len(self) -> float:
        return self.__pipe_len

    @_pipe_len.setter
    def _pipe_len(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, pipe lenght must be a float: {value}")
        if value < 0.0:
            raise ValueError(f"Branch {self._idx}, lenght negative: {value} [m]")
        self.__pipe_len = value

    @property
    def _roughness(self) -> float:
        return self.__roughness

    @_roughness.setter
    def _roughness(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe roughness must be a float: {value}"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, roughness very high: {value} [-]")
        self.__roughness = value

    @property
    def _previous_temp(self) -> float:
        return self.__previous_temp

    @_previous_temp.setter
    def _previous_temp(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, starting temperature must be a float: {value}"
            )
        if value > 200.0:
            logging.warning(
                f"Branch {self._idx}, starting temperature very high: {value} [°C]"
            )
        self.__previous_temp = value

    # # First try flow rate
    # # 0.1 m/s flow rate
    # self.first_try_flow_rate = 0.1 * 1000 * self.pipe_diameter**2 / 4
    def get_hydraulic_resistance(self, flow_rate=None):
        # TODO example to try
        # example just to try
        resistance = 1
        return resistance

    def get_thermal_conductance(self):
        # TODO example to try
        U = 2  # [W/m2K]
        return self.external_area * U  # [W/K]

    # def set_timestep_flow_rate(self, flow):
    #     self.timestep_flow_rate = flow

    # def get_thermal_capacity(self):
    #     V = math.pi * self.pipe_len * self.pipe_diameter ** 2 / 4  # [m3]
    #     return V * rho * cp  # [kJ/K]

    # def get_flow_rate_cp(self):
    #     return self.timestep_flow_rate * cp
