__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math
import logging

import numpy as np

from eureca_dhcs.exceptions import DuplicateNode, WrongNodeType, WrongTemperatureMode


class Node:
    """
    Class for the network Node
    """

    _idx_list = []
    _counter = 0
    _counter_undefined_flow_rate = 0
    # Starting pressure for hydraulic balance.
    # This is used just for the first timestep,
    # then pressure of the previous timestep is used
    _starting_pressure = 10.0  # [Pa]
    _cooling_starting_temperature = 15  # [°C]
    _heating_starting_temperature = 55  # [°C]

    def __init__(
        self,
        idx: str,
        node_type: str,
        supply_nodes: list,
        demand_nodes: list,
        x: float,
        y: float,
        starting_temperature=None,  # [°C]
        temperature_mode="Heating",
    ):
        """


        Parameters
        ----------
        idx : str
            id of the node. Must be unique.
        node_type : str
            [supply, disp, demand].
        supply_nodes : list
            list of the id of the connected supply nodes.
        demand_nodes : list
            list of the id of the connected demand nodes.
        x : float
            x coordinate [m]
        y : float
            y coordinate [m]
        starting_temperature : float
            default: None
            possible starting temperature [°C]
        temperature_mode : float
            default: Heating
            temperature mode: choose between Heating or Cooling


        Returns
        -------
        None.

        """
        self._idx = idx
        self._node_type = node_type
        self._x = x
        self._y = y
        self._supply_nodes_list = supply_nodes
        self._demand_nodes_list = demand_nodes
        # The following dict will contain objects of type Branch and Node
        self._supply_branches_objects = {}
        self._demand_branches_objects = {}
        self._supply_nodes_objects = {}
        self._demand_nodes_objects = {}

        Node._idx_list.append(self._idx)
        self._unique_matrix_idx = Node._counter
        Node._counter += 1

        self._node_pressure_array = np.array([])
        self._node_pressure = self._starting_pressure  # [Pa]
        self._node_temperature_array = np.array([])
        if starting_temperature != None:
            # Default starting temperature
            self._node_temperature = starting_temperature
        else:
            if temperature_mode == "Heating":
                self._node_temperature = self._heating_starting_temperature
            elif temperature_mode == "Cooling":
                self._node_temperature = self._cooling_starting_temperature
            else:
                raise WrongTemperatureMode(
                    f"Node {self._idx}: temperature mode must be or Heating or Cooling. Temperature mode: {temperature_mode}"
                )

    @property
    def _idx(self) -> str:
        return self.__idx

    @_idx.setter
    def _idx(self, value: str):
        try:
            value = str(value)
        except ValueError:
            raise TypeError(f"Node class, idx must be a str: {value}")
        if value in Node._idx_list:
            print(Node._idx_list)
            raise DuplicateNode(f"Duplicate node id: {value}")
        self.__idx = value

    @property
    def _node_type(self) -> str:
        return self.__node_type

    @_node_type.setter
    def _node_type(self, value: str):
        try:
            value = str(value)
        except ValueError:
            raise TypeError(f"Node {self._idx}, node type must be a str: {value}")
        if value not in ["supply", "disp", "demand"]:
            raise WrongNodeType(f"Node {self._idx}, wrong node type: {value}")
        self.__node_type = value

    @property
    def _node_pressure(self) -> float:
        return self.__node_pressure

    @_node_pressure.setter
    def _node_pressure(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Node {self._idx}, node pressure must be a float: {value} [Pa]"
            )
        if value > 1e10:
            logging.warning(
                f"Node {self._idx} pressure over check boundaries: {value} [Pa]"
            )
        self.__node_pressure = value
        self._node_pressure_array = np.append(self._node_pressure_array, value)

    @property
    def _node_temperature(self) -> float:
        return self.__node_temperature

    @_node_temperature.setter
    def _node_temperature(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Node {self._idx}, node temperature must be a float: {value} [°C]"
            )
        if value > 95:
            logging.warning(f"Node {self._idx} temperature over 95 °C: {value} [°C]")
        if value < 5:
            logging.warning(
                f"Node {self._idx} temperature lower than 5 °C: {value} [°C]"
            )
        self.__node_temperature = value
        self._node_temperature_array = np.append(self._node_temperature_array, value)

    @property
    def _boundary_mass_flow_rate(self) -> np.array:
        return self.__boundary_mass_flow_rate

    @_boundary_mass_flow_rate.setter
    def _boundary_mass_flow_rate(self, value: np.array):
        try:
            value = np.array(value)
        except ValueError:
            raise TypeError(
                f"Node {self._idx}, _boundary_mass_flow_rate must be a np.array: {value}"
            )
        if self._node_type == "demand" and np.any(value < 0):
            logging.warning(
                f"Node {self._idx}: demand node with negative mass flow rate"
            )
        if self._node_type == "supply" and np.any(value > 0):
            logging.warning(
                f"Node {self._idx}: supply node with positive mass flow rate"
            )
        self.__boundary_mass_flow_rate = value

    @property
    def _boundary_temperature(self) -> np.array:
        return self.__boundary_temperature

    @_boundary_temperature.setter
    def _boundary_temperature(self, value: np.array):
        try:
            value = np.array(value)
        except ValueError:
            raise TypeError(
                f"Node {self._idx}, _boundary_temperature must be a np.array: {value}"
            )
        if np.any(value < 5):
            logging.warning(f"Node {self._idx}: supply node temperature less than 5 °C")
        if np.any(value > 95):
            logging.warning(
                f"Node {self._idx}: supply node temperature higher than 95 °C"
            )
        if self._node_type != "supply":
            logging.warning(f"Node {self._idx}: temperature set for a non supply node")
        self.__boundary_temperature = value
        
    @property
    def _boundary_pressure(self) -> np.array:
        return self.__boundary_pressure

    @_boundary_pressure.setter
    def _boundary_pressure(self, value: np.array):
        try:
            value = np.array(value)
        except ValueError:
            raise TypeError(
                f"Node {self._idx}, _boundary_pressure must be a np.array: {value}"
            )
        self.__boundary_pressure = value
        
    @property
    def _boundary_mass_flow_rate_undefined(self) -> bool:
        return self.__boundary_mass_flow_rate_undefined

    @_boundary_mass_flow_rate_undefined.setter
    def _boundary_mass_flow_rate_undefined(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError(
                f"Node {self._idx}, _boundary_mass_flow_rate_undefined must be a bool: {value}"
            )
        self.__boundary_mass_flow_rate_undefined = value
        if value:
            self._unique_matrix_idx_undefined_flow_rate = Node._counter_undefined_flow_rate
            Node._counter_undefined_flow_rate += 1


    def get_supply_branches_unique_idx(self):
        return [
            branch._unique_matrix_idx
            for k, branch in self._supply_branches_objects.items()
        ]

    def get_demand_branches_unique_idx(self):
        return [
            branch._unique_matrix_idx
            for k, branch in self._demand_branches_objects.items()
        ]

    def get_supply_branches_mass_flow_rates(self):

        array = np.array(
            [
                [int(branch._unique_matrix_idx), branch._mass_flow_rate]
                for k, branch in self._supply_branches_objects.items()
            ]
        )
        return array[:, 0].astype(int), array[:, 1]

    def distance_from_node(self, node_2):
        """
        Takes a second Node object to calculate the distance [m]

        Parameters
        ----------
        node_2 : Node
            Second node.

        Returns
        -------
        float
            distance [m].

        """
        return math.sqrt((self._x - node_2._x) ** 2 + (self._y - node_2._y) ** 2)

    # def get_supply_temperature(self):
    #     return random.random() * 10 + 70

    # def add_supply_branch(self, branch):
    #     self.supply_branches[branch.idx] = branch

    # def get_entering_flow_rate(self):
    #     total_enetering_flow_rate = 0.0
    #     for supply_branch in self.supply_branches.values:
    #         total_enetering_flow_rate += supply_branch.timestep_flow_rate
