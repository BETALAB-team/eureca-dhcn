__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math

from eureca_dhcs.exceptions import DuplicateNode, WrongNodeType


class Node:
    """
    Class for the network Node
    """

    _idx_list = []
    _counter = 0

    def __init__(
        self,
        idx: str,
        node_type: str,
        supply_nodes: list,
        demand_nodes: list,
        x: float,
        y: float,
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
