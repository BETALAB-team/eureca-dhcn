__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
import logging

import numpy as np
import geopandas as gpd

from eureca_dhcs.node import Node
from eureca_dhcs.branch import Branch
from eureca_dhcs.exceptions import EmptyNetworkNodes, DuplicateNode, WrongNodeType


class Network:
    def __init__(self, nodes_dict: dict, branches_dict: dict):
        """
        Creates a district water network starting from nodes and branches dictionaries
        See the example below
        kwarg passed to each node and branch are passed to the Node and Branch constructor
        Parameters
        ----------
        nodes_file : dict
            Dict with nodes and kwargs:
                {
                "0":
                    {"id": "0",
                    "type": "supply",
                    "x": -0.5,
                    "y": 1.0,
                    "supply nodes": [],
                    "demand nodes": ["2"]
                    },
                "1":
                    {...
                     ...}
                .
                .
                .
        branches_file : dict
            Dict with branches and kwargs.
            {
                "1":
                {"id": "Z",
                "supply node": "0",
                "demand node": "2",
                "pipe diameter [m]": 0.3,
                "depth [m]": Na
                },
            "2":
                {...
                 ...}
            .
            .
            .

        Returns
        -------
        None.

        """
        # This two contains only the dictionary, with str and integer
        self._nodes_json_dict = nodes_dict
        self._branch_json_dict = branches_dict
        # These two are use to include the Node and Branch objects
        self._nodes_object_dict = {}
        self._branches_object_dict = {}

        # Creation of the nodes objects... THIS MUST BE DONE BEFORE THE CREATION OF BRANCHES
        self._create_nodes()
        # Creation of the branch objects... THIS MUST BE DONE BEFORE THE CREATION OF BRANCHES
        self._create_branches()
        # Couples nodes and branches objects
        self._couple_nodes_to_nodes()
        # DO IN ORDER
        self._couple_branches_to_nodes()

    def _create_nodes(self):
        """
        Creation and filling of the nodes object dictionary

        Returns
        -------
        None.

        """
        for idx, node in self._nodes_json_dict.items():
            self._nodes_object_dict[node["id"]] = Node(
                idx=node["id"],
                node_type=node["node type"],
                supply_nodes=node["supply nodes"],
                demand_nodes=node["demand nodes"],
                x=node["x"],
                y=node["y"],
            )

    def _create_branches(self):
        """
        Creation and filling of the nodes object dictionary

        Returns
        -------
        None.

        """
        if not self._nodes_object_dict:
            try:
                self._create_nodes()
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        for idx, branch in self._branch_json_dict.items():
            self._branches_object_dict[branch["id"]] = Branch(
                idx=branch["id"],
                supply_node=branch["supply node"],
                demand_node=branch["demand node"],
                pipe_diameter=branch["pipe diameter [m]"],
                nodes_objects_dict=self._nodes_object_dict,
            )

            if "pipe length [m]" in branch.keys():
                self._branches_object_dict[branch["id"]]._pipe_len = branch[
                    "pipe length [m]"
                ]
            if "roughness [-]" in branch.keys():
                self._branches_object_dict[branch["id"]]._roughness = branch[
                    "roughness [-]"
                ]
            if "starting temperature [°C]" in branch.keys():
                self._branches_object_dict[branch["id"]]._previous_temp = branch[
                    "starting temperature [°C]"
                ]

    def _couple_nodes_to_nodes(self):
        """
        Couples nodes to supply and demand nodes

        Returns
        -------
        None.

        """
        if not self._nodes_object_dict:
            try:
                self._create_nodes()
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        if not self._branches_object_dict:
            try:
                self._create_branches()
                logging.warning(
                    "Network creation. Trying to create the branches dictionary before coupling. To create couple nodes and branches you must have already create the branches' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create couple nodes and branches you must have already create the branches' dict"
                )
        for node_idx, node in self._nodes_object_dict.items():
            node._supply_nodes_objects = {
                s_node: self._nodes_object_dict[s_node]
                for s_node in node._supply_nodes_list
            }
            node._demand_nodes_objects = {
                s_node: self._nodes_object_dict[s_node]
                for s_node in node._demand_nodes_list
            }

    def _couple_branches_to_nodes(self):
        """
        Creation and filling of the nodes object dictionary

        Returns
        -------
        None.

        """
        if not self._nodes_object_dict:
            try:
                self._create_nodes()
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        if not self._branches_object_dict:
            try:
                self._create_branches()
                logging.warning(
                    "Network creation. Trying to create the branches dictionary before coupling. To create couple nodes and branches you must have already create the branches' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create couple nodes and branches you must have already create the branches' dict"
                )
        for branch_idx, branch in self._branches_object_dict.items():
            branch._supply_node_object = self._nodes_object_dict[
                branch._supply_node_idx
            ]
            branch._demand_node_object = self._nodes_object_dict[
                branch._demand_node_idx
            ]

            # Adding branch to the node objects
            # Add to the supply node the branch (it is considered as demand branch for the node)
            self._nodes_object_dict[branch._supply_node_idx]._demand_branches_objects[
                branch_idx
            ] = branch
            # Add to the demand node the branch (it is considered as supply branch for the node)
            self._nodes_object_dict[branch._demand_node_idx]._supply_branches_objects[
                branch_idx
            ] = branch

    @classmethod
    def from_shapefiles(cls, nodes_file: str, branches_file: str):
        """
        Creates a district water network starting from GIS nodes and branches shapefile


        Parameters
        ----------
        nodes_file : str
            Path to the nodes shapefile.
            The file must be a point shapefile with the following attributes:
                id: unique integer id
                node_type: choice from [supply, disp, demand]
                supply_nodes: string with list of nodes that supplies
                              the current node, e.g.: "1;3;6"
                demand_nodes: string with list of nodes that are supplied
                              the current node, e.g.: "2;7"
        branches_file : str
            Path to the branches shapefile.
            The file must be a line shapefile with the following attributes:
                id: unique integer id
                pipe_d [m]: diameter of the pipe
                roughness [-]: relative roughness for pressure losses
                supply_nod: id of the supply node for the branch
                demand_nod: id of the demand node for the branch


        Returns
        -------
        None.

        """

        nodes_df = gpd.read_file(nodes_file)
        branches_df = gpd.read_file(branches_file)

        nodes_dict = {}
        branches_dict = {}

        for nodes_df_idx in nodes_df.index:
            node = nodes_df.loc[nodes_df_idx]
            nodes_dict[node["id"]] = {
                "id": str(node["id"]),
                "node type": str(node["node_type"]),
                "x": -0.5,
                "y": 1.0,
            }

            try:
                supply_nodes = [str(n) for n in node["supply_nod"].split(";")]
            except AttributeError:
                supply_nodes = []
            nodes_dict[node["id"]]["supply nodes"] = supply_nodes
            try:
                demand_nodes = [str(n) for n in node["demand_nod"].split(";")]
            except AttributeError:
                demand_nodes = []
            nodes_dict[node["id"]]["demand nodes"] = demand_nodes
        for branch_df_idx in branches_df.index:
            branch = branches_df.loc[branch_df_idx]
            branches_dict[branch["id"]] = {
                "id": str(branch["id"]),
                "supply node": str(branch["supply_nod"]),
                "demand node": str(branch["demand_nod"]),
                "pipe diameter [m]": str(branch["pipe_d [m]"]),
            }
        return cls(nodes_dict=nodes_dict, branches_dict=branches_dict)
