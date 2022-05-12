"""
Tests
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os

import pytest
import numpy as np

from eureca_dhcs.node import Node
from eureca_dhcs.branch import Branch
from eureca_dhcs.exceptions import DuplicateNode, WrongNodeType


class TestNodesBranches:
    """
    This is a test class for the pytest module.
    It tests Nodes and Branch class and its property
    """

    def test_node(self):
        # Standard Node creation
        Node(
            idx="5",
            node_type="disp",
            supply_nodes=["3", "2"],
            demand_nodes=["1"],
            x=0.5,
            y=0.8,
        )

    def test_node_2(self):
        # Standard Node creation
        with pytest.raises(DuplicateNode):
            Node(
                idx="5",
                node_type="disp",
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=0.5,
                y=0.8,
            )

    def test_node_3(self):
        # Standard Node creation
        with pytest.raises(WrongNodeType):
            Node(
                idx="8",
                node_type=True,
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=0.5,
                y=0.8,
            )

    def test_node_4(self):
        # Standard Node creation
        # with pytest.raises(TypeError):
        Node(
            idx="9",
            node_type="disp",
            supply_nodes=["3", "2"],
            demand_nodes=["1"],
            x="A",
            y=0.8,
        )

    def test_branch(self):

        Branch(
            idx="A",
            supply_node="3",
            demand_node="4",
            pipe_diameter=0.5,  # [m]
            pipe_len=10,  # [m]
            roughness=0.2,  # [-]
            starting_temp=50.0,  # [°C]
            nodes_objects_dict=None,
        )

    def test_branch_2(self):
        nodes_objects_dict = {
            "8": Node(
                idx="8",
                node_type="disp",
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=0.5,
                y=0.8,
            ),
            "10": Node(
                idx="10",
                node_type="disp",
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=10.5,
                y=8.8,
            ),
        }

        Branch(
            idx="A",
            supply_node="8",
            demand_node="10",
            pipe_diameter=0.5,  # [m]
            roughness=0.2,  # [-]
            starting_temp=50.0,  # [°C]
            nodes_objects_dict=nodes_objects_dict,
        )
