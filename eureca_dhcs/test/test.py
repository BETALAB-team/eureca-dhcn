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


class TestNodesBranches:
    """
    This is a test class for the pytest module.
    It tests Nodes and Branch class and its property
    """

    def test_node(self):
        # Standard Node creation
        Node(
            self,
            idx="5",
            node_type="disp",
            supply_nodes=["3", "2"],
            demand_nodes=["1"],
            x=0.5,
            y=0.8,
        )
