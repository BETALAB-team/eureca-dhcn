"""
Tests
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
import json
import logging

import pytest
import numpy as np

from eureca_dhcs.node import Node
from eureca_dhcs.branch import Branch
from eureca_dhcs.network import Network
from eureca_dhcs.soil import Soil
from eureca_dhcs.exceptions import DuplicateNode, WrongNodeType


logging.basicConfig(
    filename=os.path.join(".", "eureca_dhcs", "test", "test_log.log"),
    encoding="utf-8",
    level=logging.DEBUG,
)


class TestNodesBranches:
    """
    This is a test class for the pytest module.
    It tests Nodes and Branch class and its property
    """

    def test_node(self):
        # Standard Node creation
        Node(
            idx="80000",
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
                idx="80000",
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
                idx="10000000",
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
            idx="900000",
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
            pipe_ext_diameter=0.5,  # [m]
            pipe_len=10,  # [m]
            roughness=0.2,  # [-]
            pipe_thickness=0.01,  # [m]
            pipe_depth=0.8,  # [m]
            insulation_thickness=0.01,  # [m]
            pipe_conductivity=50.0,  # [W/(m/K)]
            insulation_conductivity=0.1,  # [m]
            starting_temperature=50.0,  # [°C]
            nodes_objects_dict=None,
        )

    def test_branch_2(self):
        nodes_objects_dict = {
            "800000": Node(
                idx="800000",
                node_type="disp",
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=0.5,
                y=0.8,
            ),
            "100000000000": Node(
                idx="100000000000",
                node_type="disp",
                supply_nodes=["3", "2"],
                demand_nodes=["1"],
                x=10.5,
                y=8.8,
            ),
        }

        Branch(
            idx="A",
            supply_node="800000",
            demand_node="100000000000",
            pipe_ext_diameter=0.5,  # [m]
            roughness=0.2,  # [-]
            pipe_thickness=0.01,  # [m]
            pipe_depth=0.8,  # [m]
            insulation_thickness=0.01,  # [m]
            pipe_conductivity=50.0,  # [W/(m/K)]
            insulation_conductivity=0.1,  # [m]
            starting_temperature=50.0,  # [°C]
            nodes_objects_dict=nodes_objects_dict,
        )


class TestSoil:
    """
    This is a test class for the pytest module.
    It tests Soil class and its property
    """

    def test_soil(self):
        soil = Soil()
        t = np.zeros([365, 4])
        for c, depth in enumerate([0, 0.5, 2, 4]):
            for day in range(365):
                t[day, c] = soil.get_soil_temperature(day, depth)
        with open(
            os.path.join("eureca_dhcs", "test", "ground", "t_ground.csv"), "w"
        ) as t_ground:
            t_ground.write("0 m, 0.5 m, 2 m, 4 m\n")
            np.savetxt(t_ground, t, delimiter=",", fmt="%.2f")


class TestNetwork:
    """
    This is a test class for the pytest module.
    It tests Nodes and Branch class and its property
    """

    def test_network(self):
        # Standard Node creation
        Node._counter = 0
        Branch._counter = 0

        with open(
            os.path.join("eureca_dhcs", "test", "network_config", "nodes.json"), "r"
        ) as outfile:
            nodes = json.load(outfile)
        with open(
            os.path.join("eureca_dhcs", "test", "network_config", "branches.json"), "r"
        ) as outfile:
            branches = json.load(outfile)
        soil = Soil()
        Network(nodes_dict=nodes, branches_dict=branches, soil_obj=soil)

    def test_network_from_shape(self):
        Node._counter = 0
        Branch._counter = 0
        # Standard Node creation

        path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines.shp")
        path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes.shp")
        soil = Soil()
        network = Network.from_shapefiles(
            path_nodes,
            path_lines,
            soil,
            output_path=os.path.join("eureca_dhcs", "test", "output_tests"),
        )

        boundaries = os.path.join(
            "eureca_dhcs", "test", "input_tests", "conditions.xlsx"
        )
        network.load_boundary_conditions_from_excel(boundaries, 10)
        for iteration in range(10):
            network.solve_hydraulic_balance(iteration)
            # sol = network.solve_thermal_balance(iteration)
        network.save_results()

    def test_network_hydro(self):
        Node._counter = 0
        Branch._counter = 0
        # Standard Node creation

        path_lines = os.path.join(
            "eureca_dhcs", "test", "input_tests", "test_hydro_lines.shp"
        )
        path_nodes = os.path.join(
            "eureca_dhcs", "test", "input_tests", "test_hydro_points.shp"
        )
        soil = Soil()
        network = Network.from_shapefiles(
            path_nodes,
            path_lines,
            soil,
            output_path=os.path.join("eureca_dhcs", "test", "output_tests_hydro"),
        )

        boundaries = os.path.join(
            "eureca_dhcs", "test", "input_tests", "conditions_hydro_test.xlsx"
        )
        network.load_boundary_conditions_from_excel(boundaries, 2)
        for iteration in range(2):
            network.solve_hydraulic_balance(iteration)
            # sol = network.solve_thermal_balance(iteration)
        network.save_results()
