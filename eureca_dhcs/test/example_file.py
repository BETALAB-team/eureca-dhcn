"""
This is an example file to use the network class and solve water DHC networks
"""
__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import os

from eureca_dhcs.network import Network
from eureca_dhcs.soil import Soil

path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines.shp")
path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes.shp")
# Boundary condition
boundaries = os.path.join("eureca_dhcs", "test", "input_tests", "conditions.xlsx")
soil = Soil()
network = Network.from_shapefiles(path_nodes, path_lines, soil, output_path="output")
network.load_boundary_conditions_from_excel(boundaries, 100)


# Risolvere problema roughness!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# %%
for iteration in range(10):
    x = network.solve_hydraulic_balance(iteration)
network.save_hydraulic_results()
