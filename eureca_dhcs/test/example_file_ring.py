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

path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines_ring.shp")
path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes_ring_point.shp")
# Boundary condition
boundaries = os.path.join("eureca_dhcs", "test", "input_tests", "conditions_ring.xlsx")

network = Network.from_shapefiles(path_nodes, path_lines, output_path="output_ring")
network.load_boundary_conditions_from_excel(boundaries, 100)

# %%
for iteration in range(20):
    x = network.solve_hydraulic_balance(iteration)
    print(f"##### timestep {iteration} ######")
    for b_k, b in network._branches_object_dict.items():
        print(f"Branch {b_k}: mass {b._mass_flow_rate}\t\t{b._roughness}")
    for n_k, n in network._nodes_object_dict.items():
        print(f"Node {n_k}: pressure {n._node_pressure}")
network.save_hydraulic_results()
