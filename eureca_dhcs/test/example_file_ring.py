"""
This is an example file to use the network class and solve water DHC networks
"""
__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import os
import time

from eureca_dhcs.network import Network
from eureca_dhcs.soil import Soil

sim = "ring"
if sim == "ring":
    path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines_ring.shp")
    path_nodes = os.path.join(
        "eureca_dhcs", "test", "input_tests", "nodes_ring_point.shp"
    )
    # Boundary condition
    boundaries = os.path.join(
        "eureca_dhcs", "test", "input_tests", "conditions_ring.xlsx"
    )
    soil = Soil()
    network = Network.from_shapefiles(
        path_nodes, path_lines, soil, output_path="output_ring"
    )
    network.load_boundary_conditions_from_excel(boundaries, 8760)
else:
    path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines.shp")
    path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes.shp")
    # Boundary condition
    boundaries = os.path.join("eureca_dhcs", "test", "input_tests", "conditions.xlsx")
    soil = Soil()
    network = Network.from_shapefiles(
        path_nodes, path_lines, soil, output_path="output"
    )
    network.load_boundary_conditions_from_excel(boundaries, 8760)
# %%
start = time.time()
for iteration in range(8760):
    xh, qh, x0h = network.solve_hydraulic_balance(iteration)

    xt, At, qt = network.solve_thermal_balance(iteration)
    # print(f"##### timestep {iteration} ######")
    # for b_k, b in network._branches_object_dict.items():
    #     print(
    #         f"Branch {b_k}: mass {b._mass_flow_rate:.1f}\t\t{b._roughness}\t\tTemp {b._branch_temperature:.3f}"
    #     )
    # for n_k, n in network._nodes_object_dict.items():
    #     print(
    #         f"Node {n_k}: pressure {n._node_pressure:6.1f}\t\t\tTemp {n._node_temperature:.3f}"
    #     )
network.save_results()
print(f"{time.time()-start:.1f} s")

#%%
path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "test_hydro_lines.shp")
path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "test_hydro_points.shp")
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
network.load_boundary_conditions_from_excel(boundaries, 3)
for iteration in range(3):
    network.solve_hydraulic_balance(iteration)
    # sol = network.solve_thermal_balance(iteration)
