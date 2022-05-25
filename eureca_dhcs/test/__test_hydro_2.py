"""
This is an example file to test the hydraulic network number 2
"""
__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


import os

from eureca_dhcs.network import Network
from eureca_dhcs.soil import Soil

path_lines = os.path.join(
    "eureca_dhcs", "test", "input_tests", "test_hydro_lines_2.shp"
)
path_nodes = os.path.join(
    "eureca_dhcs", "test", "input_tests", "test_hydro_points_2.shp"
)
# Boundary condition
boundaries = os.path.join(
    "eureca_dhcs", "test", "input_tests", "conditions_hydro_test_2.xlsx"
)
# Output
out_path = os.path.join("eureca_dhcs", "test", "output_test_hydro_2")

soil = Soil()
network = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path, temperature_mode="Cooling"
)
network.load_boundary_conditions_from_excel(boundaries, 51)

# %%
for iteration in range(51):
    network.solve_hydraulic_balance(iteration)
    # network.solve_thermal_balance(iteration, time_interval=60)
network.save_hydraulic_results()
print(network._branches_object_dict["2"]._mass_flow_rate_array.transpose())
