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
from eureca_dhcs._thermal_system_function import (
    thermal_balance_system_optimization,
    thermal_balance_system_inverse,
)

path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines.shp")
path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes.shp")
soil = Soil()
network = Network.from_shapefiles(
    path_nodes,
    path_lines,
    soil,
    output_path=os.path.join("eureca_dhcs", "test", "output_tests"),
)

boundaries = os.path.join("eureca_dhcs", "test", "input_tests", "conditions.xlsx")
network.load_boundary_conditions_from_excel(boundaries, 30)


#%%
for iteration in range(30):
    network.solve_hydraulic_balance(iteration)
    sol = network.solve_thermal_balance(iteration)
network.save_hydraulic_results()
for branch in network._branches_object_ordered_list:
    print(f"Branch {branch._idx}")
    print(f"\text diameter [m]: {branch._pipe_ext_diameter}")
    print(f"\tint diameter [m]: {branch._pipe_int_diameter}")
# Risolvere problema roughness!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # %%
# for iteration in range(20):
#     x = network.solve_hydraulic_balance(iteration)

out_path = os.path.join("eureca_dhcs", "test", "output_tests")
import pandas as pd
import matplotlib.pyplot as plt

bp1 = pd.read_csv(os.path.join(out_path, "BranchMassFlowRates.csv"))
# nt1 = pd.read_csv(os.path.join(out_path, "NodesTemperatures.csv"))
# bt1 = pd.read_csv(os.path.join(out_path, "BranchTemperatures.csv"))

# bp1 = bp1[[col for col in bp1.columns if col.endswith("33")]]
# nt1 = nt1[[col for col in nt1.columns if col.endswith("33")]]
# bt1 = bt1[[col for col in bt1.columns if col.endswith("33")]]

fig, ax4 = plt.subplots(nrows=1, ncols=1, figsize=(15, 5))
ax4.plot(
    bp1.values,
    linestyle="-",
)

ax4.set_title("HydroTest")
ax4.legend(loc="lower right")
ax4.grid()
ax4.set_ylabel("Mass flow rate [kg/s]")
