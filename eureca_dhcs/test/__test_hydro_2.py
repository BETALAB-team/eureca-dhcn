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
    "eureca_dhcs", "test", "input_tests", "conditions_hydro_test_stability.xlsx"
)
# Output
out_path = os.path.join("eureca_dhcs", "test", "output_test_hydro_2")

soil = Soil()
network = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path, temperature_mode="Cooling"
)
network.load_boundary_conditions_from_excel(boundaries, 135)

# %%
for iteration in range(135):
    network.solve_hydraulic_balance(iteration)
    network.solve_thermal_balance(iteration, time_interval=60)
network.save_results()
print("Branch 2: ")
print("\tExercise: 40.2")
print("\t", network._branches_object_dict["2"]._mass_flow_rate_array[102])
print("\t", network._branches_object_dict["2"]._mass_flow_rate_array[116])
print("\t", network._branches_object_dict["2"]._mass_flow_rate_array[131])

print("Branch 7: ")
print("\tExercise: 15.0")
print("\t", network._branches_object_dict["7"]._mass_flow_rate_array[102])
print("\t", network._branches_object_dict["7"]._mass_flow_rate_array[116])
print("\t", network._branches_object_dict["7"]._mass_flow_rate_array[131])

#%%
out_path = os.path.join("eureca_dhcs", "test", "output_test_hydro_2")
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
