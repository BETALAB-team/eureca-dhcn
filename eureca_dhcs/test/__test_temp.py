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

path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "test_temp_1_lines.shp")
path_nodes = os.path.join(
    "eureca_dhcs", "test", "input_tests", "test_temp_1_points.shp"
)
# Boundary condition
boundaries = os.path.join(
    "eureca_dhcs", "test", "input_tests", "conditions_temp_test.xlsx"
)
# Output
out_path = os.path.join("eureca_dhcs", "test", "output_test_temp")

soil = Soil()
network = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path, temperature_mode="Cooling"
)
network.load_boundary_conditions_from_excel(boundaries, 600)

for iteration in range(600):
    network.solve_hydraulic_balance(iteration)
    network.solve_thermal_balance(iteration, time_interval=3600)
network.save_results()

# %% Test 2

path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "test_temp_2_lines.shp")
path_nodes = os.path.join(
    "eureca_dhcs", "test", "input_tests", "test_temp_2_points.shp"
)
# Boundary condition
boundaries = os.path.join(
    "eureca_dhcs", "test", "input_tests", "conditions_temp_test.xlsx"
)
# Output
out_path2 = os.path.join("eureca_dhcs", "test", "output_test_2_temp")

soil = Soil()
network2 = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path2, temperature_mode="Cooling"
)
network2.load_boundary_conditions_from_excel(boundaries, 600)

for iteration in range(600):
    network2.solve_hydraulic_balance(iteration)
    network2.solve_thermal_balance(iteration, time_interval=3600)
network2.save_results()

#%% Result

import pandas as pd
import matplotlib.pyplot as plt

out_path = os.path.join("eureca_dhcs", "test", "output_test_temp")
out_path2 = os.path.join("eureca_dhcs", "test", "output_test_2_temp")

bp1 = pd.read_csv(os.path.join(out_path, "BranchMassFlowRates.csv"))
nt1 = pd.read_csv(os.path.join(out_path, "NodesTemperatures.csv"))
bt1 = pd.read_csv(os.path.join(out_path, "BranchTemperatures.csv"))

bp2 = pd.read_csv(os.path.join(out_path2, "BranchMassFlowRates.csv"))
nt2 = pd.read_csv(os.path.join(out_path2, "NodesTemperatures.csv"))
bt2 = pd.read_csv(os.path.join(out_path2, "BranchTemperatures.csv"))

n_timestep = 100

bp1 = bp1[[col for col in bp1.columns if col.endswith("33")]].iloc[:n_timestep]
bp2 = bp2[[col for col in bp2.columns if col.endswith("44")]].iloc[:n_timestep]

nt1 = nt1[[col for col in nt1.columns if col.endswith("33")]].iloc[:n_timestep]
nt2 = nt2[[col for col in nt2.columns if col.endswith("44")]].iloc[:n_timestep]

bt1 = bt1[[col for col in bt1.columns if col.endswith("33")]].iloc[:n_timestep]
bt2 = bt2[[col for col in bt2.columns if col.endswith("44")]].iloc[:n_timestep]

# line
fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=(15, 15))
ax1.plot(
    nt1[["133", "633", "1133"]].values,
    label=["supply node", "central node", "demand node"],
)
ax1.plot(
    bt1[["133", "533", "1033"]].values,
    linestyle="-.",
    label=["supply branch", "central branch", "demand branch"],
)
ax1.legend()
ax1.grid()

ax2.plot(
    nt2.values,
    label=["supply node 1", "supply node 2", "central node", "demand node"],
)
ax2.plot(
    bt2.values,
    linestyle="-.",
    label=["supply branch 1", "supply branch 2", "demand branch"],
)
ax2.legend()
ax2.grid()
