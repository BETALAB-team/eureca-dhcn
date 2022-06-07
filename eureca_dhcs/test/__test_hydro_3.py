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

# from ..network import Network
# from ..soil import Soil

time_int = 3600

path_nodes = "C:\\Users\\pratenr15640\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\dwg_network\\Eureca_CWNorthQuadPoints__.shp"
path_lines = "C:\\Users\\pratenr15640\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\dwg_network\\Eureca_CWNorthQuadLines_.shp"

# Boundary condition
boundaries = os.path.join("eureca_dhcs", "test", "input_tests", "ConditionsWithP.xlsx")
# Output
out_path = os.path.join("eureca_dhcs", "test", "output_test_with_P")

soil = Soil()
network = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path, temperature_mode="Cooling"
)
network.load_boundary_conditions_from_excel(boundaries, 330)

for iteration in range(330):
    print(iteration)
    # x, q, x0 = network.solve_hydraulic_balance(iteration)
    # network.solve_thermal_balance(iteration, time_interval=time_int)
network.save_results()


#%% Result
out_path = os.path.join("eureca_dhcs", "test", "output_test_with_P")

import pandas as pd
import matplotlib.pyplot as plt

bp1 = pd.read_csv(os.path.join(out_path, "BranchMassFlowRates.csv"))
nt1 = pd.read_csv(os.path.join(out_path, "NodesTemperatures.csv"))
bt1 = pd.read_csv(os.path.join(out_path, "BranchTemperatures.csv"))

fig, ax1 = plt.subplots(figsize=(15, 15))

bp1[["94", "95", "106"]].plot(ax=ax1)
for ax in [ax1]:
    ax.legend(loc="upper right")
    ax.grid()
    ax.set_ylabel("Mass flow rate [kg/s]")
    ax.set_xlabel(f"Timestep [{time_int} s]")
