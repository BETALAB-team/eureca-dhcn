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
network.load_boundary_conditions_from_excel(boundaries, 2)


# %%
import time

start = time.time()
for iteration in range(2):
    print(iteration)
    network.solve_hydraulic_balance_SIMPLE(iteration)
    # network.solve_hydraulic_balance(iteration)
    network.solve_thermal_balance(iteration, time_interval=60)
print(f"Sim time: {(time.time() - start):.1f} s")
network.save_results()
#%%
out_path = os.path.join("eureca_dhcs", "test", "output_tests_hydro")
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
