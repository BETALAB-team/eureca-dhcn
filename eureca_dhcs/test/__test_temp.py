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

time_int = 300
n_timestep = int(600)

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
    network.solve_hydraulic_balance_SIMPLE(iteration)
    network.solve_thermal_balance(iteration, time_interval=time_int)
network.save_results()

# %% Test 2
import os

from eureca_dhcs.network import Network
from eureca_dhcs.soil import Soil

# from ..network import Network
# from ..soil import Soil

time_int = 300
n_timestep = int(600)
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
    network2.solve_thermal_balance(iteration, time_interval=time_int)
network2.save_results()

#%% Result
out_path = os.path.join("eureca_dhcs", "test", "output_test_temp")
out_path2 = os.path.join("eureca_dhcs", "test", "output_test_2_temp")

import pandas as pd
import matplotlib.pyplot as plt

bp1 = pd.read_csv(os.path.join(out_path, "BranchMassFlowRates.csv"))
nt1 = pd.read_csv(os.path.join(out_path, "NodesTemperatures.csv"))
bt1 = pd.read_csv(os.path.join(out_path, "BranchTemperatures.csv"))

bp2 = pd.read_csv(os.path.join(out_path2, "BranchMassFlowRates.csv"))
nt2 = pd.read_csv(os.path.join(out_path2, "NodesTemperatures.csv"))
bt2 = pd.read_csv(os.path.join(out_path2, "BranchTemperatures.csv"))

last_t_nt1 = nt1.iloc[-1, -2]
last_t_nt2 = nt2.iloc[-1, -2]

bp1 = bp1[[col for col in bp1.columns if col.endswith("33")]].iloc[:n_timestep]
bp2 = bp2[[col for col in bp2.columns if col.endswith("44")]].iloc[:n_timestep]

nt1 = nt1[[col for col in nt1.columns if col.endswith("33")]].iloc[:n_timestep]
nt2 = nt2[[col for col in nt2.columns if col.endswith("44")]].iloc[:n_timestep]

bt1 = bt1[[col for col in bt1.columns if col.endswith("33")]].iloc[:n_timestep]
bt2 = bt2[[col for col in bt2.columns if col.endswith("44")]].iloc[:n_timestep]

# line
fig, [[ax1, ax3], [ax2, ax4]] = plt.subplots(nrows=2, ncols=2, figsize=(15, 15))
ax1.plot(
    nt1[["133", "233", "633", "1133"]].values,
    label=["supply node", "second node", "central node", "demand node"],
)
ax1.plot(
    bt1[["133", "533", "1033"]].values,
    linestyle="-.",
    label=["supply branch", "central branch", "demand branch"],
)
bt2.rename(
    {
        "344": "demand branch",
        "144": "supply branch 1",
        "244": "supply branch 2",
    },
    axis=1,
    inplace=True,
)
bp2.rename(
    {
        "344": "demand branch",
        "144": "supply branch 1",
        "244": "supply branch 2",
    },
    axis=1,
    inplace=True,
)
nt2.rename(
    {
        "144": "supply node 1",
        "244": "supply node 2",
        "344": "central node",
        "444": "demand node",
    },
    axis=1,
    inplace=True,
)
nt2.plot(
    ax=ax2,
)

bt2.plot(
    ax=ax2,
    linestyle="-.",
)

ax3.plot(
    bp1[["133", "533", "1033"]].values,
    linestyle="-.",
    label=["supply branch", "central branch", "demand branch"],
)

bp2.plot(
    ax=ax4,
    linestyle="-.",
)

ax1.set_title("One line 10 segments test")
ax2.set_title("Two branches in one test")
ax3.set_title("One line 10 segments test")
ax4.set_title("Two branches in one test")
ax1.text(
    0.7,
    0.1,
    f"Final temperature last node (after 600 time steps):\n{last_t_nt1} °C",
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax1.transAxes,
    bbox=dict(facecolor="white", alpha=0.99),
)
ax2.text(
    0.7,
    0.1,
    f"Final temperature last node (after 600 time steps):\n{last_t_nt2} °C",
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax2.transAxes,
    bbox=dict(facecolor="white", alpha=0.99),
)
# ax3.text(
#     0.65,
#     0.25,
#     f"Branches mass flow rate (after 600 time steps)",
#     horizontalalignment="center",
#     verticalalignment="center",
#     transform=ax1.transAxes,
#     bbox=dict(facecolor="white", alpha=0.99),
# )
# ax2.text(
#     0.65,
#     0.25,
#     f"Branches mass flow rate (after 600 time steps)",
#     horizontalalignment="center",
#     verticalalignment="center",
#     transform=ax2.transAxes,
#     bbox=dict(facecolor="white", alpha=0.99),
# )

for ax in [ax1, ax2]:
    # ax.legend(loc="lower left")
    ax.grid()
    ax.set_ylabel("Temperature [°C]")
    ax.set_xlabel(f"Timestep [{time_int} s]")
    ax.set_ylim([40, 52])
    ax.set_xlim([0, 600])
    ax.legend()
for ax in [ax3, ax4]:
    ax.legend(loc="upper right")
    ax.grid()
    ax.set_ylabel("Mass flow rate [kg/s]")
    ax.set_xlabel(f"Timestep [{time_int} s]")
