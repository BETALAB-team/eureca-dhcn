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

# Boundary conditions
boundaries = os.path.join(".", "eureca_dhcs", "test", "CWN", "ConditionsWithMFRT.xlsx")
# Output
out_path = os.path.join(".", "eureca_dhcs", "test", "CWN", "results")

soil = Soil()
network = Network.from_shapefiles(
    path_nodes, path_lines, soil, output_path=out_path, temperature_mode="Cooling"
)
network.load_boundary_conditions_from_excel(boundaries, 330)

for iteration in range(330):
    print(iteration)
    x, q, x0 = network.solve_hydraulic_balance(iteration)
    # network.solve_thermal_balance(iteration, time_interval=time_int)
network.save_results()
