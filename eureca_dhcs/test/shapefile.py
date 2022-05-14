from IPython import get_ipython

get_ipython().magic("reset -sf")

import geopandas as gpd
import matplotlib.pyplot as plt
import momepy
import networkx as nx
import os

from eureca_dhcs.network import Network
from eureca_dhcs.node import Node
from eureca_dhcs.branch import Branch


path_lines = os.path.join("eureca_dhcs", "test", "input_tests", "lines.shp")
path_nodes = os.path.join("eureca_dhcs", "test", "input_tests", "nodes.shp")

# nodes_df = gpd.read_file(path_nodes)
# branches_df = gpd.read_file(path_lines)

# for nodes_df_idx in nodes_df.index:
#     node = nodes_df.loc[nodes_df_idx]
#     print(node["geometry"].x)

network = Network.from_shapefiles(path_nodes, path_lines)
# %%
x = network.solve_hydraulic_balance()
