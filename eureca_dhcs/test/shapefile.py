import geopandas as gpd
import matplotlib.pyplot as plt
import momepy
import networkx as nx

from eureca_dhcs.network import Network

path_lines = "C:\\Users\\pratenr15640\\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\lines\\lines.shp"
path_nodes = "C:\\Users\\pratenr15640\\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\lines\\nodes.shp"

nodes_df = gpd.read_file(path_nodes)
branches_df = gpd.read_file(path_lines)

nodes_dict = {}
branches_dict = {}


for nodes_df_idx in nodes_df.index:
    node = nodes_df.loc[nodes_df_idx]
    nodes_dict[node["id"]] = {
        "id": str(node["id"]),
        "node type": str(node["node_type"]),
        "x": -0.5,
        "y": 1.0,
    }

    try:
        supply_nodes = [str(n) for n in node["supply_nod"].split(";")]
    except AttributeError:
        supply_nodes = []
    nodes_dict[node["id"]]["supply nodes"] = supply_nodes
    try:
        demand_nodes = [str(n) for n in node["demand_nod"].split(";")]
    except AttributeError:
        demand_nodes = []
    nodes_dict[node["id"]]["demand nodes"] = demand_nodes
for branch_df_idx in branches_df.index:
    branch = branches_df.loc[branch_df_idx]
    branches_dict[branch["id"]] = {
        "id": str(branch["id"]),
        "supply node": str(branch["supply_nod"]),
        "demand node": str(branch["demand_nod"]),
        "pipe diameter [m]": str(branch["pipe_d [m]"]),
    }
network = Network(nodes_dict=nodes_dict, branches_dict=branches_dict)
