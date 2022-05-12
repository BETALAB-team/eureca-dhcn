import geopandas
import matplotlib.pyplot as plt
import momepy
import networkx as nx

path_lines = "C:\\Users\\pratenr15640\\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\lines\\lines.shp"
path_nodes = "C:\\Users\\pratenr15640\\OneDrive - Università degli Studi di Padova\\AttivitaDavis\\QGIS\\lines\\nodes.shp"

gdf_branches = geopandas.read_file(path_lines)
gdf_nodes = geopandas.read_file(path_lines)
