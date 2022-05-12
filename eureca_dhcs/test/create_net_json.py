# nodes info

import pandas as pd
import numpy as np
import os
import json

nodes = {
    "10": {
        "id": "10",
        "node type": "supply",
        "x": -0.5,
        "y": 1.0,
        "supply nodes": [],
        "demand nodes": ["12"],
    },
    "11": {
        "id": "11",
        "node type": "supply",
        "x": -1,
        "y": 0.5,
        "supply nodes": [],
        "demand nodes": ["12"],
    },
    "12": {
        "id": "12",
        "node type": "disp",
        "x": 0,
        "y": 0,
        "supply nodes": ["11"],
        "demand nodes": ["13", "14"],
    },
    "13": {
        "id": "13",
        "node type": "demand",
        "x": 1,
        "y": 0,
        "supply nodes": ["12"],
        "demand nodes": [],
    },
    "14": {
        "id": "14",
        "node type": "disp",
        "x": -0.5,
        "y": -0.5,
        "supply nodes": ["12"],
        "demand nodes": ["15"],
    },
    "15": {
        "id": "15",
        "node type": "demand",
        "x": 0,
        "y": -1,
        "supply nodes": ["15"],
        "demand nodes": [],
    },
}


branches = {
    "Z": {
        "id": "Z",
        "supply node": "10",
        "demand node": "12",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "A": {
        "id": "A",
        "supply node": "11",
        "demand node": "12",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "B": {
        "id": "B",
        "supply node": "12",
        "demand node": "13",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "C": {
        "id": "C",
        "supply node": "12",
        "demand node": "14",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "D": {
        "id": "D",
        "supply node": "14",
        "demand node": "15",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
}

with open(
    os.path.join("eureca_dhcs", "test", "network_config", "nodes.json"), "w"
) as outfile:
    json.dump(nodes, outfile)
with open(
    os.path.join("eureca_dhcs", "test", "network_config", "branches.json"), "w"
) as outfile:
    json.dump(branches, outfile)
