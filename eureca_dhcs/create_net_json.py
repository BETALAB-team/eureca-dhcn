# nodes info

import pandas as pd
import numpy as np
import os
import json

nodes = {
    "0": {
        "id": "0",
        "node_type": "supply",
        "x": -0.5,
        "y": 1.0,
        "supply nodes": [],
        "demand nodes": ["2"],
    },
    "1": {
        "id": "1",
        "node_type": "supply",
        "x": -1,
        "y": 0.5,
        "supply nodes": [],
        "demand nodes": ["2"],
    },
    "2": {
        "id": "2",
        "node_type": "disp",
        "x": 0,
        "y": 0,
        "supply nodes": ["1"],
        "demand nodes": ["3", "4"],
    },
    "3": {
        "id": "3",
        "node_type": "demand",
        "x": 1,
        "y": 0,
        "supply nodes": ["2"],
        "demand nodes": [],
    },
    "4": {
        "id": "4",
        "node_type": "disp",
        "x": -0.5,
        "y": -0.5,
        "supply nodes": ["2"],
        "demand nodes": ["5"],
    },
    "5": {
        "id": "5",
        "node_type": "demand",
        "x": 0,
        "y": -1,
        "supply nodes": ["5"],
        "demand nodes": [],
    },
}


branches = {
    "Z": {
        "id": "Z",
        "supply node": "0",
        "demand node": "2",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "A": {
        "id": "A",
        "supply node": "1",
        "demand node": "2",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "B": {
        "id": "B",
        "supply node": "2",
        "demand node": "3",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "C": {
        "id": "C",
        "supply node": "2",
        "demand node": "4",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
    "D": {
        "id": "D",
        "supply node": "4",
        "demand node": "5",
        "pipe diameter [m]": 0.3,
        "depth [m]": np.nan,
    },
}

with open(os.path.join("network_config", "nodes.json"), "w") as outfile:
    json.dump(nodes, outfile)
with open(os.path.join("network_config", "branches.json"), "w") as outfile:
    json.dump(branches, outfile)
