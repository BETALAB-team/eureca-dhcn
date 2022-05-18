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
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "A": {
        "id": "A",
        "supply node": "11",
        "demand node": "12",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "B": {
        "id": "B",
        "supply node": "12",
        "demand node": "13",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "C": {
        "id": "C",
        "supply node": "12",
        "demand node": "14",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "D": {
        "id": "D",
        "supply node": "14",
        "demand node": "15",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
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
