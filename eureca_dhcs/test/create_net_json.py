# nodes info

import pandas as pd
import numpy as np
import os
import json

nodes = {
    "101": {
        "id": "101",
        "node type": "supply",
        "x": -0.5,
        "y": 1.0,
        "supply nodes": [],
        "demand nodes": ["121"],
    },
    "111": {
        "id": "111",
        "node type": "supply",
        "x": -1,
        "y": 0.5,
        "supply nodes": [],
        "demand nodes": ["121"],
    },
    "121": {
        "id": "121",
        "node type": "disp",
        "x": 0,
        "y": 0,
        "supply nodes": ["111"],
        "demand nodes": ["131", "141"],
    },
    "131": {
        "id": "131",
        "node type": "demand",
        "x": 1,
        "y": 0,
        "supply nodes": ["121"],
        "demand nodes": [],
    },
    "141": {
        "id": "141",
        "node type": "disp",
        "x": -0.5,
        "y": -0.5,
        "supply nodes": ["121"],
        "demand nodes": ["151"],
    },
    "151": {
        "id": "151",
        "node type": "demand",
        "x": 0,
        "y": -1,
        "supply nodes": ["151"],
        "demand nodes": [],
    },
}


branches = {
    "Z": {
        "id": "Z",
        "supply node": "101",
        "demand node": "121",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "A": {
        "id": "A",
        "supply node": "111",
        "demand node": "121",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "B": {
        "id": "B",
        "supply node": "121",
        "demand node": "131",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "C": {
        "id": "C",
        "supply node": "121",
        "demand node": "141",
        "pipe ext diameter [m]": 0.3,
        "depth [m]": 0.8,
        "pipe thickness [m]": 0.02,
        "insulation thickness [m]": 0.03,
        "pipe conductivity [W/(mK)]": 50.0,
        "insulation conductivity [W/(mK)]": 0.1,
    },
    "D": {
        "id": "D",
        "supply node": "141",
        "demand node": "151",
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
