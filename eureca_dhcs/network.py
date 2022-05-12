class Network:
    def __init__(self, nodes_dict: dict, branches_dict: dict):
        """
        Creates a district water network starting from nodes and branches dictionaries
        See the example below
        kwarg passed to each node and branch are passed to the Node and Branch constructor
        Parameters
        ----------
        nodes_file : dict
            Dict with nodes and kwargs:
                {
                "0":
                    {"id": "0",
                    "type": "supply",
                    "x": -0.5,
                    "y": 1.0,
                    "supply nodes": [],
                    "demand nodes": ["2"]
                    },
                "1":
                    {...
                     ...}
                .
                .
                .
        branches_file : dict
            Dict with branches and kwargs.
            {
                "1":
                {"id": "Z",
                "supply node": "0",
                "demand node": "2",
                "pipe diameter [m]": 0.3,
                "depth [m]": Na
                },
            "2":
                {...
                 ...}
            .
            .
            .

        Returns
        -------
        None.

        """
        # This two contains only the dictionary, with str and integer
        self._nodes_json_dict = nodes_dict
        self._branch_json_dict = branches_dict
        # These two are use to include the Node and Branch objects
        self._nodes_json_dict
        self._nodes_json_dict

    @classmethod
    def from_shapefiles(self, nodes_file: str, branches_file: str):
        """
        Creates a district water network starting from GIS nodes and branches shapefile


        Parameters
        ----------
        nodes_file : str
            Path to the nodes shapefile.
            The file must be a point shapefile with the following attributes:
                id: unique integer id
                node_type: choice from [supply, disp, demand]
                supply_nodes: string with list of nodes that supplies
                              the current node, e.g.: "1;3;6"
                demand_nodes: string with list of nodes that are supplied
                              the current node, e.g.: "2;7"
        branches_file : str
            Path to the branches shapefile.
            The file must be a line shapefile with the following attributes:
                id: unique integer id
                pipe_d [m]: diameter of the pipe
                roughness [-]: relative roughness for pressure losses
                supply_nod: id of the supply node for the branch
                demand_nod: id of the demand node for the branch


        Returns
        -------
        None.

        """
