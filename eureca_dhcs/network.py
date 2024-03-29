__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os
import logging

import numpy as np
import geopandas as gpd
import pandas as pd
from scipy.optimize import root

from eureca_dhcs.node import Node
from eureca_dhcs.branch import Branch
from eureca_dhcs.soil import Soil
from eureca_dhcs._hydraulic_system_function import hydraulic_balance_system
from eureca_dhcs._thermal_system_function import (
    thermal_balance_system_optimization,
    thermal_balance_system_inverse,
)
from eureca_dhcs.exceptions import (
    EmptyNetworkNodes,
    DuplicateNode,
    WrongNodeType,
    BoundaryConditionNotProvided,
)


class Network:
    def __init__(
        self,
        nodes_dict: dict,
        branches_dict: dict,
        soil_obj: Soil,
        output_path=None,
        temperature_mode="Heating",
    ):
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
                "id": "D",
                "supply node": "14",
                "demand node": "15",
                "pipe ext diameter [m]": 0.3,
                "depth [m]": 0.8,
                "pipe thickness [m]": 0.02,
                "insulation thickness [m]": 0.03,
                "pipe conductivity [W/(mK)]": 50.0,
                "insulation conductivity [W/(mK)]": 0.1,
            "2":
                {...
                 ...}
            .
            .
            .

        soil_obj: Soil
            Soil object to associate to the network (with soil properties)

        output_path: str
            path to where results are save. If it does not exist the tool creates it

        temperature_mode: str
            Choose between Heating and Cooling.

        Returns
        -------
        None.

        """
        # Standard Node creation
        Node._counter = 0
        Branch._counter = 0
        # This two contains only the dictionary, with str and integer
        self._nodes_json_dict = nodes_dict
        self._branch_json_dict = branches_dict
        # These two are use to include the Node and Branch objects
        self._nodes_object_dict = {}
        self._branches_object_dict = {}
        # Needed to vectorize hydraulic resistances
        self._branches_object_ordered_list = []
        self._nodes_object_ordered_list = []
        self._soil_obj = soil_obj
        # Creation of the nodes objects... THIS MUST BE DONE BEFORE THE CREATION OF BRANCHES
        self._create_nodes(temperature_mode)
        # Creation of the branch objects... THIS MUST BE DONE BEFORE THE CREATION OF BRANCHES
        self._create_branches(temperature_mode)
        # Couples nodes and branches objects
        self._couple_nodes_to_nodes()
        # DO IN ORDER
        self._couple_branches_to_nodes()
        # function to calc adjacency matrix
        self._calc_adjacency_matrix()

        # Output file
        self.output_path = None
        if output_path != None:
            self.output_path = str(output_path)
            self._create_output_folder()

    @property
    def _soil_obj(self) -> Soil:
        return self.__soil_obj

    @_soil_obj.setter
    def _soil_obj(self, value: Soil):
        if not isinstance(value, Soil):
            raise TypeError(
                f"Network soil object, a soil object must be passed. Type: {type(value)}"
            )
        self.__soil_obj = value

    def _create_output_folder(self):
        if self.output_path == None:
            raise AttributeError(f"To create output file you must set an output folder")
        if not os.path.isdir(self.output_path):
            os.mkdir(self.output_path)
        # Logging file
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)  # or whatever
        handler = logging.FileHandler(
            os.path.join(self.output_path, "logging.log"), "w", "utf-8"
        )  # or whatever
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )  # or whatever
        handler.setFormatter(formatter)  # Pass handler as a parameter, not assign
        root_logger.addHandler(handler)

    def _create_nodes(self, temperature_mode):
        """
        Creation and filling of the nodes object dictionary

        Returns
        -------
        None.

        """
        for idx, node in self._nodes_json_dict.items():
            self._nodes_object_dict[node["id"]] = Node(
                idx=node["id"],
                node_type=node["node type"],
                supply_nodes=node["supply nodes"],
                demand_nodes=node["demand nodes"],
                x=node["x"],
                y=node["y"],
                temperature_mode=temperature_mode,
            )
            self._nodes_object_ordered_list.append(self._nodes_object_dict[node["id"]])
        self._nodes_number = int(len(self._nodes_object_dict.keys()))

    def _create_branches(self, temperature_mode):
        """
        Creation and filling of the nodes object dictionary


        Parameters
        ----------
        temperature_mode : str
            Choose between Heating and Cooling.

        Raises
        ------
        EmptyNetworkNodes
            raises if the network has no nodes.

        Returns
        -------
        None.

        """

        if not self._nodes_object_dict:
            try:
                self._create_nodes(temperature_mode)
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        for idx, branch in self._branch_json_dict.items():
            self._branches_object_dict[branch["id"]] = Branch(
                idx=branch["id"],
                supply_node=branch["supply node"],
                demand_node=branch["demand node"],
                pipe_ext_diameter=branch["pipe ext diameter [m]"],
                pipe_thickness=branch["pipe thickness [m]"],
                pipe_depth=branch["depth [m]"],
                insulation_thickness=branch["insulation thickness [m]"],
                pipe_conductivity=branch["pipe conductivity [W/(mK)]"],
                insulation_conductivity=branch["insulation conductivity [W/(mK)]"],
                nodes_objects_dict=self._nodes_object_dict,
                soil_obj=self._soil_obj,
                temperature_mode=temperature_mode,
            )

            if "pipe length [m]" in branch.keys():
                self._branches_object_dict[branch["id"]]._pipe_len = branch[
                    "pipe length [m]"
                ]
            if "roughness [-]" in branch.keys():
                self._branches_object_dict[branch["id"]]._roughness = branch[
                    "roughness [-]"
                ]
            if "starting temperature [°C]" in branch.keys():
                self._branches_object_dict[branch["id"]]._previous_temp = branch[
                    "starting temperature [°C]"
                ]
            self._branches_object_ordered_list.append(
                self._branches_object_dict[branch["id"]]
            )
        self._branches_number = int(len(self._branches_object_dict.keys()))

    def _couple_nodes_to_nodes(self):
        """
        Couples nodes to supply and demand nodes

        Returns
        -------
        None.

        """
        if not self._nodes_object_dict:
            try:
                self._create_nodes()
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        if not self._branches_object_dict:
            try:
                self._create_branches()
                logging.warning(
                    "Network creation. Trying to create the branches dictionary before coupling. To create couple nodes and branches you must have already create the branches' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create couple nodes and branches you must have already create the branches' dict"
                )
        for node_idx, node in self._nodes_object_dict.items():
            node._supply_nodes_objects = {
                s_node: self._nodes_object_dict[s_node]
                for s_node in node._supply_nodes_list
            }
            node._demand_nodes_objects = {
                s_node: self._nodes_object_dict[s_node]
                for s_node in node._demand_nodes_list
            }

    def _couple_branches_to_nodes(self):
        """
        Creation and filling of the nodes object dictionary

        Returns
        -------
        None.

        """
        if not self._nodes_object_dict:
            try:
                self._create_nodes()
                logging.warning(
                    "Network creation. Trying to create the nodes dictionary before the branches dictionary. To create branches dict you must have already create the nodes' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create branches dict you must have already create the nodes' dict"
                )
        if not self._branches_object_dict:
            try:
                self._create_branches()
                logging.warning(
                    "Network creation. Trying to create the branches dictionary before coupling. To create couple nodes and branches you must have already create the branches' dict"
                )
            except Exception:
                raise EmptyNetworkNodes(
                    f"Network creation. To create couple nodes and branches you must have already create the branches' dict"
                )
        for branch_idx, branch in self._branches_object_dict.items():
            branch._supply_node_object = self._nodes_object_dict[
                branch._supply_node_idx
            ]
            branch._demand_node_object = self._nodes_object_dict[
                branch._demand_node_idx
            ]

            # Adding branch to the node objects
            # Add to the supply node the branch (it is considered as demand branch for the node)
            self._nodes_object_dict[branch._supply_node_idx]._demand_branches_objects[
                branch_idx
            ] = branch
            # Add to the demand node the branch (it is considered as supply branch for the node)
            self._nodes_object_dict[branch._demand_node_idx]._supply_branches_objects[
                branch_idx
            ] = branch

    def _calc_adjacency_matrix(self):
        # calc the coupling nodes/branches matrix
        # Matrix -> n_nodes x n_branches
        connection_matrix = np.zeros([self._nodes_number, self._branches_number])
        for branch_id, branch in self._branches_object_dict.items():
            column = branch._unique_matrix_idx
            connection_matrix[
                branch._supply_node_object._unique_matrix_idx, column
            ] = -1
            connection_matrix[branch._demand_node_object._unique_matrix_idx, column] = 1
        self._adjacency_matrix = connection_matrix

    def load_boundary_conditions_from_excel(
        self,
        excel_path: str,
        number_of_timesteps: int,
    ):
        boundaries = pd.read_excel(
            excel_path,
            sheet_name=["Hydraulic", "Thermal"],
            index_col=0,
            header=[0, 1, 2],
            parse_dates=True,
        )

        nodes = boundaries["Hydraulic"]["Node"]["Mass flow rate [kg/s]"]
        nodes_dict = nodes.to_dict(orient="List")
        nodes_dict = {str(k): np.array(node) for k, node in nodes_dict.items()}
        try:
            branches = boundaries["Hydraulic"]["Branch"]["Pump pressure raise [Pa]"]
            branches_dict = branches.to_dict(orient="List")
            # Just to convert in str and np.array
            branches_dict = {
                str(k): np.array(branch) for k, branch in branches_dict.items()
            }
        except KeyError:
            branches_dict = {}
        self.load_hydraulic_boundary_conditions(
            nodes_dict, branches_dict, number_of_timesteps
        )

        nodes = boundaries["Thermal"]["Node"]["Temperature [°C]"]
        nodes_dict = nodes.to_dict(orient="List")
        nodes_dict = {str(k): np.array(node) for k, node in nodes_dict.items()}
        self.load_thermal_boundary_conditions(nodes_dict, number_of_timesteps)

        # DateTime Index
        self.set_timestep_serie(boundaries["Hydraulic"].index)

    def set_timestep_serie(self, series: pd.Series):
        """
        This method sets the timesteps array:
        Please load a pandas series with datetime object.
        It is necessary for ground calculation temeprature

        Parameters
        ----------
        series: pd.Series

        Returns
        -------
        None.
        """
        if not isinstance(series, pd.core.indexes.datetimes.DatetimeIndex):
            raise TypeError(
                f"You must provide a pandas Datetime index! Type provided: {type(series)}"
            )
        self.timestep_array = series

    def load_thermal_boundary_conditions(
        self,
        nodes_boundary_conditions: dict,
        number_of_timesteps: int,
    ):
        """
        This method execute:
            self.load_nodes_temperature_boundary_condition(nodes_boundary_conditions, number_of_timesteps)

        Parameters
        ----------
        nodes_boundary_conditions : dict
            Dictionary with the following sintax (temperature of the supply nodes):
                "1" : np.array([50,55,55,....]),
                "3" : np.array([70,75,78,....]),
                .
                .
                ..

        number_of_timesteps : int
            Number of timesteps to be considered.

        Returns
        -------
        None.
        """
        self.load_nodes_temperature_boundary_condition(
            nodes_boundary_conditions, number_of_timesteps
        )

    def load_nodes_temperature_boundary_condition(
        self, boundary_conditions: dict, number_of_timesteps: int
    ):
        """
        This method loads the temperature boundary conditions of the nodes.
        Use °C as reference unit.

        Parameters
        ----------
        nodes_boundary_conditions : dict
            Dictionary with the following sintax (temperature of the supply nodes):
                "1" : np.array([50,55,55,....]),
                "3" : np.array([70,75,78,....]),
                .
                .
                ..

        number_of_timesteps : int
            Number of timesteps to be considered.

        Returns
        -------
        None.
        """

        for node_k, node in self._nodes_object_dict.items():
            if node._node_type in ["supply"]:
                try:
                    if len(boundary_conditions[node_k]) < number_of_timesteps:
                        raise ValueError(
                            f"Node {node._idx}: the boundary condition for the node is shorter than the number of timesteps. Provide a longer boundary condition"
                        )
                    node._boundary_temperature = boundary_conditions[node_k][
                        :number_of_timesteps
                    ]
                except KeyError:
                    raise BoundaryConditionNotProvided(
                        f"Node {node._idx}: the node is a {node._node_type} node, but no boundary temperature is provided.\nPlease provide a boundary condition"
                    )
            else:
                # node._boundary_temperature = np.zeros(number_of_timesteps)
                pass

    def load_hydraulic_boundary_conditions(
        self,
        nodes_boundary_conditions: dict,
        branches_boundary_conditions: dict,
        number_of_timesteps: int,
    ):
        """
        This method execute:
            self.load_nodes_mass_flow_rate_boundary_condition(nodes_boundary_conditions, number_of_timesteps)
            self.load_branches_pump_pressure_raise_boundary_condition(branches_boundary_conditions, number_of_timesteps)


        Parameters
        ----------
        nodes_boundary_conditions : dict
            Dictionary with the following sintax (mass flow rates of the supply/demend nodes):
                "1" : np.array([2,2.5,2.5,....]),
                "3" : np.array([-2,-5,-3,....]),
                "7" : np.array([2.3,2.6,2,....]),
                .
                .
                ..
        branches_boundary_conditions : dict
            Dictionary with the following sintax (pump pressure raise for the branches):
                "1" : np.array([2,2.5,2.5,....]),
                "3" : np.array([-2,-5,-3,....]),
                "7" : np.array([2.3,2.6,2,....]),
                .
                .
                ..
        number_of_timesteps : int
            Number of timesteps to be considered.

        Returns
        -------
        None.

        """
        self.load_nodes_mass_flow_rate_boundary_condition(
            nodes_boundary_conditions, number_of_timesteps
        )
        self.load_branches_pump_pressure_raise_boundary_condition(
            branches_boundary_conditions, number_of_timesteps
        )

    def load_nodes_mass_flow_rate_boundary_condition(
        self, boundary_conditions: dict, number_of_timesteps: int
    ):
        """
        This method loads the mass flow rate boundary conditions of the nodes.
        Use kg/s as reference unit.
        Dispatch nodes are those without an entering/exiting mass flow rate. Their boundary condition is zero.
        Provide negative flow rates for supply nodes and positive for demand nodes.

        Parameters
        ----------
        boundary_conditions : dict
            Dictionary with the following sintax:
                "1" : np.array([2,2.5,2.5,....]),
                "3" : np.array([-2,-5,-3,....]),
                "7" : np.array([2.3,2.6,2,....]),
                .
                .
                .

        number_of_timesteps: int
            Number of timesteps to be considered

        Returns
        -------
        None.

        """

        for node_k, node in self._nodes_object_dict.items():
            if node._node_type in ["demand", "supply"]:
                try:
                    if len(boundary_conditions[node_k]) < number_of_timesteps:
                        raise ValueError(
                            f"Node {node._idx}: the boundary condition for the node is shorter than the number of timesteps. Provide a longer boundary condition"
                        )
                    node._boundary_mass_flow_rate = boundary_conditions[node_k][
                        :number_of_timesteps
                    ]
                except KeyError:
                    raise BoundaryConditionNotProvided(
                        f"Node {node._idx}: the node is a {node._node_type} node, but no boundary mass flow rate is provided.\nPlease provide a boundary condition"
                    )
            else:
                node._boundary_mass_flow_rate = np.zeros(number_of_timesteps)

    def load_branches_pump_pressure_raise_boundary_condition(
        self, boundary_conditions: dict, number_of_timesteps: int
    ):
        """
        This method loads the pump pressure rise boundary conditions of the branches.
        Use Pa as reference unit.
        Provide positive values

        Parameters
        ----------
        boundary_conditions : dict
            Dictionary with the following sintax:
                "1" : np.array([2,2.5,2.5,....]),
                "3" : np.array([-2,-5,-3,....]),
                "7" : np.array([2.3,2.6,2,....]),
                .
                .
                .

        number_of_timesteps: int
            Number of timesteps to be considered

        Returns
        -------
        None.

        """

        for branch_k, branch in self._branches_object_dict.items():
            try:
                if len(boundary_conditions[branch_k]) < number_of_timesteps:
                    raise ValueError(
                        f"Branch {branch._idx}: the boundary condition for the branch is shorter than the number of timesteps. Provide a longer boundary condition"
                    )
                branch._pump_pressure_raise = boundary_conditions[branch_k][
                    :number_of_timesteps
                ]
            except KeyError:
                branch._pump_pressure_raise = np.zeros(number_of_timesteps)

    def solve_hydraulic_balance(self, timestep: int):
        # Boundary condition
        q = self._generate_hydraulic_balance_boundary_condition(timestep)
        if q[0 : self._nodes_number].sum() > 1e-10:
            raise ValueError(
                f"Timestep {timestep}: input - output mass flow rates not equal. Mass balance cannot be solved"
            )
        # First try vector
        x0 = self._generate_hydraulic_balance_starting_vector()
        x = root(hydraulic_balance_system, x0, args=(q, self), method="lm")
        # print(f"\n############## timestep {timestep} ###########")
        # print("m0: ", [f"{m:.2f}" for m in x0[: self._branches_number]])
        # # print("f0: ", [f"{f:.4f}" for f in x0[-self._branches_number :]])
        # # print(
        # #     "p0: ",
        # #     [f"{f:.0f}" for f in x0[self._nodes_number : -self._branches_number]],
        # # )
        # print(q[: self._nodes_number])
        # print("m: ", [f"{m:.2f}" for m in x.x[: self._branches_number]])
        # print("f: ", [f"{f*100:.4f}" for f in x.x[-self._branches_number :]])

        # print(
        #     "p: ",
        #     [f"{f:.3f}" for f in x.x[self._nodes_number : -self._branches_number]],
        # )
        self._set_hydraulic_balance_results_vector(x.x)
        return x.x, q, x0

    def solve_thermal_balance(self, timestep: int, time_interval: int = 3600):
        # Boundary condition
        q = self._generate_thermal_balance_boundary_condition(timestep, time_interval)
        # First try vector
        ################################################################
        # DO NOT USE thermal_balance_system_optimization: WRONG RESULTS
        # x0 = self._generate_thermal_balance_starting_vector()
        # x = root(
        #     thermal_balance_system_optimization,
        #     x0,
        #     args=(q, self, time_interval),
        #     method="hybr",
        # )
        ################################################################
        x, A, q = thermal_balance_system_inverse(self, q, time_interval)
        self._set_thermal_balance_results_vector(x)

        return x, A, q

    def _generate_hydraulic_balance_boundary_condition(self, timestep):
        q = []
        for node in self._nodes_object_ordered_list:
            try:
                q.append(node._boundary_mass_flow_rate[timestep])
            except AttributeError:
                raise BoundaryConditionNotProvided(
                    f"Node {node._idx}: boundary condition not provided. Be sure to load the boundaries before simulation. "
                )
            except IndexError:
                raise IndexError(
                    f"Node {node._idx}: selected timestep {timestep} longer than boundary conditions. "
                )
        for branch in self._branches_object_ordered_list:
            try:
                q.append(branch._pump_pressure_raise[timestep])
            except AttributeError:
                raise BoundaryConditionNotProvided(
                    f"Branch {branch._idx}: boundary condition not provided. Be sure to load the boundaries before simulation. "
                )
            except IndexError:
                raise IndexError(
                    f"Branch {branch._idx}: selected timestep {timestep} longer than boundary conditions. "
                )
        [q.append(0) for i in range(self._branches_number)]
        return np.array(q)

    def _generate_hydraulic_balance_starting_vector(self):
        """
        Generate the starting guess vector from the value of the previous timestep
        stored in the nodes and branches abjects

        Returns
        -------
        x : np.array
            array with the first try value [mass_flow_rates, pressures, f_factors].

        """
        # branches_mass_flow_rates = np.array(
        #     [
        #         Branch._starting_mass_flow_rate
        #         for branch in self._branches_object_ordered_list
        #     ]
        # )
        # nodes_pressures = np.array(
        #     [Node._starting_pressure for node in self._nodes_object_ordered_list]
        # )
        # branches_friction_factors = np.array(
        #     [
        #         Branch._starting_friction_factor
        #         for branch in self._branches_object_ordered_list
        #     ]
        # )
        r = np.random.randint(50) / 10
        branches_mass_flow_rates = np.array(
            [
                branch._mass_flow_rate + r / 10
                for branch in self._branches_object_ordered_list
            ]
        )
        nodes_pressures = np.array(
            [node._node_pressure + r * 10 for node in self._nodes_object_ordered_list]
        )
        branches_friction_factors = np.array(
            [
                Branch._starting_friction_factor + 0.0001 * r
                for branch in self._branches_object_ordered_list
            ]
        )
        return np.hstack(
            [branches_mass_flow_rates, nodes_pressures, branches_friction_factors]
        )

    def _generate_thermal_balance_boundary_condition(self, timestep, time_interval):
        """
        Generate thermal balance know term
        Equations are node mixing balances and branches energy balance

        Parameters
        ----------
        timestep : int
            This is the current time step the simulation is calculating.
        time_interval : int
            The time interval between timestep in seconds.

        Raises
        ------
        IndexError
            If the timestep is over the boundary condition.

        Returns
        -------
        np.array
            array of the known terms.

        """

        q = []
        for node in self._nodes_object_ordered_list:
            try:
                q.append(node._boundary_temperature[timestep])
            except AttributeError:
                q.append(0.0)  # This is the case of a
            except IndexError:
                raise IndexError(
                    f"Node {node._idx}: selected timestep {timestep} longer than boundary conditions. "
                )
        day = self.timestep_array[timestep].dayofyear

        for branch in self._branches_object_ordered_list:
            T_ground = self._soil_obj.get_soil_temperature(day, branch._pipe_depth)
            C = branch.get_dynamic_capacity()
            f_loss = branch.ground_loss_factor

            q.append(C * branch._branch_temperature / time_interval + f_loss * T_ground)
        return np.array(q)

    def _generate_thermal_balance_starting_vector(self):
        """
        Generate the starting guess vector from the value of the previous timestep
        stored in the nodes and branches abjects

        Returns
        -------
        x : np.array
            array with the first try value [nodes_temperatures, branches_temperatures].

        """
        nodes_temperatures = np.array(
            [node._node_temperature for node in self._nodes_object_ordered_list]
        )
        branches_temperatures = np.array(
            [
                branch._branch_temperature
                for branch in self._branches_object_ordered_list
            ]
        )
        return np.hstack([nodes_temperatures, branches_temperatures])

    def _set_hydraulic_balance_results_vector(self, x):
        """
        Set the new status variables to the nodes and branches objects

        Parameters
        ----------
        x : np.array
            array with the results of the hydraulic [mass_flow_rates, pressures, f_factors].

        Returns
        -------
        None.

        """
        new_branches_mass_flow_rates = x[: self._branches_number]
        new_nodes_pressures = x[
            self._branches_number : (self._branches_number + self._nodes_number)
        ]
        new_branches_friction_factors = x[
            (self._branches_number + self._nodes_number) :
        ]
        for mass, friction, branch in zip(
            new_branches_mass_flow_rates,
            new_branches_friction_factors,
            self._branches_object_ordered_list,
        ):
            branch._mass_flow_rate = mass
            branch._friction_factor = friction
        for pressure, node in zip(
            new_nodes_pressures,
            self._nodes_object_ordered_list,
        ):
            node._node_pressure = pressure

    def _set_thermal_balance_results_vector(self, x):
        """
        Set the new status variables to the nodes and branches objects.

        Parameters
        ----------
        x : np.array
            array with the results of the hydraulic [nodes_temperatures, branches_temperatures].

        Returns
        -------
        None.

        """
        new_nodes_temperatures = x[: self._nodes_number]
        new_branches_temperatures = x[
            self._nodes_number : (self._branches_number + self._nodes_number)
        ]
        for temp, node in zip(
            new_nodes_temperatures,
            self._nodes_object_ordered_list,
        ):
            node._node_temperature = temp
        for temp, branch in zip(
            new_branches_temperatures,
            self._branches_object_ordered_list,
        ):
            branch._branch_temperature = temp

    def save_results(self):
        self.save_hydraulic_results()
        self.save_thermal_results()

    def save_hydraulic_results(self):
        if self.output_path == None:
            logging.warning(
                f"Output folder not provided, results will be saved in the working directory"
            )
            self.output_path = os.path.join(".")
            self._create_output_folder()
        self.save_nodes_pressures()
        self.save_branches_mass_flow_rates()

    def save_thermal_results(self):
        if self.output_path == None:
            logging.warning(
                f"Output folder not provided, results will be saved in the working directory"
            )
            self.output_path = os.path.join(".")
            self._create_output_folder()
        self.save_nodes_temperatures()
        self.save_branches_temperatures()

    def save_nodes_pressures(self):
        nodes_pressures_header = ""
        matrix = []
        for node_k, node in self._nodes_object_dict.items():
            nodes_pressures_header += node_k + ","
            matrix.append(node._node_pressure_array)
        matrix = np.array(matrix).transpose()
        with open(os.path.join(self.output_path, "NodesPressures.csv"), "w") as nodes:
            nodes.write(nodes_pressures_header + "\n")
            np.savetxt(nodes, matrix, delimiter=",", fmt="%.0f")

    def save_branches_mass_flow_rates(self):
        branches_pressures_header = ""
        matrix = []
        for branch_k, branch in self._branches_object_dict.items():
            branches_pressures_header += branch_k + ","
            matrix.append(branch._mass_flow_rate_array)
        matrix = np.array(matrix).transpose()
        with open(
            os.path.join(self.output_path, "BranchMassFlowRates.csv"), "w"
        ) as branches:
            branches.write(branches_pressures_header + "\n")
            np.savetxt(branches, matrix, delimiter=",", fmt="%.2f")

    def save_nodes_temperatures(self):
        nodes_temperatures_header = ""
        matrix = []
        for node_k, node in self._nodes_object_dict.items():
            nodes_temperatures_header += node_k + ","
            matrix.append(node._node_temperature_array)
        matrix = np.array(matrix).transpose()
        with open(
            os.path.join(self.output_path, "NodesTemperatures.csv"), "w"
        ) as nodes:
            nodes.write(nodes_temperatures_header + "\n")
            np.savetxt(nodes, matrix, delimiter=",", fmt="%.2f")

    def save_branches_temperatures(self):
        branches_temperatures_header = ""
        matrix = []
        for branch_k, branch in self._branches_object_dict.items():
            branches_temperatures_header += branch_k + ","
            matrix.append(branch._branch_temperature_array)
        matrix = np.array(matrix).transpose()
        with open(
            os.path.join(self.output_path, "BranchTemperatures.csv"), "w"
        ) as branches:
            branches.write(branches_temperatures_header + "\n")
            np.savetxt(branches, matrix, delimiter=",", fmt="%.2f")

    @classmethod
    def from_shapefiles(
        cls,
        nodes_file: str,
        branches_file: str,
        soil_obj: Soil,
        output_path=None,
        temperature_mode="Heating",
    ):
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

        soil_obj: Soil
            Soil object to associate to the network (with soil properties)

        output_path: str
            path to where results are save. If it does not exist the tool creates it


        Returns
        -------
        None.

        """

        nodes_df = gpd.read_file(nodes_file)
        branches_df = gpd.read_file(branches_file)

        nodes_dict = {}
        branches_dict = {}

        for nodes_df_idx in nodes_df.index:
            node = nodes_df.loc[nodes_df_idx]
            nodes_dict[node["id"]] = {
                "id": str(node["id"]),
                "node type": str(node["node_type"]),
                "x": node["geometry"].x,
                "y": node["geometry"].y,
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
                "pipe ext diameter [m]": str(branch["p_de [m]"]),
                "roughness [-]": str(branch["roughness"]),
                "depth [m]": str(branch["depth [m]"]),
                "pipe thickness [m]": str(branch["pipe_t [m]"]),
                "insulation thickness [m]": str(branch["ins_t [m]"]),
                "pipe conductivity [W/(mK)]": str(branch["p_l [W/mK]"]),
                "insulation conductivity [W/(mK)]": str(branch["i_l [W/mK]"]),
            }
            # add possible pipe lenght
            if "lenght [m]" in branch.index:
                branches_dict[branch["id"]]["pipe length [m]"] = str(
                    branch["lenght [m]"]
                )
        return cls(
            nodes_dict=nodes_dict,
            branches_dict=branches_dict,
            soil_obj=soil_obj,
            output_path=output_path,
            temperature_mode=temperature_mode,
        )
