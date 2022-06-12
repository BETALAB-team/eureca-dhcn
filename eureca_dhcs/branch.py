__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math
import logging

import numpy as np
from scipy.optimize import root_scalar

from eureca_dhcs.exceptions import DuplicateBranch, WrongTemperatureMode
from eureca_dhcs.soil import Soil
from eureca_dhcs._hydraulic_system_function import (
    darcy_equation,
    darcy_equation_der,
)


class Branch:

    """
    Class for the network Branch
    """

    _idx_list = []
    _counter = 0
    # this values are use just for the first timestep calculation
    _cooling_starting_temperature = 15  # [°C]
    _heating_starting_temperature = 55  # [°C]
    _starting_mass_flow_rate = 50  # [kg/s]
    _starting_friction_factor = 0.001  # [-]
    # Darcy-Weissback equation
    # The friction factor MUST be around this value to make the hydraulic balanc stable

    def __init__(
        self,
        idx: str,
        supply_node: str,
        demand_node: str,
        pipe_ext_diameter: float,  # [m]
        pipe_thickness: float,  # [m]
        pipe_depth: float,  # [m]
        insulation_thickness: float,  # [m]
        pipe_conductivity: float,  # [W/(m/K)]
        insulation_conductivity: float,  # [m]
        pipe_len=None,  # [m]
        roughness=None,  # [-]
        starting_temperature=None,  # [°C]
        nodes_objects_dict=None,
        soil_obj=None,
        temperature_mode="Heating",
    ):
        self._idx = idx
        self._supply_node_idx = supply_node
        self._demand_node_idx = demand_node
        self._pipe_thickness = pipe_thickness
        self._pipe_ext_diameter = pipe_ext_diameter  ##############################
        self._pipe_int_diameter = self._pipe_ext_diameter - 2 * self._pipe_thickness
        self._pipe_depth = pipe_depth
        # Other properties
        self._insulation_thickness = insulation_thickness
        self._pipe_conductivity = pipe_conductivity
        self._insulation_conductivity = insulation_conductivity

        # Check if pipe len is passed
        if pipe_len == None:
            try:
                supply_node_obj = nodes_objects_dict[supply_node]
                demand_node_obj = nodes_objects_dict[demand_node]
                self._pipe_len = supply_node_obj.distance_from_node(demand_node_obj)
            except KeyError:
                raise KeyError(
                    f"Branch {idx}. Supply/demand node key not found: {supply_node}, {demand_node}"
                )
            except TypeError:
                raise NameError(
                    f"Branch {idx}. If pipe lenght not passed, the dictionary of the nodes objects must be passed to calculate the pipe lenght"
                )
        else:
            self._pipe_len = pipe_len
        if roughness == None:
            # Default roughness
            self._roughness = 0.01
        else:
            self._roughness = roughness
        # set a unique integer for the incidence matrix
        self._unique_matrix_idx = Branch._counter
        Branch._counter += 1
        # Other useful properties
        self._perimeter = self._pipe_ext_diameter * math.pi
        self._external_area = self._perimeter * self._pipe_len
        self._volume = self._pipe_int_diameter**2 / 4 * np.pi * self._pipe_len

        self._branch_temperature_array = np.array([])
        # Set some values for the dynamic simulation
        if starting_temperature != None:
            # Default starting temperature
            self._branch_temperature = starting_temperature
        else:
            if temperature_mode == "Heating":
                self._branch_temperature = self._heating_starting_temperature
            elif temperature_mode == "Cooling":
                self._branch_temperature = self._cooling_starting_temperature
            else:
                raise WrongTemperatureMode(
                    f"Branch {self._idx}: temperature mode must be or Heating or Cooling. Temperature mode: {temperature_mode}"
                )
        self._mass_flow_rate_array = np.array([])
        self._mass_flow_rate = self._starting_mass_flow_rate
        self._friction_factor = self._starting_friction_factor
        if soil_obj != None:
            self.calc_ground_loss_factor(soil_obj)

    @property
    def _idx(self) -> str:
        return self.__idx

    @_idx.setter
    def _idx(self, value: str):
        try:
            value = str(value)
        except ValueError:
            raise TypeError(f"Branch class, idx must be a str: {value}")
        if value in Branch._idx_list:
            print(Branch._idx_list)
            raise DuplicateBranch(f"Duplicate branch id: {value}")
        self.__idx = value

    @property
    def _pipe_ext_diameter(self) -> float:
        return self.__pipe_ext_diameter

    @_pipe_ext_diameter.setter
    def _pipe_ext_diameter(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe diameter must be a float: {value}"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, pipe diameter very high: {value} [m]")
        self.__pipe_ext_diameter = value

    @property
    def _pipe_int_diameter(self) -> float:
        return self.__pipe_int_diameter

    @_pipe_int_diameter.setter
    def _pipe_int_diameter(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe diameter must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, pipe diameter negative pipe diameter: {value} [m]"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, pipe diameter very high: {value} [m]")
        if value < 0.1:
            logging.warning(
                f"Branch {self._idx}, pipe internal diameter very low: {value} [m]. Hydraulic system can be unstable"
            )
        self.__pipe_int_diameter = value

    @property
    def _pipe_thickness(self) -> float:
        return self.__pipe_thickness

    @_pipe_thickness.setter
    def _pipe_thickness(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe thickness must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, pipe diameter negative pipe diameter: {value} [m]"
            )
        if value > 0.1:
            logging.warning(
                f"Branch {self._idx}, pipe thickness very high: {value} [m]"
            )
        self.__pipe_thickness = value

    @property
    def _pipe_depth(self) -> float:
        return self.__pipe_depth

    @_pipe_depth.setter
    def _pipe_depth(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, pipe depth must be a float: {value}")
        if value > 10.0:
            logging.warning(f"Branch {self._idx}, pipe depth very high: {value} [m]")
        self.__pipe_depth = value

    @property
    def _insulation_thickness(self) -> float:
        return self.__insulation_thickness

    @_insulation_thickness.setter
    def _insulation_thickness(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, insulation thickness must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, insulation thickness negative: {value} [m]"
            )
        if value > 0.1:
            logging.warning(
                f"Branch {self._idx}, insulation thickness very high: {value} [m]"
            )
        self.__insulation_thickness = value

    @property
    def _pipe_conductivity(self) -> float:
        return self.__pipe_conductivity

    @_pipe_conductivity.setter
    def _pipe_conductivity(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe conductivity must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, pipe conductivity negative: {value} [m]"
            )
        if value > 100:
            logging.warning(
                f"Branch {self._idx}, pipe conductivity very high: {value} [m]"
            )
        self.__pipe_conductivity = value

    @property
    def _insulation_conductivity(self) -> float:
        return self.__insulation_conductivity

    @_insulation_conductivity.setter
    def _insulation_conductivity(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, insulation conductivity must be a float: {value}"
            )
        if value < 0.0:
            raise ValueError(
                f"Branch {self._idx}, insulation conductivity negative: {value} [m]"
            )
        if value > 10:
            logging.warning(
                f"Branch {self._idx}, insulation conductivity very high: {value} [m]"
            )
        self.__insulation_conductivity = value

    @property
    def _pipe_len(self) -> float:
        return self.__pipe_len

    @_pipe_len.setter
    def _pipe_len(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, pipe lenght must be a float: {value}")
        if value < 0.0:
            raise ValueError(f"Branch {self._idx}, lenght negative: {value} [m]")
        self.__pipe_len = value

    @property
    def _roughness(self) -> float:
        return self.__roughness

    @_roughness.setter
    def _roughness(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, pipe roughness must be a float: {value}"
            )
        if value > 2.0:
            logging.warning(f"Branch {self._idx}, roughness very high: {value} [-]")
        self.__roughness = value

    @property
    def _branch_temperature(self) -> float:
        return self.__branch_temperature

    @_branch_temperature.setter
    def _branch_temperature(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, temperature must be a float: {value}")
        if value > 200.0:
            logging.warning(f"Branch {self._idx}, temperature very high: {value} [°C]")
        if value < 2.0:
            logging.warning(f"Branch {self._idx}, temperature very low: {value} [°C]")
        self.__branch_temperature = value
        self._branch_temperature_array = np.append(
            self._branch_temperature_array, value
        )

    @property
    def _mass_flow_rate(self) -> float:
        return self.__mass_flow_rate

    @_mass_flow_rate.setter
    def _mass_flow_rate(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, mass_flow_rate must be a float: {value}"
            )
        if np.abs(value) < 1e-3:
            value = 1e-3
            logging.warning(
                f"Branch {self._idx}, mass_flow_rate very small: {value}. Unstable system. Substituted to 0.00001 kg/s"
            )
        self.__mass_flow_rate = value
        self._mass_flow_rate_array = np.append(self._mass_flow_rate_array, value)
        self._reinolds_number = (
            4 * value / (np.pi * self.get_dynamic_viscosity() * self._pipe_int_diameter)
        )
        v = value / self.get_density()  # m3/s
        self._fluid_velocity = 4 * v / (self._pipe_int_diameter**2 * np.pi)  # m /s

    @property
    def _reinolds_number(self) -> float:
        return self.__reinolds_number

    @_reinolds_number.setter
    def _reinolds_number(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, reinolds must be a float: {value}")
        if np.abs(value) < 2300:
            logging.warning(
                f"Branch {self._idx}, reinolds is going under 2300. Reinolds: {value}"
            )
        self.__reinolds_number = value

    @property
    def _fluid_velocity(self) -> float:
        return self.__fluid_velocity

    @_fluid_velocity.setter
    def _fluid_velocity(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(f"Branch {self._idx}, velocity must be a float: {value}")
        self.__fluid_velocity = value

    @property
    def _friction_factor(self) -> float:
        return self.__friction_factor

    @_friction_factor.setter
    def _friction_factor(self, value: float):
        try:
            value = float(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, friction factor must be a float: {value}"
            )
        if value > 0.1:
            logging.warning(
                f"Branch {self._idx}, friction_factor over 0.1: {value} [-]. The hydraulic system can be unstable"
            )
            # value = 0.1
        if value < 0.0:
            logging.error(
                f"Branch {self._idx}, negative friction factor: {value} [-]. The hydraulic system can be unstable"
            )
        self.__friction_factor = value

    @property
    def _pump_pressure_raise(self) -> np.array:
        return self.__pump_pressure_raise

    @_pump_pressure_raise.setter
    def _pump_pressure_raise(self, value: np.array):
        try:
            value = np.array(value)
        except ValueError:
            raise TypeError(
                f"Branch {self._idx}, _pump_pressure_raise must be a np.array: {value}"
            )
        if np.any(value < 0):
            logging.warning(
                f"Branch {self._idx}: branch with negative pump pressure raise. Calculation will continue.."
            )
        self.__pump_pressure_raise = value

    def get_density(self):
        # https://www.researchgate.net/publication/222573141_Enhanced_modeling_of_moisture_equilibrium_and_transport_in_cementitious_materials_under_arbitrary_temperature_and_relative_humidity_history/figures?lo=1
        t = self._branch_temperature + 273.15
        d = 1.54e-8 * t**3 - 1.85e-5 * t**2 + 6.65e-3 * t + 0.247  # [g/cm3]
        return d * 1000

    def get_specific_heat(self):
        # TODO: put correlation for spec_heat - self._temperature
        return 4182  # J/kg K

    def get_dynamic_viscosity(self):
        # https://www.researchgate.net/publication/222573141_Enhanced_modeling_of_moisture_equilibrium_and_transport_in_cementitious_materials_under_arbitrary_temperature_and_relative_humidity_history/figures?lo=1
        t = self._branch_temperature + 273.15
        # mu = (
        #     3.38e-8 * t**4 - 4.63e-5 * t**3 + 2.37e-2 * t**2 + 5.45 * t + 470
        # )  # [kg/(m s)]
        # https://powderprocess.net/Tools_html/Data_Diagrams/Water_Properties_Correlations.html
        mu = np.exp(-3.7188 + 578.919 / (-137.546 + t)) / 1000  # [kg/(m s)]
        return mu

    def get_dynamic_capacity(self):
        return self.get_density() * self._volume * self.get_specific_heat()

    def calc_soil_linear_resistance(self, soil: Soil):
        # Based onn ASHRAE DHC Handbook
        try:
            soil_conductivity = soil._soil_conductivity
        except Exception:
            raise TypeError(
                f"Branch {self._idx}, calc_ground_resistance \n You must provide a proper soil object with a soil conductivity"
            )
        r_0 = self._pipe_ext_diameter / 2 + self._insulation_thickness
        d = self._pipe_depth
        if d / r_0 > 4:
            self._soil_linear_resistance = np.log(2 * d / r_0) / (
                2 * np.pi * soil_conductivity
            )
        else:
            self._soil_linear_resistance = np.log(
                d / r_0 + np.sqrt(d**2 / r_0**2 - 1)
            ) / (2 * np.pi * soil_conductivity)
            if d / r_0 < 2:
                logging.warning(
                    f"Branch {self._idx}, _soil_linear_resistance calculation. Ratio d/r0 lower than 2"
                )

    def calc_insulation_resistance(self):
        self._insulation_linear_resistance = np.log(
            (self._pipe_ext_diameter + self._insulation_thickness * 2)
            / self._pipe_ext_diameter
        ) / (2 * np.pi * self._insulation_conductivity)

    def calc_pipe_resistance(self):
        self._pipe_linear_resistance = np.log(
            self._pipe_ext_diameter / self._pipe_int_diameter
        ) / (2 * np.pi * self._pipe_conductivity)

    def calc_ground_loss_factor(self, soil):
        self.calc_soil_linear_resistance(soil)
        self.calc_insulation_resistance()
        self.calc_pipe_resistance()

        R = (
            self._pipe_linear_resistance
            + self._insulation_linear_resistance
            + self._soil_linear_resistance
        )  # [(m K)/W]
        U = self._pipe_len / R  # [W/K]
        self.ground_loss_factor = U  # W/k

    def check_courant_stability(self, time_interval):
        if (
            time_interval < self._pipe_len / self._fluid_velocity
            and self._mass_flow_rate > 5e-5
        ):
            logging.warning(
                f"Branch {self._idx}, the temperature profile can be unstable. Subdivide the branch or increase the timestep (at least to {self._pipe_len / self._fluid_velocity:.0f} s)"
            )

    def calc_hydraulic_resistance(self, mass_flow_rate):
        reinolds = np.abs(
            (
                4
                * mass_flow_rate
                / (np.pi * self.get_dynamic_viscosity() * self._pipe_int_diameter)
            )
        )

        if reinolds > 2300:
            sol = root_scalar(
                darcy_equation,
                x0=0.02,
                fprime=darcy_equation_der,
                method="newton",
                args=(reinolds, self._roughness, self._pipe_int_diameter),
            )
            if not sol.converged:
                logging.warning(f"Branch {self._idx}: friction factor not converged")
            friction_factor = sol.root
        elif reinolds < 640:
            friction_factor = 0.09
        else:
            friction_factor = 64 / reinolds
        resistance = (
            8
            * friction_factor
            * self._pipe_len
            / (self.get_density() * np.pi**2 * self._pipe_int_diameter**5)
        )
        return resistance, friction_factor

    # def get_ground_temperature(self):
    #     # TODO: put real_t_ground
    #     return 13.0

    # # First try flow rate
    # # 0.1 m/s flow rate
    # # self.first_try_flow_rate = 0.1 * 1000 * self.pipe_diameter**2 / 4
    # def get_hydraulic_resistance(self, flow_rate=None):
    #     # TODO example to try
    #     # example just to try
    #     resistance = 2
    #     return resistance

    # def get_thermal_conductance(self):
    #     # TODO example to try
    #     U = 2  # [W/m2K]
    #     return self.external_area * U  # [W/K]

    # def set_timestep_flow_rate(self, flow):
    #     self.timestep_flow_rate = flow

    # def get_thermal_capacity(self):
    #     V = math.pi * self.pipe_len * self.pipe_diameter ** 2 / 4  # [m3]
    #     return V * rho * cp  # [kJ/K]

    # def get_flow_rate_cp(self):
    #     return self.timestep_flow_rate * cp
