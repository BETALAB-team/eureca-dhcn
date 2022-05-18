__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math
import logging

import numpy as np

from eureca_dhcs.exceptions import DuplicateBranch, WrongBranchTemperatureMode


class Soil:
    def __init__(self):
        # Based on ASHRAE handbook district heating and cooling

        self._soil_conductivity = 1.40  # Davis approximate data
        self._soil_specific_heat = 0.73
        self._water_specific_heat = 4.18
        self._soil_density = 1450
        self._moisture_content = 10  # %
        self._number_of_days = 365

        self._soil_diffusivity = (
            24
            * 3600
            * self._soil_conductivity
            / (
                1000
                * self._soil_density
                * (
                    self._soil_specific_heat
                    + self._water_specific_heat * self._moisture_content / 100
                )
            )
        )

        # self._mean_annual_surface_temperature

    def get_soil_temperature(self):
        # TODO be Modify
        return 15.0
