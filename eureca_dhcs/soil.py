__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import math
import logging

import numpy as np

from eureca_dhcs.exceptions import DuplicateBranch, WrongTemperatureMode


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
            86.4
            * self._soil_conductivity
            / (
                self._soil_density
                * (
                    self._soil_specific_heat
                    + self._water_specific_heat * self._moisture_content / 100
                )
            )
        )
        self.day_phase = -80  # Theta lag ASHRAE handbook_soil_conductivity
        self.amplitude = 13.79
        self.average_outdoor_t = 15.9

    def get_soil_temperature(self, day, depth):
        factor = -depth * np.sqrt(np.pi / (self._soil_diffusivity * 365))
        t_ground = self.average_outdoor_t - self.amplitude * np.exp(factor) * np.sin(
            2 * np.pi * (day - self.day_phase) / 365 + factor
        )

        return t_ground
