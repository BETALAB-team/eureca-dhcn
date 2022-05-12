"""
Tests
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import os

import pytest
import numpy as np

from eureca_building.material import Material, AirGapMaterial


class TestMaterials:
    """
    This is a test class for the pytest module.
    It tests Material class and its property
    """

    def test_material(self):
        # Standard material creation
        Material("Test material")

    def test_material_comp(self):
        # Standard material creation
        Material("Test material", thick=0.100, cond=1.00, spec_heat=1000.0, dens=1000.0)

    def test_material_with_prop_wrong(self):
        # Standard material creation
        with pytest.raises(MaterialPropertyOutsideBoundaries):
            Material("Test material", cond=1000.0)
