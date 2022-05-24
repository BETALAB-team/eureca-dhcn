"""
List of custom exceptions
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"


class DuplicateNode(Exception):
    pass


class WrongNodeType(Exception):
    pass


class DuplicateBranch(Exception):
    pass


class WrongTemperatureMode(Exception):
    pass


class EmptyNetworkNodes(Exception):
    pass


class BoundaryConditionNotProvided(Exception):
    pass
