"""
This script include the function that builds the hydraulic system
to pass to  scipy.optimize.fsolve
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import numpy as np


def hydraulic_balance_system(x, q, network):
    """
    This function builds the hydraulic system for the solution with f solve:
    WARNING:
        the x0 values for the friction factor MUST be around 0.02, otherwise the system is not stable

    The equation are:
        n_nodes equations: the mass flow rate balances of each node
        n_branch equations: the Darcy-Weissback pressure loss equation of each branch
        n_branch equations: the Colebrook-White friction factor equation of each branch
        1 pressure reference for a random node

    The unknown variable are (x vector):
        n_branches mass flow rates [kg/s]
        n_nodes pressures [Pa]
        n_branches friction factors [-]

    The known values are (q vector):
        nodes mass flow rates entering (-) or exiting (+) [kg/s]
        branches pumps (+) [Pa]
        known term of Colebrook-White = 0 [-]


    Parameters
    ----------
    x : np.array
        array with the first try value [mass_flow_rates, pressures, f_factors].
    q : np.array
        boundary conditions [nodes_mass_flow_rates, branches_pumps_pressure_raise, 0].
    network : Network
        network object. This is used to manage the index of the matrix with respect to branches and nodes objects.

    Returns
    -------
    system : List
        Use this function in an fsolve.

    """

    # System id the list where equations are inserted
    system = []
    # node balances
    # This equations are the mass balance for each node
    for node in network._nodes_object_ordered_list[1:]:
        supply_idx = node.get_supply_branches_unique_idx()
        demand_idx = node.get_demand_branches_unique_idx()
        system.append(
            x[demand_idx].sum() - x[supply_idx].sum() + q[node._unique_matrix_idx]
        )
    for branch in network._branches_object_ordered_list:
        # system.append(
        #     x[branch._demand_node_object._unique_matrix_idx + network._branches_number]
        #     - x[
        #         branch._supply_node_object._unique_matrix_idx + network._branches_number
        #     ]
        #     + x[branch._unique_matrix_idx]
        #     * 2
        #     * x[
        #         network._branches_number
        #         + network._nodes_number
        #         + branch._unique_matrix_idx
        #     ]
        # )

        # Darcy–Weisbach equation for each branch
        system.append(
            x[branch._demand_node_object._unique_matrix_idx + network._branches_number]
            - x[
                branch._supply_node_object._unique_matrix_idx + network._branches_number
            ]
            - q[branch._unique_matrix_idx + network._nodes_number]
            + x[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ]
            * np.sign(x[branch._unique_matrix_idx])
            * x[branch._unique_matrix_idx] ** 2
            * 8
            * branch._pipe_len
            / (np.pi**2 * branch.get_density() * branch._pipe_int_diameter**5)  #  )
        )

        # Colebrook - White equation for each branch
        # f factor is as follow
        # x[
        #         network._branches_number
        #         + network._nodes_numbes
        #         + branch._unique_matrix_idx
        #     ]
        reinolds = (
            4
            * x[branch._unique_matrix_idx]
            / (np.pi * branch.get_dynamic_viscosity() * branch._pipe_int_diameter)
        )
        if np.abs(reinolds) > 2300:
            system.append(
                1
                / np.sqrt(
                    x[
                        network._branches_number
                        + network._nodes_number
                        + branch._unique_matrix_idx
                    ]
                )
                + 2
                * np.log(
                    branch._roughness / (3.7 * branch._pipe_int_diameter)
                    + np.pi
                    * 2.51
                    * branch._pipe_int_diameter
                    * branch.get_dynamic_viscosity()
                    / (
                        4
                        * np.abs(x[branch._unique_matrix_idx])
                        * np.sqrt(
                            x[
                                network._branches_number
                                + network._nodes_number
                                + branch._unique_matrix_idx
                            ]
                        )
                    )
                )
            )
        elif np.abs(reinolds) < 640:
            system.append(
                x[
                    network._branches_number
                    + network._nodes_number
                    + branch._unique_matrix_idx
                ]
                - 0.09
            )
        else:
            system.append(
                x[
                    network._branches_number
                    + network._nodes_number
                    + branch._unique_matrix_idx
                ]
                - 64
                / (
                    4
                    * x[branch._unique_matrix_idx]
                    / (
                        np.pi
                        * branch.get_dynamic_viscosity()
                        * branch._pipe_int_diameter
                    )
                )
            )
    system.append(
        x[network._branches_number + network._nodes_number - 1]
        - q[network._branches_number + network._nodes_number - 1]
    )
    return system
