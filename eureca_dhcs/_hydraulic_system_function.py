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
from eureca_dhcs.exceptions import TooManyBoundaryCondition, TooFewBoundaryCondition


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
    # This equation are the mass balance for each node
    for node in network._nodes_object_ordered_list:
        supply_idx = node.get_supply_branches_unique_idx()
        demand_idx = node.get_demand_branches_unique_idx()
        # This if is to skip one of the nodes
        if not node._first_supply_node:
            # system.append(
            #     x[demand_idx].sum() - x[supply_idx].sum() + q[node._unique_matrix_idx]
            # )
            if node._boundary_mass_flow_rate_undefined:
                # The mass flow rate entering in this node is unknown
                system.append(
                    x[demand_idx].sum()
                    - x[supply_idx].sum()
                    + x[
                        network._branches_number * 2
                        + network._nodes_number
                        + node._unique_matrix_idx_undefined_flow_rate
                    ]
                    + q[node._unique_matrix_idx]
                )
            else:
                system.append(
                    x[demand_idx].sum()
                    - x[supply_idx].sum()
                    + q[node._unique_matrix_idx]
                )
        # system.append(
        #     x[demand_idx].sum() - x[supply_idx].sum() + q[node._unique_matrix_idx]
        # )
        # Pressure
        # MUST BE FIXED.
        try:
            node._boundary_pressure
            system.append(
                x[network._branches_number + node._unique_matrix_idx]
                - q[
                    network._branches_number
                    + network._nodes_number
                    + node._unique_matrix_idx
                ]
            )
        except AttributeError:
            pass
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

        # Darcy–Weissbach equation for each branch
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
            * np.abs(x[branch._unique_matrix_idx])
            * x[branch._unique_matrix_idx]
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
        reinolds = np.abs(
            (
                4
                * x[branch._unique_matrix_idx]
                / (np.pi * branch.get_dynamic_viscosity() * branch._pipe_int_diameter)
            )
        )
        if reinolds > 2300:
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
        elif reinolds < 640:
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
                - 64 / reinolds
            )
    if len(system) > (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
    ):
        raise TooManyBoundaryCondition()
    elif len(system) < (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
        - 1
    ):
        raise TooFewBoundaryCondition()
    elif len(system) == (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
        - 1
    ):
        # Mass balance
        system.append(
            x[(network._branches_number * 2 + network._nodes_number) :].sum()
            + q[: network._nodes_number].sum()
        )
    return system


def hydraulic_balance_system_jac(x, q, network):
    """
    This function builds the hydraulic system jacobian for the solution with f solve:
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
    # This equation are the mass balance for each node
    for node in network._nodes_object_ordered_list:
        if not node._first_supply_node:
            line_list = np.zeros(
                network._branches_number * 2
                + network._nodes_number
                + network._nodes_number_undefined_flow_rate
            )
            # Derivative of the node mass balance for each mass flow rate
            supply_idx = node.get_supply_branches_unique_idx()
            line_list[supply_idx] = -1
            demand_idx = node.get_demand_branches_unique_idx()
            line_list[demand_idx] = 1
            if node._boundary_mass_flow_rate_undefined:
                line_list[
                    network._branches_number * 2
                    + network._nodes_number
                    + node._unique_matrix_idx_undefined_flow_rate
                ] = 1
            system.append(line_list)
            try:
                node._boundary_pressure
                line_list = np.zeros(
                    network._branches_number * 2
                    + network._nodes_number
                    + network._nodes_number_undefined_flow_rate
                )
                line_list[network._branches_number + node._unique_matrix_idx] = 1
                system.append(line_list)
            except AttributeError:
                pass
    for branch in network._branches_object_ordered_list:
        line_list = np.zeros(
            network._branches_number * 2
            + network._nodes_number
            + network._nodes_number_undefined_flow_rate
        )
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

        # Darcy–Weissbach equation for each branch
        line_list[
            branch._demand_node_object._unique_matrix_idx + network._branches_number
        ] = 1
        line_list[
            branch._supply_node_object._unique_matrix_idx + network._branches_number
        ] = -1
        # Derivative with respect to mass flow rate
        line_list[branch._unique_matrix_idx] = (
            2
            * x[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ]
            * np.sign(x[branch._unique_matrix_idx])
            * x[branch._unique_matrix_idx]
            * 8
            * branch._pipe_len
            / (np.pi**2 * branch.get_density() * branch._pipe_int_diameter**5)  #  )
        )
        # Derivative with respect to friction factor
        line_list[
            network._branches_number + network._nodes_number + branch._unique_matrix_idx
        ] = (
            np.sign(x[branch._unique_matrix_idx])
            * x[branch._unique_matrix_idx] ** 2
            * 8
            * branch._pipe_len
            / (np.pi**2 * branch.get_density() * branch._pipe_int_diameter**5)  #  )
        )

        system.append(line_list)

        # Colebrook - White equation for each branch
        # f factor is as follow
        # x[
        #         network._branches_number
        #         + network._nodes_numbes
        #         + branch._unique_matrix_idx
        #     ]
        line_list_colebrook = np.zeros(
            network._branches_number * 2
            + network._nodes_number
            + network._nodes_number_undefined_flow_rate
        )
        reinolds = np.abs(
            (
                4
                * x[branch._unique_matrix_idx]
                / (np.pi * branch.get_dynamic_viscosity() * branch._pipe_int_diameter)
            )
        )
        if reinolds > 2300:
            # Derivative with respect to mass flow rate
            line_list_colebrook[branch._unique_matrix_idx] = (
                -1
                * np.sign(x[branch._unique_matrix_idx])
                * x[branch._unique_matrix_idx] ** (-2)
                * (
                    2.51
                    * np.pi
                    * branch._pipe_int_diameter
                    * branch.get_dynamic_viscosity()
                )
                / (
                    2
                    * np.sqrt(
                        x[
                            network._branches_number
                            + network._nodes_number
                            + branch._unique_matrix_idx
                        ]
                    )
                    * (
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
            )
            # Derivative with respect to friction factor
            line_list_colebrook[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ] = -0.5 * x[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ] ** (
                -3 / 2
            )
            +(
                -2.51
                * np.pi
                * branch._pipe_int_diameter
                * branch.get_dynamic_viscosity()
                * x[
                    network._branches_number
                    + network._nodes_number
                    + branch._unique_matrix_idx
                ]
                ** (-3 / 2)
            ) / (
                4
                * np.abs(x[branch._unique_matrix_idx])
                * (
                    branch._roughness / (3.7 * branch._pipe_int_diameter)
                    + (
                        2.51
                        * np.pi
                        * branch._pipe_int_diameter
                        * branch.get_dynamic_viscosity()
                    )
                    / (
                        np.abs(x[branch._unique_matrix_idx])
                        * x[
                            network._branches_number
                            + network._nodes_number
                            + branch._unique_matrix_idx
                        ]
                        ** (1 / 2)
                    )
                )
            )
            system.append(line_list_colebrook)
        elif reinolds < 640:

            line_list_colebrook[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ] = 1
            system.append(line_list_colebrook)
        else:
            # mass flow derivative
            line_list_colebrook[branch._unique_matrix_idx] = (
                np.sign(x[branch._unique_matrix_idx])
                * 16
                * np.pi
                * branch._pipe_int_diameter
                * branch.get_dynamic_viscosity()
            ) / x[branch._unique_matrix_idx] ** 2

            # friction derivative
            line_list_colebrook[
                network._branches_number
                + network._nodes_number
                + branch._unique_matrix_idx
            ] = 1

            system.append(line_list_colebrook)
    # line_list = np.zeros(network._branches_number * 2 + network._nodes_number)
    # line_list[network._branches_number + network._nodes_number - 1] = 1
    # system.append(line_list)
    if len(system) > (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
    ):
        raise TooManyBoundaryCondition()
    elif len(system) < (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
        - 1
    ):
        raise TooFewBoundaryCondition()
    elif len(system) == (
        network._branches_number * 2
        + network._nodes_number
        + network._nodes_number_undefined_flow_rate
        - 1
    ):
        # Mass balance
        line_list = np.zeros(
            network._branches_number * 2
            + network._nodes_number
            + network._nodes_number_undefined_flow_rate
        )
        line_list[(network._branches_number * 2 + network._nodes_number) :] = 1
        system.append(line_list)
    return system
