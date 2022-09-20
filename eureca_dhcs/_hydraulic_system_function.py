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
import logging
from eureca_dhcs.exceptions import (
    HydraulicSystemNotSolved,
)


def darcy_equation(ff, reinolds, roughness, diameter):
    return 1 / np.sqrt(ff) + 2 * np.log10(
        roughness / (3.7 * diameter) + 2.51 / (reinolds * np.sqrt(ff))
    )


def darcy_equation_der(ff, reinolds, roughness, diameter):
    # https://www.chegg.com/homework-help/questions-and-answers/2-colebrook-equation-friction-factor-turbulent-pipe-flow-given-2log-e-d-251-f-air-flow-tub-q26706228
    return (
        -0.5
        * ff ** (-3 / 2)
        * (
            1
            + 2.18261
            / reinolds
            / (roughness / (3.7 * diameter) + 2.51 / (reinolds * np.sqrt(ff)))
        )
    )


def hydraulic_balance_system_SIMPLE(x0, q, network, timestep, n_iter_max=50):
    # Boundary conditions
    # h: pressures, G_ext: massflow rates
    h = q[network._nodes_number :]
    G_ext = q[: network._nodes_number]
    # Tollerances and first try vectors
    tol_P = 100
    tol_G = 1
    G0 = x0[: network._branches_number]
    P0 = x0[network._branches_number :]
    # Little changes to avoid zero divisions in the solution
    G0 = G0 + np.linspace(-1, 1, len(G0)) * 2.0
    P0 = P0 + np.linspace(-1, 1, len(P0)) * 5000.0
    n_iter = 0
    while tol_P > 50 or tol_G > 0.01:
        n_iter += 1
        if n_iter > n_iter_max:
            logging.error(
                f"Timestep {timestep}: maximum number of iteration for the solution of the hydraulic system reached"
            )
            raise HydraulicSystemNotSolved(
                f"Timestep {timestep}: maximum number of iteration for the solution of the hydraulic system reached"
            )
        # Iterative scheme
        if np.linalg.norm(G0) > 1e-5:
            # The case of non-null mass flow rates
            R = np.array(
                [
                    branch.calc_hydraulic_resistance(G0[i])[0] * np.abs(G0[i])
                    if np.abs(G0[i]) > 1e-6
                    else 1e-6
                    for i, branch in enumerate(network._branches_object_ordered_list)
                ]
            )
            Y = np.diag(1 / np.abs(R))
            P_star = P0
            if tol_G > 0.01:
                # Case where G is changed
                G_star = Y @ network._adjacency_matrix.transpose() @ P_star + Y @ h
                # new calculation of Y
                _R = np.array(
                    [
                        branch.calc_hydraulic_resistance(G_star[i])[0]
                        * np.abs(G_star[i])
                        if np.abs(G_star[i]) > 1e-6
                        else 1e-6
                        for i, branch in enumerate(
                            network._branches_object_ordered_list
                        )
                    ]
                )
                Y_star = np.diag(1 / np.abs(_R))
            else:
                # Otherwise it keeps the old G
                G_star = G0
                Y_star = Y
            ##########################################################
            # Delta calculation creating the equation 11-12 combined #
            ##########################################################
            aux_A_Q = network._adjacency_matrix
            aux_q_Q = -network._adjacency_matrix @ G_star - G_ext

            aux_A_P = -(Y_star @ network._adjacency_matrix.transpose())
            aux_q_P = np.zeros(network._branches_number)

            aux_A = np.vstack(
                [
                    np.hstack(
                        [
                            aux_A_Q,
                            np.zeros([network._nodes_number, network._nodes_number]),
                        ]
                    ),
                    np.hstack([np.eye(network._branches_number), aux_A_P]),
                ]
            )
            aux_q = np.hstack([aux_q_Q, aux_q_P])
            # One of the nodes mass balances is deleted and a node's pressure is set
            aux_A_ = np.delete(aux_A, network._first_supply_node_idx, 0)
            aux_q_ = np.delete(aux_q, network._first_supply_node_idx)
            aux_A_ = np.vstack(
                [aux_A_, np.zeros(network._branches_number + network._nodes_number)]
            )
            aux_q_ = np.hstack([aux_q_, 0])
            aux_A_[-1, network._branches_number + network._first_supply_node_idx] = 1
            # Solution of the system for the deltas
            delta = np.linalg.solve(aux_A_, aux_q_)

            delta_G_corr = delta[: network._branches_number]
            delta_P_corr = delta[network._branches_number :]
            # Correction of the G and P vectors
            G = G_star + delta_G_corr
            P = P_star + delta_P_corr
            # Only the G vector is kept
            # On the contrary the P vector is calculated again based on the new G vector and the system
            # The system is create using Gand the new resistance
            __R = np.array(
                [
                    branch.calc_hydraulic_resistance(G[i])[0] * np.abs(G[i])
                    if np.abs(G[i]) > 1e-6
                    else 1e-6
                    for i, branch in enumerate(network._branches_object_ordered_list)
                ]
            )
            R = np.diag(np.abs(__R))

            aux_A = np.vstack(
                [
                    np.hstack(
                        [
                            network._adjacency_matrix,
                            np.zeros([network._nodes_number, network._nodes_number]),
                        ]
                    ),
                    np.hstack([-R, network._adjacency_matrix.transpose()]),
                ]
            )
            aux_q = np.hstack([-G_ext, -h])
            # One of the nodes mass balances is deleted and a node's pressure is set
            aux_A_ = np.delete(aux_A, network._first_supply_node_idx, 0)
            aux_q_ = np.delete(aux_q, network._first_supply_node_idx)
            aux_A_ = np.vstack(
                [aux_A_, np.zeros(network._branches_number + network._nodes_number)]
            )
            aux_q_ = np.hstack([aux_q_, 0])
            aux_A_[-1, network._branches_number + network._first_supply_node_idx] = 1

            GP = np.linalg.solve(aux_A_, aux_q_)

            G = GP[: network._branches_number]
            P = GP[network._branches_number :]

            # # aux_A = np.insert(
            # #     self._adjacency_matrix.transpose(),
            # #     self._first_supply_node_idx,
            # #     0,
            # #     axis=0,
            # # )
            # # aux_A[self._first_supply_node_idx, self._first_supply_node_idx] = 1
            # # aux_q = np.insert(
            # #     np.dot(np.diag(__R), G) - h, self._first_supply_node_idx, 0
            # # )
            # # P = np.linalg.solve(aux_A, aux_q)

            tol_P = np.linalg.norm(P0 - P)
            tol_G = np.linalg.norm(G0 - G)

            P0 = P
            G0 = G
        else:
            # This is the case where the mass flow rates are zeros
            P0 = P0 * 0  # And the pumps?
            G0 = G0 * 0
            tol_P = 0.0001
            tol_G = 0.0001
    # In case there are some over pressures the system is solved again
    __R = np.array(
        [
            branch.calc_hydraulic_resistance(G0[i])[0] * np.abs(G0[i])
            if np.abs(G0[i]) > 1e-6
            else 1e-6
            for i, branch in enumerate(network._branches_object_ordered_list)
        ]
    )
    R = np.diag(np.abs(__R))

    aux_A = np.vstack(
        [
            np.hstack(
                [
                    network._adjacency_matrix,
                    np.zeros([network._nodes_number, network._nodes_number]),
                ]
            ),
            np.hstack([-R, network._adjacency_matrix.transpose()]),
        ]
    )
    aux_q = np.hstack([-G_ext, -h])

    aux_A_ = np.delete(aux_A, network._first_supply_node_idx, 0)
    aux_q_ = np.delete(aux_q, network._first_supply_node_idx)
    aux_A_ = np.vstack(
        [aux_A_, np.zeros(network._branches_number + network._nodes_number)]
    )
    aux_q_ = np.hstack([aux_q_, 0])
    aux_A_[-1, network._branches_number + network._first_supply_node_idx] = 1

    try:
        GP = np.linalg.solve(aux_A_, aux_q_)
    except np.linalg.LinAlgError:
        GP = np.hstack([G0, P0])
        logging.error(f"Timestep {timestep}: singular matrix.")
    G = GP[: network._branches_number]
    P = GP[network._branches_number :]
    # aux_A = np.insert(
    #     self._adjacency_matrix.transpose(),
    #     self._first_supply_node_idx,
    #     0,
    #     axis=0,
    # )
    # aux_A[self._first_supply_node_idx, self._first_supply_node_idx] = 1
    # aux_q = np.insert(
    #     np.dot(np.diag(__R), G) - h, self._first_supply_node_idx, 0
    # )
    # P = np.linalg.solve(aux_A, aux_q)
    ff = np.array(
        [
            branch.calc_hydraulic_resistance(G[i])[1]
            for i, branch in enumerate(network._branches_object_ordered_list)
        ]
    )

    x = np.hstack([G, P, ff])
    return x


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
        if not node._first_supply_node:
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
    system.append(
        x[network._branches_number + network._nodes_number - 1]
        - q[network._branches_number + network._nodes_number - 1]
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
            line_list = np.zeros(network._branches_number * 2 + network._nodes_number)
            # Derivative of the node mass balance for each mass flow rate
            supply_idx = node.get_supply_branches_unique_idx()
            line_list[supply_idx] = -1
            demand_idx = node.get_demand_branches_unique_idx()
            line_list[demand_idx] = 1
            system.append(line_list)
    for branch in network._branches_object_ordered_list:
        line_list = np.zeros(network._branches_number * 2 + network._nodes_number)
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
            * np.sqrt(x[branch._unique_matrix_idx] ** 2)
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
            network._branches_number * 2 + network._nodes_number
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
    line_list = np.zeros(network._branches_number * 2 + network._nodes_number)
    line_list[network._branches_number + network._nodes_number - 1] = 1
    system.append(line_list)
    return system
