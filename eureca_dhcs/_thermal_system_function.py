"""
This script include the function that builds the thermal system
to pass to  scipy.optimize.fsolve
or np.linalg.solve
"""

__author__ = "Enrico Prataviera"
__credits__ = ["Enrico Prataviera"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Enrico Prataviera"

import numpy as np


def thermal_balance_system_optimization(x, q, network, time_interval):
    """
    DO NOT USE! WRONG RESULTS
    This function builds the thermal system for the solution with f solve:

    The equation are:
        n_nodes equations: the mixing balances for each node
        n_branch equations: the energy balance equation of each branch

    The unknown variable are (x vector):

        n_branches temperatures [°C]n_nodes temperatures [°C]

    The known values are (q vector):
        supply nodes entering temperature [°C]
        branches capacity previous timestep term and ground loss term


    Parameters
    ----------
    x : np.array
        array with the first try value [nodes_temperatures, branches_temperatures].
    q : np.array
        boundary conditions [nodes_temperatures, branches capacity and ground loss term].
    network : Network
        network object. This is used to manage the index of the matrix with respect to branches and nodes objects.
    time_interval: int
        number of seconds per each timestep

    Returns
    -------
    system : List
        Use this function in an fsolve/root.

    """

    # System id the list where equations are inserted
    system = []
    # node balances
    # This equations are the mass balance for each node
    i = 0
    for node in network._nodes_object_ordered_list:
        if node._node_type == "supply":
            # known condition
            system.append(x[node._unique_matrix_idx] - q[node._unique_matrix_idx])
        else:
            (
                branches_idx,
                brnaches_mass_flow,
            ) = node.get_supply_branches_mass_flow_rates()
            total = brnaches_mass_flow.sum()
            vector = np.zeros(network._nodes_number + network._branches_number)
            vector[branches_idx + network._nodes_number] = brnaches_mass_flow
            vector[node._unique_matrix_idx] = -1 * total
            system.append(np.dot(vector, x))
    for branch in network._branches_object_ordered_list:
        G = branch._mass_flow_rate * branch.get_specific_heat()
        C = branch.get_dynamic_capacity()
        f_loss = branch.ground_loss_factor
        # to manage negative mass flow rates
        G = np.abs(G)
        if branch._mass_flow_rate < 0:
            supply_node = branch._demand_node_object
            demand_node = branch._supply_node_object
        else:
            supply_node = branch._supply_node_object
            demand_node = branch._demand_node_object
        system.append(
            -1 * G * x[supply_node._unique_matrix_idx]
            + x[branch._unique_matrix_idx + network._nodes_number]
            * (C / time_interval + f_loss + 1 * G)
            + 0 * G * x[demand_node._unique_matrix_idx]
            + q[branch._unique_matrix_idx + network._nodes_number]
        )
    return system


def thermal_balance_system_inverse(network, q, time_interval):
    """
    This function solve the linear system Ax=q for the thermal balance:

    The equation are:
        n_nodes equations: the mixing balances for each node
        n_branch equations: the energy balance equation of each branch

    The unknown variable are (x vector):
        n_branches temperatures [°C]n_nodes temperatures [°C]

    The known values are (q vector):
        supply nodes entering temperature [°C]
        branches capacity previous timestep term and ground loss term


    Parameters
    ----------
    q : np.array
        boundary conditions [nodes_temperatures, branches capacity and ground loss term].
    network : Network
        network object. This is used to manage the index of the matrix with respect to branches and nodes objects.
    time_interval: int
        number of seconds per each timestep

    Returns
    -------
    x : np.array
        array with the solution [nodes_temperatures, branches_temperatures].

    """

    # Problema Termico

    AT = np.zeros(
        [
            (network._branches_number * 2 + network._nodes_number),
            (network._branches_number * 2 + network._nodes_number),
        ]
    )
    # Add the branch exit temperature variable and known term
    q = np.append(q, [0] * network._branches_number)

    for node in network._nodes_object_ordered_list:
        if node._node_type == "supply":
            AT[node._unique_matrix_idx, node._unique_matrix_idx] = 1
        else:
            total_entering_flow_rate = 0.0
            for (
                supply_branch_id,
                supply_branch,
            ) in node._supply_branches_objects.items():
                total_entering_flow_rate += supply_branch._mass_flow_rate
                AT[
                    node._unique_matrix_idx,
                    network._nodes_number
                    + network._branches_number
                    + supply_branch._unique_matrix_idx,
                ] = supply_branch._mass_flow_rate
            AT[node._unique_matrix_idx, node._unique_matrix_idx] = (
                -1 * total_entering_flow_rate
            )
    for branch in network._branches_object_ordered_list:
        line = network._nodes_number + branch._unique_matrix_idx
        G = branch._mass_flow_rate * branch.get_specific_heat()
        C = branch.get_dynamic_capacity()
        f_loss = branch.ground_loss_factor
        # to manage negative mass flow rates
        G = np.abs(G)
        if branch._mass_flow_rate < 0:
            supply_node = branch._demand_node_object
            demand_node = branch._supply_node_object
        else:
            supply_node = branch._supply_node_object
            demand_node = branch._demand_node_object
        AT[line, supply_node._unique_matrix_idx] = -1 * G
        AT[
            line,
            network._nodes_number
            + network._branches_number
            + branch._unique_matrix_idx,
        ] = (
            1 * G
        )
        AT[line, line] = C / time_interval + f_loss
        # Equation for average between exit and entering temperatur
        AT[line + network._branches_number, supply_node._unique_matrix_idx] = 1
        AT[line + network._branches_number, line + network._branches_number] = 1
        AT[line + network._branches_number, line] = -2
    return np.linalg.solve(AT, q), AT, q
