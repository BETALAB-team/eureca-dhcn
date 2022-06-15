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
from scipy.optimize import lsq_linear
import logging


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


def thermal_balance_system_inverse(network, q, time_interval, timestep):
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
            (network._branches_number + network._nodes_number),
            (network._branches_number + network._nodes_number),
        ]
    )
    # q = np.append(q, [0] * network._branches_number)

    for node in network._nodes_object_ordered_list:
        if node._node_type == "supply":
            AT[node._unique_matrix_idx, node._unique_matrix_idx] = 1
        else:
            total_entering_flow_rate = 0.0
            for (
                supply_branch_id,
                supply_branch,
            ) in node._supply_branches_objects.items():
                if supply_branch._mass_flow_rate >= 0:
                    total_entering_flow_rate += supply_branch._mass_flow_rate
                    AT[
                        node._unique_matrix_idx,
                        network._nodes_number + supply_branch._unique_matrix_idx,
                    ] = supply_branch._mass_flow_rate
            for (
                demand_branch_id,
                demand_branch,
            ) in node._demand_branches_objects.items():
                if demand_branch._mass_flow_rate < 0:
                    total_entering_flow_rate += -1 * demand_branch._mass_flow_rate
                    AT[
                        node._unique_matrix_idx,
                        network._nodes_number + demand_branch._unique_matrix_idx,
                    ] = (
                        -1 * demand_branch._mass_flow_rate
                    )
            AT[node._unique_matrix_idx, node._unique_matrix_idx] = (
                -1 * total_entering_flow_rate
            )
    mass_flow_rates = []
    mass_flow_rates_AT_0 = []
    mass_flow_rates_q_0 = []
    for branch in network._branches_object_ordered_list:
        mass_flow_rates.append(branch._mass_flow_rate)
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
        if branch.check_courant_stability(time_interval) * 0.1 < np.abs(
            branch._mass_flow_rate
        ):
            AT[line, supply_node._unique_matrix_idx] = (
                C / (2 * time_interval) - G + f_loss / 2
            )
            AT[line, line] = C / (2 * time_interval) + G + f_loss / 2
        else:
            # In case the mass flow rate is really low
            AT[line, line] = C / time_interval + f_loss
        # # Equation for average between exit and entering temperatur
        # AT[line + network._branches_number, supply_node._unique_matrix_idx] = 1
        # AT[line + network._branches_number, line + network._branches_number] = 1
        # AT[line + network._branches_number, line] = -2
        mass_flow_rates_AT_0.append(C / time_interval + f_loss)
        day = network.timestep_array[timestep].dayofyear
        t_ground = network._soil_obj.get_soil_temperature(day, branch._pipe_depth)
        mass_flow_rates_q_0.append(
            +f_loss * t_ground + C * branch._branch_temperature / time_interval
        )
        if t_ground > 20:
            logging.warning(f"T ground branch {branch._idx}: {t_ground} °C")
    mass_flow_rates = np.array(mass_flow_rates)
    if np.linalg.norm(mass_flow_rates) > 1e-4:
        try:
            x = np.linalg.solve(AT, q)
            # raise np.linalg.LinAlgError
        except np.linalg.LinAlgError:
            logging.warning(f"Thermal system not solved, trying with iterative method")
            res = lsq_linear(AT, q)
            x = res.x
            logging.warning(f"Iterative method provide a {res.cost} residual")
            if not res.success:
                logging.error(f"Iterative method didn't succed")
        # Calc T average branch
        t_branch = np.zeros(network._branches_number)
        for branch in network._branches_object_ordered_list:
            if branch._mass_flow_rate < 0:
                supply_node = branch._demand_node_object
            else:
                supply_node = branch._supply_node_object
            if branch.check_courant_stability(time_interval) * 0.1 < np.abs(
                branch._mass_flow_rate
            ):
                t_branch[branch._unique_matrix_idx] = (
                    x[supply_node._unique_matrix_idx]
                    + x[network._nodes_number + branch._unique_matrix_idx]
                ) / 2
            else:
                t_branch[branch._unique_matrix_idx] = x[
                    network._nodes_number + branch._unique_matrix_idx
                ]
    else:
        logging.warning(f"Mass flow rates norm 0")
        # Case all mass flows 0
        # Low mass flow rates case
        x = np.array([np.nan] * (network._nodes_number + network._nodes_number))
        t_branch = np.linalg.solve(
            np.diag(np.array(mass_flow_rates_AT_0)), np.array(mass_flow_rates_q_0)
        )
    solution = np.hstack(
        [x[: network._nodes_number], t_branch, x[network._nodes_number :]]
    )
    return solution, AT, q
