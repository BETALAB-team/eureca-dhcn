import numpy as np
import os
import json
from node import Node
from branch import Branch

p0 = -2

p1 = -8

p3 = 7

p5 = 3

timestep = 10

temp_ground = 13

with open(os.path.join("network_config", "nodes.json"), "r") as outfile:
    nodes = json.load(outfile)
with open(os.path.join("network_config", "branches.json"), "r") as outfile:
    branches = json.load(outfile)
n_nodes = len(nodes.keys())
n_branches = len(branches.keys())

if n_branches >= n_nodes:
    raise ValueError("The network has rings inside, delete them")
# Assignement of the index
matrix_idx_nodes = 0
nodes_list = {}
for node_id, node in nodes.items():
    node_obj = Node(
        idx=node["id"],
        node_type=node["type"],
        x=node["x"],
        y=node["y"],
        supply_nodes=node["supply nodes"],
        demand_nodes=node["demand nodes"],
        matrix_idx=matrix_idx_nodes,
    )

    nodes_list[node["id"]] = node_obj
    matrix_idx_nodes += 1
# Branches
matrix_idx_branches = 0
branches_list = {}
resistances = np.zeros(n_branches)
for branch_id, branch in branches.items():
    branch_obj = Branch(
        idx=branch["id"],
        pipe_diameter=branch["pipe diameter [m]"],  # [m]
        supply_node=nodes_list[branch["supply node"]],
        demand_node=nodes_list[branch["demand node"]],
        matrix_idx=matrix_idx_branches,
    )
    branches_list[branch["id"]] = branch_obj
    resistances[matrix_idx_branches] = branch_obj.get_hydraulic_resistance()
    nodes_list[branch["demand node"]].add_supply_branch(branch_obj)
    # TODO: Add resistance first try
    matrix_idx_branches += 1  # Assignement of the index
connection_matrix = np.zeros([n_nodes, n_branches])
for branch_id, branch in branches_list.items():
    column = branch.matrix_idx
    connection_matrix[branch.supply_node.matrix_idx, column] = -1
    connection_matrix[branch.demand_node.matrix_idx, column] = 1
#%% Matrix A
A = np.vstack(
    [
        np.hstack([connection_matrix, np.zeros([n_nodes, n_nodes])]),
        np.hstack([np.diag(resistances), connection_matrix.transpose()]),
    ]
)

q = np.array([p0, p1, 0, p3, 0, p5, 0, 0, 0, 0, 0])

# delete one of the node balances
A = A[1:, :]
q = q[1:]

# P5 = 0
A = np.vstack(
    [
        A,
        np.array([0] * (n_branches + n_nodes - 1) + [1]).reshape(
            1, (n_branches + n_nodes)
        ),
    ]
)
q = np.append(q, [0])
# P5 = 0


sol = np.linalg.solve(A, q)

# %% Set mass flow rate to branches
flow_rates = sol[:n_branches]
for branch in branches_list.values():
    branch.set_timestep_flow_rate(flow_rates[branch.matrix_idx])
# %% Problema Termico

AT = np.zeros([(n_branches + n_nodes), (n_branches + n_nodes)])
qT = np.zeros([(n_branches + n_nodes)])

for node_id, node in nodes_list.items():
    if node.node_type == "supply":
        AT[node.matrix_idx, node.matrix_idx] = 1
        qT[node.matrix_idx] = node.get_supply_temperature()
    else:
        total_entering_flow_rate = 0.0
        for supply_branch_id, supply_branch in node.supply_branches.items():
            entering_flow_rate = supply_branch.timestep_flow_rate
            total_entering_flow_rate += entering_flow_rate
            AT[node.matrix_idx, n_nodes + supply_branch.matrix_idx] = entering_flow_rate
        AT[node.matrix_idx, node.matrix_idx] = -1 * total_entering_flow_rate
        qT[node.matrix_idx] = 0
for branch_id, branch in branches_list.items():
    line = n_nodes + branch.matrix_idx
    capacity = branch.get_thermal_capacity()
    g_mul_cp = branch.get_flow_rate_cp()
    loss_factor = branch.get_thermal_conductance()
    AT[line, branch.supply_node.matrix_idx] = -1 * g_mul_cp
    AT[line, line] = capacity / timestep + g_mul_cp + loss_factor
    temp_previous_timestep = branch._previous_temp
    qT[line] = capacity * temp_previous_timestep / timestep + loss_factor * temp_ground
sol_t = np.linalg.solve(AT, qT)
