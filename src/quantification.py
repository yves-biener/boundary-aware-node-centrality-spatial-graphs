import gurobipy as gp
from gurobipy import GRB
import numpy as np
import matplotlib.pyplot as plt

# NOTE:
# Thank you Anna for providing her sources: https://github.com/uderhardtlab/truncated_graphs/blob/main/src/fit.py


def fit_piece_wise_linear(d, C, M=1000):
    """
    Fits a piecewise linear function to the given data using optimization.
    
    Parameters:
    d (array-like): Distance values.
    C (array-like): Corresponding centrality values.
    M (int, optional): Large constant for big-M constraints. Default is 1000.
    
    Returns:
    tuple: Optimized slope (m), intercept (c0), and breakpoint (b).
    """
    n = len(d) 

    model = gp.Model()
    m = model.addVar(vtype=GRB.CONTINUOUS, name="m") # Slope before breakpoint

    c0 = model.addVar(vtype=GRB.CONTINUOUS, name="c0") # y-intercept
    b = model.addVar(vtype=GRB.CONTINUOUS, lb=min(d), ub=max(d), name="b") # Breakpoint

    z = model.addVars(n, vtype=GRB.BINARY, name="z")
    epsilon = model.addVars(n, vtype=GRB.CONTINUOUS, name="epsilon")

    model.setObjective(gp.quicksum(epsilon[i] * epsilon[i] for i in range(n)), GRB.MINIMIZE)

    # Setting solver parameters for precision
    model.setParam('OptimalityTol', 1e-4) 
    model.setParam('MIPGap', 0.01)  

    for i in range(n):
        # Constraints enforcing piecewise linear fit
        model.addConstr(epsilon[i] >= (m * d[i] + c0 - C[i]) - (1 - z[i]) * M)
        model.addConstr(epsilon[i] >= -(m * d[i] + c0 - C[i]) - (1 - z[i]) * M)

        model.addConstr(epsilon[i] >= (m * b + c0 - C[i]) - z[i] * M)
        model.addConstr(epsilon[i] >= -(m * b + c0 - C[i]) - z[i] * M)

        # Binary switch for piecewise behavior
        model.addConstr(z[i] * M >= b - d[i])
        model.addConstr((1 - z[i]) * M >= d[i] - b)

    model.optimize()
    return m.X, c0.X, b.X


def plot_piece_wise_linear(d, C, m_opt, c0_opt, b_opt, measure, n, t, path):
    """
    Plots the original data and optimized piecewise linear fit.
    """
    d_curve = np.linspace(min(d), max(d), 500)
    C_curve = np.piecewise(
        d_curve,
        [d_curve <= b_opt, d_curve > b_opt],
        [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
    )

    plt.scatter(d, C, color="blue", label="Original", alpha=0.5)
    plt.plot(d_curve, C_curve, color="red", label="Optimized", linewidth=2)
    plt.xlabel("Distance to border")
    plt.ylabel(f"{measure.capitalize()} centrality")
    plt.title(f"Optimized piece-wise linear fit for {t} graph and {measure} centrality")
    plt.legend()
    plt.savefig(path, format='svg')

    fig, ax = plt.subplots()
    ax.set_title(metric_name)
    ax.set_xlabel('Distance to Bounding-Box')
    ax.set_ylabel('Centrality')
    ax.scatter(d, C, color='b', s=0.2)
    fig.savefig(path, format='svg')
