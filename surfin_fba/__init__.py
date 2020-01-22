



nosolver = 0
try:
    import cplex as cp
except:
    print("CPLEX not available.")
    nosolver += 1

try:
    import gurobipy as gb
except:
    print("Gurobi not available")
    nosolver += 1

if nosolver == 2:
    print("Surfin_FBA requires either CPLEX or GUROBI python modules installed.")
    sys.exit()

from surfin_fba.Surfin_FBA import *


def version():
    return "Surfin_FBA version 0.5 1/20/2020"
