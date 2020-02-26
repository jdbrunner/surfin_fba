


import sys
nosolver = 0
try:
    import cplex as cp
except:
    nosolver += 1

try:
    import gurobipy as gb
except:
    nosolver += 1

if nosolver == 2:
    print("Surfin_FBA requires either CPLEX or GUROBI with associated python modules installed.")
    sys.exit()

from surfinFBA.Surfin_FBA import *


def version():
    return "Surfin_FBA version 0.8 2/26/2020"
