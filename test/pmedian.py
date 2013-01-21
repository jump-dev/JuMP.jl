# test_pmedian.py
# Solves the P-Median problem
from pulp import *
import random
import time

numFacility = 100
numCustomer = 100
numLocation = 5000

customerLocations = [random.randint(0,numLocation-1) for a in range(numCustomer)]

tic = time.time()
m = LpProblem("pmedian",LpMaximize)

# Facility locations
s = LpVariable.dicts("s", [i for i in range(numLocation)], 0, 1, LpInteger)

# Aux. variable: x_a,i = 1 iff the closest facility to a is at i
x = LpVariable.dicts("x", [(a,i) for a in range(numCustomer) for i in range(numLocation)], 0, 1, LpInteger)

# Objective: min distance
m += lpSum( x[(a,i)] * abs(customerLocations[a]-i) for a in range(numCustomer) for i in range(numLocation) )


for a in range(numCustomer):
  # Subject to linking x with s
  for i in range(numLocation):
    m += x[(a,i)] <= s[i]
  # Subject to one of x must be 1
  m += lpSum(x[(a,i)] for i in range(numLocation)) == 1

# Subject to must allocate all facilities
m += lpSum(s[i] for i in range(numLocation)) == 1

m.writeLP("pmedian_py.lp")
toc = time.time()

print tic,toc,toc-tic
