# test_dense.py
# Make a big dense matrix using PuLP
from pulp import *
import random
import time

N = 10000
M = 1000
A = []
for i in range(M):
  A.append([])
  for j in range(N):
    A[i].append(random.randint(1,10))
b = []
for i in range(M):
  b.append(random.randint(N,N+N*10))
c = []
for j in range(N):
  c.append(random.randint(1,10))


tic = time.time()
m = LpProblem("dense",LpMaximize)
v = LpVariable.dict("x",[j for j in range(N)],0,1)

m += lpSum(v[j] * c[j] for j in range(N))

for i in range(M):
  m += lpSum(v[j] * A[i][j] for j in range(N)) <= b[i]

m.writeMPS("dense_py.mps")
toc = time.time()

print tic,toc,toc-tic
