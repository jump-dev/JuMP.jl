using JuMP
using Gurobi

m = Model(solver=GurobiSolver())

@defVar(m, x[1:3,1:3,1:3], Bin)

# Constraint 1 - Each row...
@addConstraint(m, row[i=1:3,val=1:3], sum(x[i,:,val]) == 1)
# Constraint 2 - Each column...
@addConstraint(m, col[j=1:3,val=1:3], sum(x[:,j,val]) == 1)
# Constraint 4 - Cells...
@addConstraint(m, cells[i=1:3,j=1:3], sum(x[i,j,:]) == 1)

# 1 2 -
# - 3 -
# - - -
@addConstraint(m, x[1,1,1] == 1)
@addConstraint(m, x[1,2,2] == 1)
@addConstraint(m, x[2,2,3] == 1)

println(m)

# Solves in presolve
presolve(m)