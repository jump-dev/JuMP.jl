using JuMP
using Gurobi

m = Model(solver=GurobiSolver())

@defVar(m, 5 <= x <= 5)
@defVar(m, y >= 0)
@defVar(m, 0.7 <= z1 <= 2.1, Int)
@defVar(m, 0.7 <= z2 <= 1.3, Int)
@defVar(m, w1 >= 0)
@defVar(m, w2 >= 0)

@setObjective(m, Max, y + x + z1 + z2 + 100)
@addConstraint(m, y + x <= 10)
@addConstraint(m, w1 + w2 <= 1)

solve(m)
println("Objective value ", getObjectiveValue(m))
println("presolve")
presolve(m)