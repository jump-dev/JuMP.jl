using JuMP
using CPLEX

m = Model(solver=CplexSolver())

@defVar(m, x[1:5])
@setObjective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
@addConstraint(m, constr, sum{x[i]^2, i = 2:5} <= x[1])
@addConstraint(m, constr, x[1] <= 10)

solve(m)
@show getDual(constr)
@show m.conicconstrDuals
@show m.linconstrDuals
