using JuMP
using CPLEX

m = Model(solver=CplexSolver())

@defVar(m, x[1:5])
@setObjective(m, Max, x[1] + x[2] + x[3] + x[4] + 2*x[5])
@addConstraint(m, constrcon1, norm(x[2:5]) <= x[1])
@addConstraint(m, constrlin1, x[1] <= 5)
@addConstraint(m, constrlin2, x[1]+2*x[2] - x[3] <= 10)
@addConstraint(m, constrcon2, norm([2 3;1 1]*x[2:3]-[3;4]) <= x[5] - 2)
#@addConstraint(m, constrcon3, x[5]^2 + x[4]^2 <= 10)
#@addConstraint(m, constrcon2, norm2{x[i], i = 2:5} <= 10)

solve(m)
@show MathProgBase.getconstrduals(m.internalModel)
@show getDual(constrcon1)
@show getDual(constrcon2)
@show m.conicconstrDuals
@show m.linconstrDuals
