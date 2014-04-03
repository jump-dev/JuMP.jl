using JuMP

m = Model()
@defSDPVar(m, X[3])

addConstraint(m, X >= ones(3,3))
setObjective(m, :Min, trace(ones(3,3)*X))

stat = solve(m)
