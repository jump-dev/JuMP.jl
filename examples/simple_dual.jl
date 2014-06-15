using JuMP

n = 2

m = Model()
@defVar(m, 0 <= y <= 4)

addConstraint(m, y*ones(n,n) >= 2*ones(n,n))
setObjective(m, :Max, 1.0*y)

stat = solve(m)
