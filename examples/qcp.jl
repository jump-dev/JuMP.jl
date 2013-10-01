using JuMP

if Pkg.installed("Gurobi") == nothing
  error("Must have Gurobi available for quadratic constraints")
end

MathProgBase.setlpsolver(:Gurobi)

m = Model(:Max)

@defVar(m, x)
@defVar(m, y)
@defVar(m, z)

@setObjective(m, x)
# addConstraint(m, x + y + z == 1)
addConstraint(m, x*x + y*y - z*z <= 0)
addConstraint(m, x*x - y*z <= 0)

print(m)

status = solve(m)

println("Objective value: ", getObjectiveValue(m))
println("x = ", getValue(x))
println("y = ", getValue(y))
