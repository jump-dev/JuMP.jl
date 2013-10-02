#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/IainNZ/JuMP.jl
#############################################################################
# qcp.jl
#
# A simple quadratically constrained probgram
# Based on http://www.gurobi.com/documentation/5.5/example-tour/node25
#############################################################################

using JuMP

# Only solver supported that can solve QCPs so far is Gurobi.
if Pkg.installed("Gurobi") == nothing
  error("Must have Gurobi available for quadratic constraints")
end

# Maximization problem
m = Model(:Max, lpsolver=LPSolver(:Gurobi))

# Need nonnegativity for (rotated) second-order cone
@defVar(m, x)
@defVar(m, y >= 0)
@defVar(m, z >= 0)

# Maximize x
@setObjective(m, x)

# Subject to 1 linear and 2 nonlinear constraints
addConstraint(m, x + y + z == 1)
addConstraint(m, x*x + y*y - z*z <= 0)
addConstraint(m, x*x - y*z <= 0)

# Print the model to check correctness
print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getObjectiveValue(m))
println("x = ", getValue(x))
println("y = ", getValue(y))
