#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# qcp.jl
#
# A simple quadratically constrained probgram
# Based on http://www.gurobi.com/documentation/5.5/example-tour/node25
#############################################################################

using JuMP, Gurobi

# Will require either Gurobi.jl, CPLEX.jl, or Mosek.jl to run
m = Model(solver=GurobiSolver())

# Need nonnegativity for (rotated) second-order cone
@variable(m, x)
@variable(m, y >= 0)
@variable(m, z >= 0)

# Maximize x
@objective(m, Max, x)

# Subject to 1 linear and 2 nonlinear constraints
@constraint(m, x + y + z == 1)
@constraint(m, x*x + y*y - z*z <= 0)
@constraint(m, x*x - y*z <= 0)

# Print the model to check correctness
print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))
