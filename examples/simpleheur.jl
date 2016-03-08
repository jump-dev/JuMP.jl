#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# simpleheur.jl
#
# Solve a simple problem using user heuristic. See documentation for more
# details about callbacks.
#############################################################################

using JuMP
using Gurobi

# We will use Gurobi and disable PreSolve, Cuts, and (in-built) Heuristics so
# only our heuristic will be used
m = Model(solver=GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0))

# Define our variables to be inside a box, and integer
@defVar(m, 0 <= x <= 2, Int)
@defVar(m, 0 <= y <= 2, Int)

# Optimal solution is trying to go towards top-right corner (2.0, 2.0)
@setObjective(m, Max, x + 2y)

# We have one constraint that cuts off the top right corner
@addConstraint(m, y + x <= 3.5)

# Optimal solution of relaxed problem will be (1.5, 2.0)

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
function myheuristic(cb)
    x_val = getValue(x)
    y_val = getValue(y)
    println("In callback function, x=$x_val, y=$y_val")

    setSolutionValue!(cb, x, floor(x_val))
    # Leave y undefined - solver should handle as it sees fit. In the case
    # of Gurobi it will try to figure out what it should be.
    addSolution(cb)

    # Submit a second solution
    setSolutionValue!(cb, x, ceil(x_val))
    addSolution(cb)
end  # End of callback function

# Tell JuMP/Gurobi to use our callback function
addHeuristicCallback(m, myheuristic)

# Solve the problem
solve(m)

# Print our final solution
println("Final solution: [ $(getValue(x)), $(getValue(y)) ]")
