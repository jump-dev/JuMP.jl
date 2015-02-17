#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# simpleusercut.jl
#
# Solve a simple problem using user cuts. See documentation for more
# details about callbacks.
#############################################################################

using JuMP
using Gurobi

# We will use Gurobi, which requires that we manually set the attribute
# PreCrush to 1 if we have user cuts. We will also disable PreSolve, Cuts,
# and Heuristics so only our cut will be used
m = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0))

# Define our variables to be inside a box, and integer
@defVar(m, 0 <= x <= 2, Int)
@defVar(m, 0 <= y <= 2, Int)

# Optimal solution is trying to go towards top-right corner (2.0, 2.0)
@setObjective(m, Max, x + 2y)

# We have one constraint that cuts off the top right corner
@addConstraint(m, y + x <= 3.5)

# Optimal solution of relaxed problem will be (1.5, 2.0)
# We can add a user cut that will cut of this fractional solution.

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
function mycutgenerator(cb)
    x_val = getValue(x)
    y_val = getValue(y)
    println("In callback function, x=$x_val, y=$y_val")

    # Allow for some impreciseness in the solution
    TOL = 1e-6

    # Check top right
    if y_val + x_val > 3 + TOL
        # Cut off this solution
        println("Fractional solution was in top right, cut it off")
        # Use the original variables
        @addUserCut(cb, y + x <= 3)
    end
end  # End of callback function

# Tell JuMP/Gurobi to use our callback function
addCutCallback(m, mycutgenerator)

# Solve the problem
solve(m)

# Print our final solution
println("Final solution: [ $(getValue(x)), $(getValue(y)) ]")
