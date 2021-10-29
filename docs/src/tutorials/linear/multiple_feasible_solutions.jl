# Copyright (c) 2021 James D Foster and contribut                                              #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Finding multiple feasible solutions 

# **Originally Contributed by**: James Foster

# This tutorial demonstrates how to formulate and solve a combinatorial problem
# with multiple solutions. In fact, we will show how to find _all_ feasible solutions to our problem.

# We'll look at some objects from recreational mathematics that could be called **"symmetric square sums"**.
# Here are a few examples (inspired by [Number Squares](https://www.futilitycloset.com/2012/12/05/number-squares/)):
"""
    1529      2318      5219
    5837      3790      2384
    2340      1956      1867
  + ____    + ____    + ____
    9706      8064      9470
"""

# Notice how all the digits 0 to 9 are used at least once,
# the first three rows sum to the last row,
# and the columns in each are the same as the corresponding rows (forming a symmetric matrix).

# You might ask: how many such squares are there?
# In this tutorial, we will find out
# and show how to go about answering this and similar questions.

# This tutorial uses the following packages:
using JuMP
import Gurobi

# Here we are using [Gurobi](https://github.com/jump-dev/Gurobi.jl) 
# since it provides the required functionality
# for this example (i.e. finding multiple feasible solutions).
# Gurobi is a commercial solver, that is, a paid license is needed
# for those using the solver for commercial purposes. However there are
# trial and/or free licenses available for academic and student users.

# ## Model Specifics

solver = Gurobi.Optimizer
model = JuMP.Model(solver)

# We need to set specific Gurobi parameters to enable the
# [multiple solution functionality](https://www.gurobi.com/documentation/9.0/refman/finding_multiple_solutions.html).

# The first setting turns on the exhaustive search mode for multiple solutions:
JuMP.set_optimizer_attribute(model, "PoolSearchMode", 2)

# The second sets a limit for the number of solutions found. 
# Here the value 100 is an arbitrary but large enough whole number
# for our particular model (but will depend on the application in general).
JuMP.set_optimizer_attribute(model, "PoolSolutions", 100)

# ## Setting up the model:

# We are going to uses 4-digit numbers:
number_of_digits = 4

# So we define the index sets for our variables and constraints.
# We keep track of each "place" (units, tens, one-hundreds, one-thousands):
PLACES = 0:(number_of_digits-1)
# The number of rows of the symmetric square sums are the same are the number of digits:
ROWS = 1:number_of_digits

## Let's define the model's core variables.
# Here a given digit between 0 and 9 is found in the `i`-th row at the `j`-th place:
@variable(model, 0 <= Digit[i = ROWS, j = PLACES] <= 9, Int)
# We also need a higher level "term" variable that represents the actual number in each row:
@variable(model, Term[ROWS] >= 1, Int)
# The lower bound of 1 is because we require non-zero solutions.

# Now define the objective function
# (since this is a feasibility problem, it can just be 0):
@objective(model, Max, 0)
# @objective(model, Max, Digit[1,0]) #src

# Make sure the leading digit of each row is not zero:
@constraint(model, NonZeroLead[i in ROWS], Digit[i, (number_of_digits-1)] >= 1)

# Define the terms from the digits:
@constraint(
    model,
    TermDef[i in ROWS],
    Term[i] == sum((10^j) * Digit[i, j] for j in PLACES)
)

# The sum of the first three terms equals the last term:
@constraint(
    model,
    SumHolds,
    Term[number_of_digits] == sum(Term[i] for i in 1:(number_of_digits-1))
)

# The square is symmetric, that is, the sum should work either row-wise or column-wise:
@constraint(
    model,
    Symmetry[i in ROWS, j in PLACES; i + j <= (number_of_digits - 1)],
    Digit[i, j] == Digit[number_of_digits-j, number_of_digits-i]
)

# We also want to make sure we use each digit exactly once on the diagonal or upper triangular region.
# The following set, with the collection of binary variables and constraints ensures this property.
COMPS = [
    (i, j, k, m) for i in ROWS for j in PLACES for k in ROWS for
    m in PLACES if (
        i + j <= number_of_digits &&
        k + m <= number_of_digits &&
        (i, j) < (k, m)
    )
]

@variable(model, BinDiffs[COMPS], Bin)

@constraint(
    model,
    AllDiffLo[(i, j, k, m) in COMPS],
    Digit[i, j] <= Digit[k, m] - 1 + 50 * BinDiffs[(i, j, k, m)]
)

@constraint(
    model,
    AllDiffHi[(i, j, k, m) in COMPS],
    Digit[i, j] >= Digit[k, m] + 1 - 50 * (1 - BinDiffs[(i, j, k, m)])
)

# Note that the constant 50 is a "big enough" number to make these valid constraints; see 
# [this paper](https://doi.org/10.1287/ijoc.13.2.96.10515) and 
# [blog](https://yetanothermathprogrammingconsultant.blogspot.com/2016/05/all-different-and-mixed-integer.html)
# for more information.

# We can then call `optimize!` and view the results.
optimize!(model)

# Let's check it worked:
@assert JuMP.termination_status(model) == MOI.OPTIMAL
@assert JuMP.primal_status(model) == MOI.FEASIBLE_POINT

JuMP.objective_value(model)
JuMP.value.(Digit,) |> show

# ### Viewing the Results

# Now that we have results, we can print out the feasible solutions found:
TermInt = Dict()
for i in 1:result_count(model)
    TermInt[i] = convert.(Int64, round.(value.(Term; result = i).data))
end

an_optimal_solution = display(TermInt[1])
optimal_objective = objective_value(model; result = 1)
for i in 1:result_count(model)
    @assert has_values(model; result = i)
    println("Solution $(i): ")
    display(TermInt[i])
    print("\n")
end

# The result is the full list of feasible solutions.
# The answer to "how many such squares are there?" turns out to be 20.
