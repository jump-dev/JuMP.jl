# Copyright (c) 2021 James D Foster, and contributors                            #src
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

# _Author: James Foster (@jd-foster)_

# This tutorial demonstrates how to formulate and solve a combinatorial problem
# with multiple feasible solutions. In fact, we will see how to find _all_
# feasible solutions to our problem. We will also see how to enforce an
# "all-different" constraint on a set of integer variables.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import Gurobi
import Test

# !!! warning
#     This tutorial uses [Gurobi.jl](@ref) as the solver because it supports
#     returning multiple feasible solutions, something that open-source MIP
#     solvers such as HiGHS do not currently support. Gurobi is a commercial
#     solver and requires a paid license. However, there are free licenses
#     available for academic and student users. See [Gurobi.jl](@ref) for more
#     details.

# ## Symmetric number squares

# Symmetric [number squares](https://www.futilitycloset.com/2012/12/05/number-squares/)
# and their sums often arise in recreational mathematics. Here are a few
# examples:
# ```
#    1 5 2 9       2 3 1 8        5 2 1 9
#    5 8 3 7       3 7 9 0        2 3 8 4
# +  2 3 4 0     + 1 9 5 6      + 1 8 6 7
# =  9 7 0 6     = 8 0 6 4      = 9 4 7 0
# ```

# Notice how all the digits 0 to 9 are used at least once, the first three rows
# sum to the last row, the columns in each are the same as the corresponding
# rows (forming a symmetric matrix), and `0` does not appear in the first
# column.

# We will answer the question: how many such squares are there?

# ## JuMP model

# We now encode the symmetric number square as a JuMP model. First, we need a
# symmetric matrix of decision variables between `0` and `9` to represent each
# number:

n = 4
model = Model()
set_silent(model)
@variable(model, 0 <= x_digits[row in 1:n, col in 1:n] <= 9, Int, Symmetric)

# We modify the lower bound to ensure that the first column cannot contain `0`:

set_lower_bound.(x_digits[:, 1], 1)

# Then, we need a constraint that the sum of the first three rows equals the
# last row:

@expression(model, x_base_10, x_digits * [1_000, 100, 10, 1]);
@constraint(model, sum(x_base_10[i] for i in 1:n-1) == x_base_10[n])

# And we use [`MOI.AllDifferent`](@ref) to ensure that each digit is used
# exactly once in the upper triangle matrix of `x_digits`:

x_digits_upper = [x_digits[i, j] for j in 1:n for i in 1:j]
@constraint(model, x_digits_upper in MOI.AllDifferent(length(x_digits_upper)));

# If we optimize this model, we find that Gurobi has returned one solution:

set_optimizer(model, Gurobi.Optimizer)
optimize!(model)
Test.@test is_solved_and_feasible(model)
Test.@test result_count(model) == 1
solution_summary(model)

# We need to set specific Gurobi parameters to enable the
# [multiple solution functionality](https://www.gurobi.com/documentation/9.0/refman/finding_multiple_solutions.html).

# The first setting turns on the exhaustive search mode for multiple solutions:

set_optimizer(model, Gurobi.Optimizer)
set_attribute(model, "PoolSearchMode", 2)

# The second sets a limit for the number of solutions found:

set_attribute(model, "PoolSolutions", 100)

# Here the value 100 is an "arbitrary but large enough" whole number
# for our particular model (and in general will depend on the application).

# We can then call `optimize!` and view the results.

optimize!(model)
Test.@test is_solved_and_feasible(model)
solution_summary(model)

# Now Gurobi has found 20 solutions:

Test.@test result_count(model) == 20

# ### Viewing the Results

# Access the various feasible solutions by using the [`value`](@ref) function
# with the `result` keyword:

solutions =
    [round.(Int, value.(x_digits; result = i)) for i in 1:result_count(model)];

# Here we have converted the solution to an integer after rounding off very
# small numerical tolerances.

# An example of one feasible solution is:

solutions[1]

# and we can nicely print out all the feasible solutions with

function solution_string(x::Matrix)
    header = [" ", " ", "+", "="]
    return join([join(vcat(header[i], x[i, :]), " ") for i in 1:4], "\n")
end

for i in 1:result_count(model)
    println("Solution $i: \n", solution_string(solutions[i]), "\n")
end

# The result is the full list of feasible solutions. So the answer to "how many
# such squares are there?" turns out to be 20.