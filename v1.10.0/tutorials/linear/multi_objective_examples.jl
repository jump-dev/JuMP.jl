# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Simple multi-objective examples

# This tutorial contains a number of examples of multi-objective programs from
# the literature.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import HiGHS
import MultiObjectiveAlgorithms as MOA

# ## Bi-objective linear problem

# This is example is taken from Example 6.3 (from Steuer, 1985), page 154 of
# Ehrgott, M. (2005). _Multicriteria Optimization_. Springer, Berlin. The code
# was adapted from an example in [vOptGeneric](https://github.com/vOptSolver/vOptGeneric.jl)
# by [@xgandibleux](https://github.com/xgandibleux).

model = Model()
set_silent(model)
@variable(model, x1 >= 0)
@variable(model, 0 <= x2 <= 3)
@objective(model, Min, [3x1 + x2, -x1 - 2x2])
@constraint(model, 3x1 - x2 <= 6)
set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_attribute(
    model,
    MOA.Algorithm(),
    MOA.Lexicographic(; all_permutations = true),
)
optimize!(model)
solution_summary(model)

#-

for i in 1:result_count(model)
    print(i, ": z = ", round.(Int, objective_value(model; result = i)), " | ")
    println("x = ", value.([x1, x2]; result = i))
end

# ## Bi-objective linear assignment problem

# This is example is taken from Example 9.38 (from Ulungu and Teghem, 1994),
# page 255 of Ehrgott, M. (2005). _Multicriteria Optimization_. Springer, Berlin.
# The code was adapted from an example in [vOptGeneric](https://github.com/vOptSolver/vOptGeneric.jl)
# by [@xgandibleux](https://github.com/xgandibleux).

C1 = [5 1 4 7; 6 2 2 6; 2 8 4 4; 3 5 7 1]
C2 = [3 6 4 2; 1 3 8 3; 5 2 2 3; 4 2 3 5]
n = size(C2, 1)
model = Model()
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@objective(model, Min, [sum(C1 .* x), sum(C2 .* x)])
@constraint(model, [i = 1:n], sum(x[i, :]) == 1)
@constraint(model, [j = 1:n], sum(x[:, j]) == 1)
set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
optimize!(model)
solution_summary(model)

#-

for i in 1:result_count(model)
    print(i, ": z = ", round.(Int, objective_value(model; result = i)), " | ")
    println("x = ", round.(Int, value.(x; result = i)))
end

# ## Bi-objective shortest path problem

# This is example is taken from Exercise 9.5 page 269 of Ehrgott, M. (2005).
# _Multicriteria Optimization_. Springer, Berlin. The code was adapted from an
# example in [vOptGeneric](https://github.com/vOptSolver/vOptGeneric.jl) by
# [@xgandibleux](https://github.com/xgandibleux).

M = 50
C1 = [
    M 4 5 M M M
    M M 2 1 2 7
    M M M 5 2 M
    M M 5 M M 3
    M M M M M 4
    M M M M M M
]
C2 = [
    M 3 1 M M M
    M M 1 4 2 2
    M M M 1 7 M
    M M 1 M M 2
    M M M M M 2
    M M M M M M
]
n = size(C2, 1)
model = Model()
set_silent(model)
@variable(model, x[1:n, 1:n], Bin)
@objective(model, Min, [sum(C1 .* x), sum(C2 .* x)])
@constraint(model, sum(x[1, :]) == 1)
@constraint(model, sum(x[:, n]) == 1)
@constraint(model, [i = 2:n-1], sum(x[i, :]) - sum(x[:, i]) == 0)
set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
optimize!(model)
solution_summary(model)

#-

for i in 1:result_count(model)
    print(i, ": z = ", round.(Int, objective_value(model; result = i)), " | ")
    X = round.(Int, value.(x; result = i))
    print("Path:")
    for ind in findall(val -> val ≈ 1, X)
        i, j = ind.I
        print(" $i->$j")
    end
    println()
end
