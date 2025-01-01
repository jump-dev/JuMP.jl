# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Arbitrary precision arithmetic

# The purpose of this tutorial is to explain how to use a solver which supports
# arithmetic using a number type other than `Float64`.

# This tutorial uses the following packages:

using JuMP
import CDDLib
import Clarabel

# ## Higher-precision arithmetic

# To create a model with a number type other than `Float64`, use [`GenericModel`](@ref)
# with an optimizer which supports the same number type:

model = GenericModel{BigFloat}(Clarabel.Optimizer{BigFloat})

# The syntax for adding decision variables is the same as a normal JuMP
# model, except that values are converted to `BigFloat`:

@variable(model, -1 <= x[1:2] <= sqrt(big"2"))

# Note that each `x` is now a [`GenericVariableRef{BigFloat}`](@ref), which
# means that the value of `x` in a solution will be a `BigFloat`.

# The lower and upper bounds of the decision variables are also `BigFloat`:

lower_bound(x[1])

#-

typeof(lower_bound(x[1]))

#-

upper_bound(x[2])

#-

typeof(upper_bound(x[2]))

# The syntax for adding constraints is the same as a normal JuMP model, except
# that coefficients are converted to `BigFloat`:

@constraint(model, c, x[1] == big"2" * x[2])

# The function is a [`GenericAffExpr`](@ref) with `BigFloat` for the
# coefficient and variable types;

constraint = constraint_object(c)
typeof(constraint.func)

# and the set is a [`MOI.EqualTo{BigFloat}`](@ref):

typeof(constraint.set)

# The syntax for adding and objective is the same as a normal JuMP model, except
# that coefficients are converted to `BigFloat`:

@objective(model, Min, 3x[1]^2 + 2x[2]^2 - x[1] - big"4" * x[2])

#-

typeof(objective_function(model))

# Here's the model we have built:

print(model)

# Let's solve and inspect the solution:

optimize!(model)
@assert is_solved_and_feasible(model; dual = true)
solution_summary(model)

# The value of each decision variable is a `BigFloat`:

value.(x)

# as well as other solution attributes like the objective value:

objective_value(model)

# and dual solution:

dual(c)

# This problem has an analytic solution of `x = [3//7, 3//14]`. Currently, our
# solution has an error of approximately `1e-9`:

value.(x) .- [3 // 7, 3 // 14]

# But by reducing the tolerances, we can obtain a more accurate solution:

set_attribute(model, "tol_gap_abs", 1e-32)
set_attribute(model, "tol_gap_rel", 1e-32)
optimize!(model)
@assert is_solved_and_feasible(model)
value.(x) .- [3 // 7, 3 // 14]

# ## Rational arithmetic

# In addition to higher-precision floating point number types like `BigFloat`,
# JuMP also supports solvers with exact rational arithmetic. One example is
# CDDLib.jl, which supports the `Rational{BigInt}` number type:

model = GenericModel{Rational{BigInt}}(CDDLib.Optimizer{Rational{BigInt}})

# As before, we can create variables using rational bounds:

@variable(model, 1 // 7 <= x[1:2] <= 2 // 3)

#-

lower_bound(x[1])

#-

typeof(lower_bound(x[1]))

# As well as constraints:

@constraint(model, c1, (2 // 1) * x[1] + x[2] <= 1)

#-

@constraint(model, c2, x[1] + 3x[2] <= 9 // 4)

# and objective functions:

@objective(model, Max, sum(x))

# Here's the model we have built:

print(model)

# Let's solve and inspect the solution:

optimize!(model)
@assert is_solved_and_feasible(model)
solution_summary(model)

# The optimal values are given in exact rational arithmetic:

value.(x)

#-

objective_value(model)

#-

value(c2)
