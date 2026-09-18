# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Sensitivity analysis of a linear program

# This tutorial explains how to perform sensitivity analysis on a linear program
# using JuMP's [`lp_sensitivity_report`](@ref), producing tables similar to those in
# Excel Solver. Sensitivity analysis answers how much problem data can change
# before the current optimal basis—and therefore the current solution—
# changes.
#
# In brief, sensitivity analysis of a linear program is about asking two
# questions:
#  1) Given an optimal solution, how much can the objective coefficients change
#     by before a different solution becomes optimal?
#  2) Given an optimal solution, how much can the right-hand side of a linear
#     constraint change by before a different solution becomes optimal?
#
# **Learning intentions:**
# * Understand what sensitivity analysis tells you: how far an objective
#   coefficient or constraint right-hand side can move before the optimal basis
#   changes
# * Use [`lp_sensitivity_report`](@ref) to retrieve allowable ranges, and
#   organize them alongside reduced costs, shadow prices, and slacks into
#   readable DataFrames
# * Read the resulting tables to identify which variables are basic or
#   non-basic and which constraints are binding

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import HiGHS
import DataFrames

# ## Setup

# As an example, we use this small linear program:

model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@variable(model, z <= 1)
@objective(model, Min, 12x + 20y - z)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
@constraint(model, c3, x + y <= 20)
optimize!(model)
assert_is_solved_and_feasible(model)
solution_summary(model; verbose = true)

# Can you identify:
#  * The objective coefficient of each variable?
#  * The right-hand side of each constraint?
#  * The optimal primal and dual solutions?

# ## Sensitivity reports

# Now let's call [`lp_sensitivity_report`](@ref):

report = lp_sensitivity_report(model)

# It returns a [`SensitivityReport`](@ref) object, which maps:
#
#  - Every variable reference to a tuple `(d_lo, d_hi)::Tuple{Float64,Float64}`,
#    explaining how much the objective coefficient of the corresponding variable
#    can change by, such that the original basis remains optimal.
#  - Every constraint reference to a tuple `(d_lo, d_hi)::Tuple{Float64,Float64}`,
#    explaining how much the right-hand side of the corresponding constraint can
#    change by, such that the basis remains optimal.

# Both tuples are relative, rather than absolute. So, given an objective
# coefficient of `1.0` and a tuple `(-0.5, 0.5)`, the objective coefficient can
# range between `1.0 - 0.5` an `1.0 + 0.5`.

# For example:

report[x]

# indicates that the objective coefficient on `x`, that is, `12`, can
# decrease by `-0.333` or increase by `3.0` and the primal solution
# `(15, 1.25)` will remain optimal. In addition:

report[c1]

# means that the right-hand side of the `c1` constraint (100), can decrease
# by 4 units, or increase by 2.85 units, and the primal solution `(15, 1.25)`
# will remain optimal.

# ## Variable sensitivity

# By themselves, the tuples aren't informative. Let's put them in context by
# collating a range of other information about a variable:

function variable_report(xi)
    return (
        name = name(xi),
        lower_bound = has_lower_bound(xi) ? lower_bound(xi) : -Inf,
        value = value(xi),
        upper_bound = has_upper_bound(xi) ? upper_bound(xi) : Inf,
        reduced_cost = reduced_cost(xi),
        obj_coefficient = coefficient(objective_function(model), xi),
        allowed_decrease = report[xi][1],
        allowed_increase = report[xi][2],
    )
end

# Calling our function on `x`:

x_report = variable_report(x)

# That's a bit hard to read, so let's call this on every variable in the model
# and put things into a DataFrame:

variable_df =
    DataFrames.DataFrame(variable_report(xi) for xi in all_variables(model))

# ## Constraint sensitivity

# We can do something similar with constraints:

function constraint_report(c::ConstraintRef)
    return (
        name = name(c),
        value = value(c),
        rhs = normalized_rhs(c),
        slack = normalized_rhs(c) - value(c),
        shadow_price = shadow_price(c),
        allowed_decrease = report[c][1],
        allowed_increase = report[c][2],
    )
end

c1_report = constraint_report(c1)

# That's a bit hard to read, so let's call this on every variable in the model
# and put things into a DataFrame:

constraint_df = DataFrames.DataFrame(
    constraint_report(ci) for (F, S) in list_of_constraint_types(model) for
    ci in all_constraints(model, F, S) if F == AffExpr
)

# ## Analysis questions

# Now we can use these dataframes to ask questions of the solution.

# For example, we can find basic variables by looking for variables with a
# reduced cost of 0:

basic = filter(row -> iszero(row.reduced_cost), variable_df)

# and non-basic variables by looking for non-zero reduced costs:

non_basic = filter(row -> !iszero(row.reduced_cost), variable_df)

# we can also find constraints that are binding by looking for zero slacks:

binding = filter(row -> iszero(row.slack), constraint_df)

# or non-zero shadow prices:

binding2 = filter(row -> !iszero(row.shadow_price), constraint_df)
