# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Computing the duals of a mixed-integer program

# This tutorial explains how to compute the duals of a mixed-integer linear
# program by fixing the discrete variables to their optimal solution and
# resolving as a linear program.

# This tutorial uses the following packages:

using JuMP
import HiGHS

# ## The model

# Our example model is the unit commitment example from [Unit commitment](@ref).
# The details are unimportant, other than to note that there are two types of
# continuous variables, `g` and `w`, representing the quantity of generation
# from thermal and wind plants, and a discrete variable `dispatch`, which is
# `1` if plant `i` is operating, and `0` if not.

# We are interested in the "dual" of the `power_balance` constraint, because it
# represents the marginal price of electricity that consumers should pay for
# their consumption.

generators = [
    (min = 0.0, max = 1000.0, fixed_cost = 1000.0, variable_cost = 50.0),
    (min = 300.0, max = 1000.0, fixed_cost = 0.0, variable_cost = 100.0),
]
N = length(generators)
model = Model(HiGHS.Optimizer)
set_silent(model)
@variables(model, begin
    generators[i].min <= g[i = 1:N] <= generators[i].max
    0 <= w <= 200
    dispatch[i = 1:N], Bin
end)
@constraints(model, begin
    power_balance, sum(g[i] for i in 1:N) + w == 1500
    [i = 1:N], g[i] <= generators[i].max * dispatch[i]
    [i = 1:N], g[i] >= generators[i].min * dispatch[i]
end)
@objective(
    model,
    Min,
    sum(
        generators[i].fixed_cost * dispatch[i] +
        generators[i].variable_cost * g[i] for i in 1:N
    )
)
print(model)

# ## Manually fix the variables

# If we optimize this model, we obtain a [`dual_status`](@ref) of `NO_SOLUTION`:

optimize!(model)
dual_status(model)

# This is because HiGHS cannot compute the duals of a mixed-integer program. We
# can work around this problem by fixing the integer variables to their optimal
# solution, relaxing integrality, and re-solving as a linear program.

discrete_values = value.(dispatch)
fix.(dispatch, discrete_values; force = true)
unset_binary.(dispatch)
print(model)

# Now if we re-solve the problem, we obtain a `FEASIBLE_POINT` for the dual:

optimize!(model)
dual_status(model)

# and a marginal price of electricity of \$100/MWh:

dual(power_balance)

# To reset the problem back to a mixed-integer linear program, we need to
# [`unfix`](@ref) and call [`set_binary`](@ref):

unfix.(dispatch)
set_binary.(dispatch)
print(model)

# ## Use `fix_discrete_variables`

# Manually choosing the variables to relax and fix works for our small example,
# but it becomes more difficult in problems with a larger number of binary and
# integer variables. To automate the process we just did manually, JuMP provides
# the [`fix_discrete_variables`](@ref) function:

optimize!(model)
dual_status(model)

#-

undo = fix_discrete_variables(model);

# Here `undo` is a function that, when called with no arguments, returns the
# model to the original mixed-integer formulation.

# !!! tip
#     After calling [`fix_discrete_variables`](@ref), you can set a new solver
#     with [`set_optimizer`](@ref) if your mixed-integer solver does not support
#     computing a dual solution.

print(model)

#-

optimize!(model)
dual_status(model)

#-

dual(power_balance)

# Finally, call `undo` to revert the reformulation

undo()
print(model)
