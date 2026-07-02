# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Primal and dual warm-starts

# Some conic solvers have the ability to set warm-starts for the primal and dual
# solution. This can improve performance, particularly if you are repeatedly
# solving a sequence of related problems.

# !!! tip
#     See [`set_start_values`](@ref) for a generic implementation of this
#     function that was added to JuMP after this tutorial was written.

# In this tutorial, we demonstrate how to write a function that sets the primal
# and dual starts as the optimal solution stored in a model. It is intended to
# be a starting point for which you can modify if you want to do something
# similar in your own code.

# !!! warning
#     This tutorial does not set start values for nonlinear models.

# This tutorial uses the following packages:

using JuMP
import SCS

# The main component of this tutorial is the following function. The most
# important observation is that we cache all of the solution values first, and
# then we modify the model second. (Alternating between querying a value and
# modifying the model is not allowed in JuMP.)

function set_optimal_start_values(model::Model)
    ## Store a mapping of the variable primal solution
    variable_primal = Dict(x => value(x) for x in all_variables(model))
    ## In the following, we loop through every constraint and store a mapping
    ## from the constraint index to a tuple containing the primal and dual
    ## solutions.
    constraint_solution = Dict()
    for (F, S) in list_of_constraint_types(model)
        ## We add a try-catch here because some constraint types might not
        ## support getting the primal or dual solution.
        try
            for ci in all_constraints(model, F, S)
                constraint_solution[ci] = (value(ci), dual(ci))
            end
        catch
            @info("Something went wrong getting $F-in-$S. Skipping")
        end
    end
    ## Now we can loop through our cached solutions and set the starting values.
    for (x, primal_start) in variable_primal
        set_start_value(x, primal_start)
    end
    for (ci, (primal_start, dual_start)) in constraint_solution
        set_start_value(ci, primal_start)
        set_dual_start_value(ci, dual_start)
    end
    return
end

# To test our function, we use the following linear program:

model = Model(SCS.Optimizer)
@variable(model, x[1:3] >= 0)
@constraint(model, sum(x) <= 1)
@objective(model, Max, sum(i * x[i] for i in 1:3))
optimize!(model)

# By looking at the log (not shown in Documenter due to a bug), we can see that
# SCS took 100 iterations to find the optimal solution. Now we set the optimal
# solution as our starting point:

set_optimal_start_values(model)

# and we re-optimize:

optimize!(model)

# Now the optimization terminates after 0 iterations because our starting point
# is already optimal.

# Note that some solvers do not support setting some parts of the starting
# solution, for example, they may support only `set_start_value` for variables.
# If you encounter an `UnsupportedSupported` attribute error for
# [`MOI.VariablePrimalStart`](@ref), [`MOI.ConstraintPrimalStart`](@ref), or
# [`MOI.ConstraintDualStart`](@ref), comment out the corresponding part of the
# `set_optimal_start_values` function.
