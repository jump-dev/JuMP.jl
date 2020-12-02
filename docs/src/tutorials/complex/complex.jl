# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # A basic example of complex model

# This example uses the following packages:

using JuMP
import GLPK
import ComplexOptInterface
const COI = ComplexOptInterface
import Test  #src

# Build the model:

model = Model(GLPK.Optimizer)
add_bridge(model, COI.Bridges.Constraint.SplitZeroBridge)

@variable(model, 0 <= x[1:3] <= 1)

@objective(model, Max, x[1] + 2x[2] - x[3])

@constraint(
    model,
    [(1 + im) * x[1] + (1 + im) * x[2] + (1 - 3im) * x[3] - 1.0] in
    MOI.Zeros(1)
)

# Optimize the model:

optimize!(model)

# Check the termination and primal status to see if we have a solution:

println("Termination status : ", termination_status(model))
println("Primal status      : ", primal_status(model))

# Print the solution:

obj_value = objective_value(model)
x_values = value.(x)

println("Objective value : ", obj_value)
println("x values        : ", x_values)

Test.@test obj_value ≈ 1.25             #src
Test.@test x_values ≈ [0.0, 0.75, 0.25] #src
