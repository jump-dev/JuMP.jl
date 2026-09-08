# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # A basic example

# This example uses the following packages:

using JuMP
import GLPK
import Test  #src

# Build the model:

model = Model(GLPK.Optimizer)

@variable(model, 0 <= x <= 2)
@variable(model, 0 <= y <= 30)

@objective(model, Max, 5x + 3y)

@constraint(model, 1x + 5y <= 3.0)

print(model)

# Optimize the model:

optimize!(model)

# Check the termination and primal status to see if we have a solution:

println("Termination status : ", termination_status(model))
println("Primal status      : ", primal_status(model))

# Print the solution:

obj_value = objective_value(model)
x_value = value(x)
y_value = value(y)

println("Objective value : ", obj_value)
println("x value         : ", x_value)
println("y value         : ", y_value)

Test.@test obj_value ≈ 10.6  #src
Test.@test x_value ≈ 2       #src
Test.@test y_value ≈ 0.2     #src
