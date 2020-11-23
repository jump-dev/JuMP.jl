# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # The urban planning problem

# An "urban planning" problem based on an [example from puzzlor](http://www.puzzlor.com/2013-08_UrbanPlanning.html).

using JuMP
import GLPK
import Test

function example_urban_plan()
    model = Model(GLPK.Optimizer)
    ## x is indexed by row and column
    @variable(model, 0 <= x[1:5, 1:5] <= 1, Int)
    ## y is indexed by R or C, the points, and an index in 1:5. Note how JuMP
    ## allows indexing on arbitrary sets.
    rowcol = ["R", "C"]
    points = [5, 4, 3, -3, -4, -5]
    @variable(model, 0 <= y[rowcol, points, 1:5] <= 1, Int)
    ## Objective - combine the positive and negative parts
    @objective(model, Max, sum(
          3 * (y["R", 3, i] + y["C", 3, i])
        + 1 * (y["R", 4, i] + y["C", 4, i])
        + 1 * (y["R", 5, i] + y["C", 5, i])
        - 3 * (y["R", -3, i] + y["C", -3, i])
        - 1 * (y["R", -4, i] + y["C", -4, i])
        - 1 * (y["R", -5, i] + y["C", -5, i])
        for i in 1:5)
    )
    ## Constrain the number of residential lots
    @constraint(model, sum(x) == 12)
    ## Add the constraints that link the auxiliary y variables to the x variables
    for i = 1:5
        @constraints(model, begin
            ## Rows
            y["R", 5, i] <= 1 / 5 * sum(x[i, :]) # sum = 5
            y["R", 4, i] <= 1 / 4 * sum(x[i, :]) # sum = 4
            y["R", 3, i] <= 1 / 3 * sum(x[i, :]) # sum = 3
            y["R", -3, i] >= 1 - 1 / 3 * sum(x[i, :]) # sum = 2
            y["R", -4, i] >= 1 - 1 / 2 * sum(x[i, :]) # sum = 1
            y["R", -5, i] >= 1 - 1 / 1 * sum(x[i, :]) # sum = 0
            ## Columns
            y["C", 5, i] <= 1 / 5 * sum(x[:, i]) # sum = 5
            y["C", 4, i] <= 1 / 4 * sum(x[:, i]) # sum = 4
            y["C", 3, i] <= 1 / 3 * sum(x[:, i]) # sum = 3
            y["C", -3, i] >= 1 - 1 / 3 * sum(x[:, i]) # sum = 2
            y["C", -4, i] >= 1 - 1 / 2 * sum(x[:, i]) # sum = 1
            y["C", -5, i] >= 1 - 1 / 1 * sum(x[:, i]) # sum = 0
        end)
    end
    ## Solve it
    optimize!(model)
    Test.@test termination_status(model) == MOI.OPTIMAL
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 14.0
    return
end

example_urban_plan()
