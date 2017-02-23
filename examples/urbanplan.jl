#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# urbanplan.jl
#
# An "urban planning" problem.
# Based on
#  http://puzzlor.editme.com/urbanplanning
#############################################################################

using JuMP, Cbc

solver = CbcSolver()

function SolveUrban()

    m = Model(solver=solver)

    # x is indexed by row and column
    @variable(m, 0 <= x[1:5,1:5] <= 1, Int)

    # y is indexed by R or C, and the points
    # JuMP allows indexing on arbitrary sets
    rowcol = ["R","C"]
    points = [+5,+4,+3,-3,-4,-5]
    @variable(m, 0 <= y[rowcol,points,i=1:5] <= 1, Int)

    # Objective - combine the positive and negative parts
    @objective(m, Max, sum(
      3*(y["R", 3,i] + y["C", 3,i])
    + 1*(y["R", 4,i] + y["C", 4,i])
    + 1*(y["R", 5,i] + y["C", 5,i])
    - 3*(y["R",-3,i] + y["C",-3,i])
    - 1*(y["R",-4,i] + y["C",-4,i])
    - 1*(y["R",-5,i] + y["C",-5,i]) for i=1:5))

    # Constrain the number of residential lots
    @constraint(m, sum(x[i,j] for i=1:5, j=1:5) == 12)

    # Add the constraints that link the auxiliary y variables
    # to the x variables
    # Rows
    for i = 1:5
        @constraint(m, y["R", 5,i] <=   1/5*sum(x[i,j] for j=1:5)) # sum = 5
        @constraint(m, y["R", 4,i] <=   1/4*sum(x[i,j] for j=1:5)) # sum = 4
        @constraint(m, y["R", 3,i] <=   1/3*sum(x[i,j] for j=1:5)) # sum = 3
        @constraint(m, y["R",-3,i] >= 1-1/3*sum(x[i,j] for j=1:5)) # sum = 2
        @constraint(m, y["R",-4,i] >= 1-1/2*sum(x[i,j] for j=1:5)) # sum = 1
        @constraint(m, y["R",-5,i] >= 1-1/1*sum(x[i,j] for j=1:5)) # sum = 0
    end
    # Columns
    for j = 1:5
        @constraint(m, y["C", 5,j] <=   1/5*sum(x[i,j] for i=1:5)) # sum = 5
        @constraint(m, y["C", 4,j] <=   1/4*sum(x[i,j] for i=1:5)) # sum = 4
        @constraint(m, y["C", 3,j] <=   1/3*sum(x[i,j] for i=1:5)) # sum = 3
        @constraint(m, y["C",-3,j] >= 1-1/3*sum(x[i,j] for i=1:5)) # sum = 2
        @constraint(m, y["C",-4,j] >= 1-1/2*sum(x[i,j] for i=1:5)) # sum = 1
        @constraint(m, y["C",-5,j] >= 1-1/1*sum(x[i,j] for i=1:5)) # sum = 0
    end

    # Solve it with the default solver (CBC)
    status = solve(m)
    if status == :Infeasible
        error("Solver couldn't find solution!")
    end

    # Print results
    println("Best objective: $(round(getobjectivevalue(m)))")
    println(getvalue(x))
end

SolveUrban()
