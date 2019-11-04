#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

using JuMP, GLPK, Test
const MOI = JuMP.MathOptInterface

"""
    example_sudoku(filepath)

A sudoku solver that uses a MIP to find solutions.

We have binary variables x[i, j, k] which, if = 1, say that cell (i, j) contains
the number k. The constraints are:
 1 - Each cell has one value only
 2 - Each row contains each number exactly once
 3 - Each column contains each number exactly once
 4 - Each 3x3 subgrid contains each number exactly once
We will take the initial grid as a CSV file at `filepath`, where 0s are blanks.
"""
function example_sudoku(filepath)
    initial_grid = zeros(Int, 9, 9)
    open(filepath, "r") do fp
        for row in 1:9
            line = readline(fp)
            initial_grid[row, :] .= parse.(Int, split(line, ","))
        end
    end

    model = Model(GLPK.Optimizer)

    @variable(model, x[1:9, 1:9, 1:9], Bin)

    @constraints(model, begin
        # Constraint 1 - Only one value appears in each cell
        cell[i in 1:9, j in 1:9], sum(x[i, j, :]) == 1
        # Constraint 2 - Each value appears in each row once only
        row[i in 1:9, k in 1:9], sum(x[i, :, k]) == 1
        # Constraint 3 - Each value appears in each column once only
        col[j in 1:9, k in 1:9], sum(x[:, j, k]) == 1
        # Constraint 4 - Each value appears in each 3x3 subgrid once only
        subgrid[i=1:3:7, j=1:3:7, val=1:9], sum(x[i:i + 2, j:j + 2, val]) == 1
    end)

    # Initial solution
    for row in 1:9, col in 1:9
        if initial_grid[row, col] != 0
            @constraint(model, x[row, col, initial_grid[row, col]] == 1)
        end
    end

    # Solve it
    JuMP.optimize!(model)

    term_status = JuMP.termination_status(model)
    primal_status = JuMP.primal_status(model)
    is_optimal = term_status == MOI.OPTIMAL

    # Check solution
    if is_optimal
        mip_solution = JuMP.value.(x)
        sol = zeros(Int, 9, 9)
        for row in 1:9, col in 1:9, val in 1:9
            if mip_solution[row, col, val] >= 0.9
                sol[row, col] = val
            end
        end
        return sol
    else
        error("The solver did not find an optimal solution.")
    end
end

function print_sudoku_solution(solution)
    println("Solution:")
    println("[-----------------------]")
    for row in 1:9
        print("[ ")
        for col in 1:9
            print(solution[row, col], " ")
            if col % 3 == 0 && col < 9
                print("| ")
            end
        end
        println("]")
        if row % 3 == 0
            println("[-----------------------]")
        end
    end
end

if length(ARGS) == 1
    sol = example_sudoku(ARGS[1])
    print_sudoku_solution(sol)
else
    solution = example_sudoku(joinpath(@__DIR__, "data/sudoku.csv"))
    @test solution == [
        3 1 7 9 5 8 2 6 4;
        4 6 9 3 2 7 8 1 5;
        8 2 5 1 6 4 7 9 3;
        2 7 1 6 4 5 3 8 9;
        9 3 8 7 1 2 5 4 6;
        5 4 6 8 3 9 1 7 2;
        1 8 4 2 9 3 6 5 7;
        6 9 2 5 7 1 4 3 8;
        7 5 3 4 8 6 9 2 1]
end
