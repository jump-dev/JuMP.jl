# # MIP: knapsack

# Formulate and solve a simple knapsack problem:
#
#     max sum(p_j x_j)
#      st sum(w_j x_j) <= C
#         x binary

using JuMP, GLPK, Test

function example_knapsack(; verbose = true)
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:5], Bin)
    ## Objective: maximize profit
    @objective(model, Max, profit' * x)
    ## Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    ## Solve problem using MIP solver
    optimize!(model)
    if verbose
        println("Objective is: ", objective_value(model))
        println("Solution is:")
        for i in 1:5
            print("x[$i] = ", value(x[i]))
            println(", p[$i]/w[$i] = ", profit[i] / weight[i])
        end
    end
    @test termination_status(model) == MOI.OPTIMAL
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test objective_value(model) == 16.0
    return
end

example_knapsack()
