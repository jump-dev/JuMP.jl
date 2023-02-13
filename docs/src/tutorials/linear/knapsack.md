```@meta
EditURL = "<unknown>/docs/src/tutorials/linear/knapsack.jl"
```

# The knapsack problem

Formulate and solve a simple knapsack problem:

    max sum(p_j x_j)
     st sum(w_j x_j) <= C
        x binary

````@example knapsack
using JuMP
import HiGHS
import Test

function example_knapsack()
    profit = [5, 3, 2, 7, 4]
    weight = [2, 8, 4, 2, 5]
    capacity = 10
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:5], Bin)
    # Objective: maximize profit
    @objective(model, Max, profit' * x)
    # Constraint: can carry all
    @constraint(model, weight' * x <= capacity)
    # Solve problem using MIP solver
    optimize!(model)
    println("Objective is: ", objective_value(model))
    println("Solution is:")
    for i in 1:5
        print("x[$i] = ", value(x[i]))
        println(", p[$i]/w[$i] = ", profit[i] / weight[i])
    end
    Test.@test termination_status(model) == OPTIMAL
    Test.@test primal_status(model) == FEASIBLE_POINT
    Test.@test objective_value(model) == 16.0
    return
end

example_knapsack()
````

---

!!! tip
    This tutorial was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl). [View the source `.jl` file on GitHub](https://github.com/jump-dev/JuMP.jl/tree/master/docs/src/tutorials/linear/knapsack.jl).
