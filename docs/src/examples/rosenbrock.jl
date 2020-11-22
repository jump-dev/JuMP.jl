# # The Rosenbrock function

# A nonlinear example of the classical Rosenbrock function.

using JuMP
import Ipopt
import Test

function example_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)

    Test.@test termination_status(model) == MOI.LOCALLY_SOLVED
    Test.@test primal_status(model) == MOI.FEASIBLE_POINT
    Test.@test objective_value(model) ≈ 0.0 atol = 1e-10
    Test.@test value(x) ≈ 1.0
    Test.@test value(y) ≈ 1.0
    return
end

example_rosenbrock()
