# # NLP: Rosenbrock

# A nonlinear example of the classical Rosenbrock function.

using JuMP, Ipopt, Test

function example_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @NLobjective(model, Min, (1 - x)^2 + 100 * (y - x^2)^2)
    optimize!(model)

    @test termination_status(model) == MOI.LOCALLY_SOLVED
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test objective_value(model) ≈ 0.0 atol = 1e-10
    @test value(x) ≈ 1.0
    @test value(y) ≈ 1.0
    return
end

example_rosenbrock()
