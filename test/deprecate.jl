module TestDeprecate

using Test
using JuMP

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_NonlinearExpression()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sin(x))
    @test_logs (:warn,) expr.m === model
    return
end

function test_NonlinearParameter()
    model = Model()
    @variable(model, x)
    @NLparameter(model, p == 1)
    @test_logs (:warn,) p.m === model
    return
end

function test_NLPEvaluator()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sin(x))
    d = JuMP.NLPEvaluator(model)
    @test_logs (:warn,) d.m === model
    return
end

function test_value()
    model = Model()
    @variable(model, x)
    @test_logs (:warn,) value(x, i -> 1.0) == 1.0
    return
end

function test_Model_caching_mode()
    @test_logs (:warn,) Model(caching_mode = MOIU.MANUAL)
    return
end

function test_Model_solver()
    optimizer =
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    @test_throws ErrorException Model(solver = optimizer)
    return
end

function test_set_optimizer()
    model = Model()
    optimizer =
        () -> MOI.Utilities.MockOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        )
    @test_logs(
        (:warn,),
        set_optimizer(model, optimizer; bridge_constraints = false),
    )
    return
end

end

TestDeprecate.runtests()
