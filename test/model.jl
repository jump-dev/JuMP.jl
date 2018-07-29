@testset "optimizehook" begin
    m = Model()
    @test m.optimizehook === nothing
    called = false
    function hook(m)
        called = true
    end
    JuMP.setoptimizehook(m, hook)
    @test !called
    optimize(m)
    @test called
end
@testset "UniversalFallback" begin
    m = Model()
    MOI.set!(m, MOIT.UnknownModelAttribute(), 1)
    @test MOI.get(m, MOIT.UnknownModelAttribute()) == 1
end
# Simple LP model not supporting Interval
@MOIU.model LPModel () (EqualTo, GreaterThan, LessThan) () () (SingleVariable,) (ScalarAffineFunction,) () ()
@testset "Bridges" begin
    @testset "Automatic bridging" begin
        # optimizer not supporting Interval
        optimizer = MOIU.MockOptimizer(LPModel{Float64}());
        model = Model(optimizer=optimizer)
        @variable model x
        cref = @constraint model 0 <= x + 1 <= 1
        @test cref isa JuMP.ConstraintRef{JuMP.Model,MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}}}
        JuMP.optimize(model)
    end
    @testset "Automatic bridging disabled with `bridge_constraints` keyword" begin
        optimizer = MOIU.MockOptimizer(LPModel{Float64}());
        model = Model(optimizer=optimizer, bridge_constraints=false)
        @test model.moibackend isa MOIU.CachingOptimizer
        @test model.moibackend === JuMP.caching_optimizer(model)
        @variable model x
        @test_throws ErrorException @constraint model 0 <= x + 1 <= 1
    end
    @testset "No bridge automatically added in Direct mode" begin
        optimizer = MOIU.MockOptimizer(LPModel{Float64}());
        model = Model(backend=optimizer, mode=JuMP.Direct)
        @variable model x
        @test_throws ErrorException @constraint model 0 <= x + 1 <= 1
    end
end
