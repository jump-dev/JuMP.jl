# TODO: Copy over tests that are still relevant from old/macros.jl.

@testset "Macros" begin
    @testset "Nested tuple destructuring" begin
        m = Model()
        d = Dict((1,2) => 3)
        ex = @expression(m, sum(i+j+k for ((i,j),k) in d))
        @test ex == 6
    end

    @testset "constructconstraint! on variable" begin
        m = Model()
        @variable(m, x)
        @test JuMP.constructconstraint!(x, MOI.GreaterThan(0.0)) isa JuMP.SingleVariableConstraint{MOI.GreaterThan{Float64}}
        @test JuMP.constructconstraint!(x, MOI.LessThan(0.0)) isa JuMP.SingleVariableConstraint{MOI.LessThan{Float64}}
        @test JuMP.constructconstraint!(x, MOI.EqualTo(0)) isa JuMP.SingleVariableConstraint{MOI.EqualTo{Int}}
    end
end
