# TODO: Copy over tests that are still relevant from old/macros.jl.

mutable struct MyVariable
    info::JuMP.VariableInfo
    test_kw::Int
end

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

    @testset "Extension of @variable with constructvariable! #1029" begin
        JuMP.variabletype(m::Model, ::Type{MyVariable}) = MyVariable
        function JuMP.constructvariable!(m::Model, ::Type{MyVariable}, _error::Function, args...; test_kw::Int = 0)
            MyVariable(args..., test_kw)
        end
        m = Model()
        @variable(m, 1 <= x <= 2, MyVariable, binary = true, test_kw = 1, start = 3)
        @test isa(x, MyVariable)
        @test x.info.haslb
        @test x.info.lowerbound == 1
        @test x.info.hasub
        @test x.info.upperbound == 2
        @test !x.info.hasfix
        @test isnan(x.info.fixedvalue)
        @test x.info.binary
        @test !x.info.integer
        @test x.info.hasstart
        @test x.info.start == 3
        @test x.info.name == "x"
        @test x.test_kw == 1
        @variable(m, y[1:3] >= 0, MyVariable, test_kw = 2)
        @test isa(y, Vector{MyVariable})
        for i in 1:3
            @test y[i].info.haslb
            @test y[i].info.lowerbound == 0
            @test !y[i].info.hasub
            @test y[i].info.upperbound == Inf
            @test !y[i].info.hasfix
            @test isnan(y[i].info.fixedvalue)
            @test !y[i].info.binary
            @test !y[i].info.integer
            @test !y[i].info.hasstart
            @test isnan(y[i].info.start)
            @test y[i].info.name == "y[$i]"
            @test y[i].test_kw == 2
        end
    end

end
