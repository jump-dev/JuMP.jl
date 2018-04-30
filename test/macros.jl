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

    @testset "buildconstraint on variable" begin
        m = Model()
        @variable(m, x)
        @test JuMP.buildconstraint(error, x, MOI.GreaterThan(0.0)) isa JuMP.SingleVariableConstraint{MOI.GreaterThan{Float64}}
        @test JuMP.buildconstraint(error, x, MOI.LessThan(0.0)) isa JuMP.SingleVariableConstraint{MOI.LessThan{Float64}}
        @test JuMP.buildconstraint(error, x, MOI.EqualTo(0)) isa JuMP.SingleVariableConstraint{MOI.EqualTo{Int}}
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

    @testset "Check @constraint basics" begin
        m = Model()
        @variable(m, w)
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        t = 10.0

        cref = @constraint(m, 3x - y == 3.3(w + 2z) + 5)
        c = JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, 3x - y == (w + 2z)*3.3 + 5)
        c = JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, (x+y)/2 == 1)
        c = JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 0.5*x + 0.5*y)
        @test c.set == MOI.EqualTo(1.0)

        cref = @constraint(m, -1 <= x-y <= t)
        c = JuMP.constraintobject(cref, AffExpr, MOI.Interval)
        @test JuMP.isequal_canonical(c.func, x - y)
        @test c.set == MOI.Interval(-1.0, t)

        cref = @constraint(m, -1 <= x+1 <= 1)
        c = JuMP.constraintobject(cref, AffExpr, MOI.Interval)
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-2.0, 0.0)

        cref = @constraint(m, -1 <= x <= 1)
        c = JuMP.constraintobject(cref, Variable, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, -1 <= x <= sum(0.5 for i = 1:2))
        c = JuMP.constraintobject(cref, Variable, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, 1 >= x >= 0)
        c = JuMP.constraintobject(cref, Variable, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(0.0, 1.0)

        @test_throws ErrorException @constraint(m, x <= t <= y)
        @test_throws ErrorException @constraint(m, 0 <= nothing <= 1)
        @test_macro_throws @constraint(1 <= x <= 2, foo=:bar)

        @test JuMP.isequal_canonical(@expression(m, 3x - y - 3.3(w + 2z) + 5), 3*x - y - 3.3*w - 6.6*z + 5)
        @test JuMP.isequal_canonical(@expression(m, quad, (w+3)*(2x+1)+10), 2*w*x + 6*x + w + 13)

        cref = @constraint(m, 3 + 5*7 <= 0)
        c = JuMP.constraintobject(cref, AffExpr, MOI.LessThan)
        @test JuMP.isequal_canonical(c.func, zero(AffExpr))
        @test c.set == MOI.LessThan(-38.0)
    end
end
