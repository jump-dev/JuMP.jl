# TODO: Copy over tests that are still relevant from old/macros.jl.

mutable struct MyVariable
    test_kw::Int
    info::JuMP.VariableInfo
end

@testset "Extension of @variable with buildvariable #1029" begin
    const MyVariable = Tuple{JuMP.VariableInfo, Int}
    JuMP.variabletype(m::Model, ::Type{MyVariable}) = MyVariable
    names = Dict{MyVariable, String}()
    function JuMP.addvariable(m::Model, v::MyVariable, name::String)
        names[v] = name
        v
    end
    function JuMP.buildvariable(_error::Function, info::JuMP.VariableInfo, ::Type{MyVariable}; test_kw::Int = 0)
        (info, test_kw)
    end
    m = Model()
    @variable(m, 1 <= x <= 2, MyVariable, binary = true, test_kw = 1, start = 3)
    @test isa(x, MyVariable)
    info = x[1]
    test_kw = x[2]
    @test info.lowerbound == 1
    @test info.hasub
    @test info.upperbound == 2
    @test !info.hasfix
    @test isnan(info.fixedvalue)
    @test info.binary
    @test !info.integer
    @test info.hasstart
    @test info.start == 3
    @test names[x] == "x"
    @test test_kw == 1
    @variable(m, y[1:3] >= 0, MyVariable, test_kw = 2)
    @test isa(y, Vector{MyVariable})
    for i in 1:3
        info = y[i][1]
        test_kw = y[i][2]
        @test info.haslb
        @test info.lowerbound == 0
        @test !info.hasub
        @test isnan(info.upperbound)
        @test !info.hasfix
        @test isnan(info.fixedvalue)
        @test !info.binary
        @test !info.integer
        @test !info.hasstart
        @test isnan(info.start)
        @test names[y[i]] == "y[$i]"
        @test test_kw == 2
    end
end

function macros_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}

    @testset "Nested tuple destructuring" begin
        m = ModelType()
        d = Dict((1,2) => 3)
        ex = @expression(m, sum(i+j+k for ((i,j),k) in d))
        @test ex == 6
    end

    @testset "buildconstraint on variable" begin
        m = ModelType()
        @variable(m, x)
        @test JuMP.buildconstraint(error, x, MOI.GreaterThan(0.0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.GreaterThan{Float64}}
        @test JuMP.buildconstraint(error, x, MOI.LessThan(0.0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.LessThan{Float64}}
        @test JuMP.buildconstraint(error, x, MOI.EqualTo(0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.EqualTo{Int}}
    end

    @testset "Check @constraint basics" begin
        m = ModelType()
        @variable(m, w)
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        t = 10.0

        cref = @constraint(m, 3x - y == 3.3(w + 2z) + 5)
        c = JuMP.constraintobject(cref, AffExprType, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, 3x - y == (w + 2z)*3.3 + 5)
        c = JuMP.constraintobject(cref, AffExprType, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, (x+y)/2 == 1)
        c = JuMP.constraintobject(cref, AffExprType, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 0.5*x + 0.5*y)
        @test c.set == MOI.EqualTo(1.0)

        cref = @constraint(m, -1 <= x-y <= t)
        c = JuMP.constraintobject(cref, AffExprType, MOI.Interval)
        @test JuMP.isequal_canonical(c.func, x - y)
        @test c.set == MOI.Interval(-1.0, t)

        cref = @constraint(m, -1 <= x+1 <= 1)
        c = JuMP.constraintobject(cref, AffExprType, MOI.Interval)
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-2.0, 0.0)

        cref = @constraint(m, -1 <= x <= 1)
        c = JuMP.constraintobject(cref, VariableRefType, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, -1 <= x <= sum(0.5 for i = 1:2))
        c = JuMP.constraintobject(cref, VariableRefType, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, 1 >= x >= 0)
        c = JuMP.constraintobject(cref, VariableRefType, MOI.Interval)
        @test c.func == x
        @test c.set == MOI.Interval(0.0, 1.0)

        @test_throws ErrorException @constraint(m, x <= t <= y)
        @test_throws ErrorException @constraint(m, 0 <= Dict() <= 1)
        @test_macro_throws ErrorException @constraint(1 <= x <= 2, foo=:bar)

        @test JuMP.isequal_canonical(@expression(m, 3x - y - 3.3(w + 2z) + 5), 3*x - y - 3.3*w - 6.6*z + 5)
        @test JuMP.isequal_canonical(@expression(m, quad, (w+3)*(2x+1)+10), 2*w*x + 6*x + w + 13)

        @test_throws ErrorException @constraint(m, 3 + 5*7 <= 0)
    end

    @testset "Helpful error for variable declaration with misplaced constant" begin
        m = ModelType()
        @test_macro_throws ErrorException @variable m 0 <= x
        @test_macro_throws ErrorException @variable m 0 >= x
    end
end

@testset "Macros for JuMP.Model" begin
    macros_test(Model, VariableRef{Model{JuMP.NonDirectBackendType}})
end

@testset "Macros for JuMPExtension.MyModel" begin
    macros_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
