# TODO: Copy over tests that are still relevant from old/macros.jl.

mutable struct MyVariable
    test_kw::Int
    info::JuMP.VariableInfo
end

@testset "Extension of @variable with build_variable #1029" begin
    local MyVariable = Tuple{JuMP.VariableInfo, Int}
    JuMP.variable_type(m::Model, ::Type{MyVariable}) = MyVariable
    names = Dict{MyVariable, String}()
    function JuMP.add_variable(m::Model, v::MyVariable, name::String="")
        names[v] = name
        v
    end
    function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, ::Type{MyVariable}; test_kw::Int = 0)
        (info, test_kw)
    end
    m = Model()
    @variable(m, 1 <= x <= 2, MyVariable, binary = true, test_kw = 1, start = 3)
    @test isa(x, MyVariable)
    info = x[1]
    test_kw = x[2]
    @test info.has_lb
    @test info.lower_bound == 1
    @test info.has_ub
    @test info.upper_bound == 2
    @test !info.has_fix
    @test isnan(info.fixed_value)
    @test info.binary
    @test !info.integer
    @test info.has_start
    @test info.start == 3
    @test names[x] == "x"
    @test test_kw == 1

    @variable(m, y[1:3] >= 0, MyVariable, test_kw = 2)
    @test isa(y, Vector{MyVariable})
    for i in 1:3
        info = y[i][1]
        test_kw = y[i][2]
        @test info.has_lb
        @test info.lower_bound == 0
        @test !info.has_ub
        @test isnan(info.upper_bound)
        @test !info.has_fix
        @test isnan(info.fixed_value)
        @test !info.binary
        @test !info.integer
        @test !info.has_start
        @test isnan(info.start)
        @test names[y[i]] == "y[$i]"
        @test test_kw == 2
    end

    z = @variable(m, variable_type=MyVariable, upper_bound=3, test_kw=5)
    info = z[1]
    test_kw = z[2]
    @test isa(z, MyVariable)
    @test !info.has_lb
    @test isnan(info.lower_bound)
    @test info.has_ub
    @test info.upper_bound == 3
    @test !info.has_fix
    @test isnan(info.fixed_value)
    @test !info.binary
    @test !info.integer
    @test !info.has_start
    @test isnan(info.start)
    @test names[z] == ""
    @test test_kw == 5
end

function macros_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    @testset "build_constraint on variable" begin
        m = ModelType()
        @variable(m, x)
        @test JuMP.build_constraint(error, x, MOI.GreaterThan(0.0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.GreaterThan{Float64}}
        @test JuMP.build_constraint(error, x, MOI.LessThan(0.0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.LessThan{Float64}}
        @test JuMP.build_constraint(error, x, MOI.EqualTo(0)) isa JuMP.SingleVariableConstraint{VariableRefType, MOI.EqualTo{Int}}
    end

    @testset "Check @constraint basics" begin
        m = ModelType()
        @variable(m, w)
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        t = 10.0

        cref = @constraint(m, 3x - y == 3.3(w + 2z) + 5)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, 3x - y == (w + 2z)*3.3 + 5)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 3*x - y - 3.3*w - 6.6*z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, (x+y)/2 == 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 0.5*x + 0.5*y)
        @test c.set == MOI.EqualTo(1.0)

        cref = @constraint(m, -1 <= x-y <= t)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x - y)
        @test c.set == MOI.Interval(-1.0, t)

        cref = @constraint(m, -1 <= x+1 <= 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-2.0, 0.0)

        cref = @constraint(m, -1 <= x <= 1)
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, -1 <= x <= sum(0.5 for i = 1:2))
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, 1 >= x >= 0)
        c = JuMP.constraint_object(cref)
        @test c.func == x
        @test c.set == MOI.Interval(0.0, 1.0)

        @test_throws ErrorException @constraint(m, x <= t <= y)
        @test_throws ErrorException @constraint(m, 0 <= Dict() <= 1)
        @test_macro_throws ErrorException @constraint(1 <= x <= 2, foo=:bar)

        @test JuMP.isequal_canonical(@expression(m, 3x - y - 3.3(w + 2z) + 5), 3*x - y - 3.3*w - 6.6*z + 5)
        @test JuMP.isequal_canonical(@expression(m, quad, (w+3)*(2x+1)+10), 2*w*x + 6*x + w + 13)

        cref = @constraint(m, 3 + 5*7 <= 0)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, zero(AffExpr))
        @test c.set == MOI.LessThan(-38.0)
    end

end

@testset "Macros for JuMP.Model" begin
    macros_test(Model, VariableRef)

    @testset "Nested tuple destructuring" begin
        m = Model()
        d = Dict((1,2) => 3)
        ex = @expression(m, sum(i+j+k for ((i,j),k) in d))
        @test ex == 6
    end

    @testset "Error on unexpected comparison" begin
        m = Model()
        @variable(m, x)
        @test_macro_throws ErrorException @expression(m, x <= 1)
    end

    @testset "Warn on unexpected assignment" begin
        m = Model()
        @variable(m, x)
        @static if VERSION >= v"1.0"
            # function getindex does not accept keyword arguments
            @test_throws ErrorException x[i=1]
            @test_throws ErrorException @constraint(m, x[i=1] <= 1)
        else
            @static if VERSION >= v"0.7-"
                # https://github.com/JuliaLang/julia/issues/25612
                @test_logs((:warn, r"Unexpected assignment"),
                           macroexpand(JuMP, :(@constraint(m, x[i=1] <= 1))))
            else
                @test_warn "Unexpected assignment" macroexpand(:(@constraint(m, x[i=1] <= 1)))
            end
        end
    end
end

@testset "Macros for JuMPExtension.MyModel" begin
    macros_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
