#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/macros.jl
# Testing for macros
#############################################################################

using JuMP
using Test
using Base.Meta

const MA = JuMP._MA

include(joinpath(@__DIR__, "utilities.jl"))

@static if !(:JuMPExtension in names(Main))
    include(joinpath(@__DIR__, "JuMPExtension.jl"))
end

@testset "Check Julia generator expression parsing" begin
    sumexpr = :(sum(x[i, j] * y[i, j] for i = 1:N, j in 1:M if i != j))
    @test sumexpr.head == :call
    @test sumexpr.args[1] == :sum
    @test sumexpr.args[2].head == :generator
    @test sumexpr.args[2].args[1] == :(x[i, j] * y[i, j])
    @test sumexpr.args[2].args[2].head == :filter
    @test sumexpr.args[2].args[2].args[1] == :(i != j)
    @test sumexpr.args[2].args[2].args[2] == :(i = 1:N)
    @test sumexpr.args[2].args[2].args[3] == :(j = 1:M)

    sumexpr = :(sum(x[i, j] * y[i, j] for i in 1:N, j in 1:M))
    @test sumexpr.head == :call
    @test sumexpr.args[1] == :sum
    @test sumexpr.args[2].head == :generator
    @test sumexpr.args[2].args[1] == :(x[i, j] * y[i, j])
    @test sumexpr.args[2].args[2] == :(i = 1:N)
    @test sumexpr.args[2].args[3] == :(j = 1:M)
end

@testset "Check Julia condition expression parsing" begin
    ex = :(x[12; 3])
    @test ex.head == :typed_vcat
    @test ex.args == [:x, 12, 3]

    ex = :(x[i = 1:3, j = S; isodd(i) && i + j >= 2])
    @test ex.head == :ref
    @test ex.args == [
        :x,
        Expr(:parameters, Expr(:&&, :(isodd(i)), :(i + j >= 2))),
        Expr(:kw, :i, :(1:3)),
        Expr(:kw, :j, :S),
    ]
end

@testset "Test _add_positional_args" begin
    call = :(f(1, a = 2))
    @test JuMP._add_positional_args(call, [:(MyObject)]) isa Nothing
    @test call == :(f(1, $(Expr(:escape, :MyObject)), a = 2))
end

@testset "MutableArithmetics.Zero (Issue #2187)" begin
    model = Model()
    c = @constraint(model, sum(1 for _ in 1:0) == sum(1 for _ in 1:0))
    @test constraint_object(c).func == AffExpr(0.0)
    @test constraint_object(c).set == MOI.EqualTo(0.0)
end

@testset "MutableArithmetics.Zero (Issue #2087)" begin
    model = Model()
    @objective(model, Min, sum(1 for _ in 1:0))
    @test objective_function(model) == AffExpr(0.0)
    c = @constraint(model, sum(1 for _ in 1:0) in MOI.EqualTo(0.0))
    @test constraint_object(c).func == AffExpr(0.0)
    @test constraint_object(c).set == MOI.EqualTo(0.0)
end

struct NewVariable <: JuMP.AbstractVariable
    info::JuMP.VariableInfo
end

@testset "Extension variables constrained on creation #2594" begin
    function JuMP.build_variable(
        _error::Function,
        info::JuMP.VariableInfo,
        ::Type{NewVariable},
    )
        return NewVariable(info)
    end
    function JuMP.add_variable(model::Model, v::NewVariable, name::String = "")
        return JuMP.add_variable(
            model,
            ScalarVariable(v.info),
            name * "_normal_add",
        )
    end
    function JuMP.add_variable(
        model::Model,
        v::VariablesConstrainedOnCreation{
            MOI.SecondOrderCone,
            VectorShape,
            NewVariable,
        },
        names,
    )
        vs = map(i -> ScalarVariable(i.info), v.scalar_variables)
        new_v = VariablesConstrainedOnCreation(vs, v.set, v.shape)
        names .*= "_constr_add"
        return JuMP.add_variable(model, new_v, names)
    end

    model = Model()
    @variable(model, 0 <= x <= 1, NewVariable, Bin)
    @test lower_bound(x) == 0
    @test upper_bound(x) == 1
    @test is_binary(x)
    @test name(x) == "x_normal_add"

    @variable(model, y[1:3] in SecondOrderCone(), NewVariable)
    @test name.(y) == ["y[$i]_constr_add" for i in 1:3]
    @test num_constraints(model, Vector{VariableRef}, MOI.SecondOrderCone) == 1
end

mutable struct MyVariable
    test_kw::Int
    info::JuMP.VariableInfo
end

@testset "Extension of @variable with build_variable #1029" begin
    local MyVariable{S,T,U,V} = Tuple{JuMP.VariableInfo{S,T,U,V},Int,Int}
    names = Dict{MyVariable,String}()
    function JuMP.add_variable(::Model, v::MyVariable, name::String = "")
        names[v] = name
        return v
    end
    # Since `VariableInfo` is an immutable struct, two objects with the same
    # fields have the same hash hence we need to add an id to distinguish
    # variables in the `names` dictionary.
    id = 0
    function JuMP.build_variable(
        _error::Function,
        info::JuMP.VariableInfo,
        ::Type{MyVariable};
        test_kw::Int = 0,
    )
        return (info, test_kw, id += 1)
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
    @test isa(y, Vector{<:MyVariable})
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

    z = @variable(m, variable_type = MyVariable, upper_bound = 3, test_kw = 5)
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

struct PowerCone{T}
    exponent::T
end
function JuMP.build_constraint(
    _error::Function,
    f,
    set::PowerCone;
    dual = false,
)
    moi_set =
        dual ? MOI.DualPowerCone(set.exponent) : MOI.PowerCone(set.exponent)
    return JuMP.build_constraint(_error, f, moi_set)
end
function build_constraint_keyword_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "build_constraint with keyword arguments" begin
        model = ModelType()
        @variable(model, x)
        cref1 = @constraint(model, [1, x, x] in PowerCone(0.5))
        @test JuMP.constraint_object(cref1).set isa MOI.PowerCone{Float64}
        cref2 = @constraint(model, [1, x, x] in PowerCone(0.5), dual = true)
        @test JuMP.constraint_object(cref2).set isa MOI.DualPowerCone{Float64}
    end
end

struct MyConstrType end
struct BadPosArg end
function JuMP.build_constraint(
    _error::Function,
    f::GenericAffExpr,
    set::MOI.EqualTo,
    extra::Type{MyConstrType};
    d = 0,
)
    new_set = MOI.LessThan{Float64}(set.value + d)
    return JuMP.build_constraint(_error, f, new_set)
end
function build_constraint_extra_arg_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "build_constraint with extra positional arguments" begin
        model = ModelType()
        @variable(model, x)
        cref = @constraint(model, x == 0, MyConstrType)
        @test JuMP.constraint_object(cref).set isa MOI.LessThan{Float64}
        cref = @constraint(model, c1, x == 0, MyConstrType, d = 1)
        @test JuMP.constraint_object(cref).set == MOI.LessThan{Float64}(1)
        @test_throws_strip ErrorException @constraint(model, x == 0, BadPosArg)
        @test_throws_strip ErrorException @constraint(
            model,
            x == 0,
            BadPosArg,
            d = 1
        )
        @test_macro_throws ErrorException @constraint(
            model,
            x == 0,
            MyConstrType,
            BadPosArg
        )
    end
end

struct CustomType end
function JuMP.parse_constraint_head(_error::Function, ::Val{:(:=)}, lhs, rhs)
    return false, :(), :(build_constraint($_error, $(esc(lhs)), $(esc(rhs))))
end
struct CustomSet <: MOI.AbstractScalarSet end
Base.copy(set::CustomSet) = set
function JuMP.build_constraint(_error::Function, func, ::CustomType)
    return JuMP.build_constraint(_error, func, CustomSet())
end
function custom_expression_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "Custom expression" begin
        model = ModelType()
        @variable(model, x)
        @constraint(model, con_ref, x := CustomType())
        con = JuMP.constraint_object(con_ref)
        @test jump_function(con) == x
        @test moi_set(con) isa CustomSet
    end
end

function JuMP.parse_one_operator_constraint(
    _error::Function,
    ::Bool,
    ::Val{:f},
    x,
)
    return :(), :(build_constraint($_error, $(esc(x)), $(esc(CustomType()))))
end
function custom_function_test(ModelType::Type{<:JuMP.AbstractModel})
    @testset "Custom function" begin
        model = ModelType()
        @variable(model, x)
        @constraint(model, con_ref, f(x))
        con = JuMP.constraint_object(con_ref)
        @test jump_function(con) == x
        @test moi_set(con) isa CustomSet

        @test_macro_throws ErrorException @constraint(model, g(x))
    end
end

function macros_test(
    ModelType::Type{<:JuMP.AbstractModel},
    VariableRefType::Type{<:JuMP.AbstractVariableRef},
)
    @testset "build_constraint on variable" begin
        m = ModelType()
        @variable(m, x)
        @test JuMP.build_constraint(error, x, MOI.GreaterThan(0.0)) isa
              JuMP.ScalarConstraint{VariableRefType,MOI.GreaterThan{Float64}}
        @test JuMP.build_constraint(error, x, MOI.LessThan(0.0)) isa
              JuMP.ScalarConstraint{VariableRefType,MOI.LessThan{Float64}}
        @test JuMP.build_constraint(error, x, MOI.EqualTo(0)) isa
              JuMP.ScalarConstraint{VariableRefType,MOI.EqualTo{Int}}
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
        @test JuMP.isequal_canonical(c.func, 3 * x - y - 3.3 * w - 6.6 * z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, 3x - y == (w + 2z) * 3.3 + 5)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 3 * x - y - 3.3 * w - 6.6 * z)
        @test c.set == MOI.EqualTo(5.0)

        cref = @constraint(m, (x + y) / 2 == 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 0.5 * x + 0.5 * y)
        @test c.set == MOI.EqualTo(1.0)

        cref = @constraint(m, -1 <= x - y <= t)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x - y)
        @test c.set == MOI.Interval(-1.0, t)

        cref = @constraint(m, -1 <= x + 1 <= 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-2.0, 0.0)

        cref = @constraint(m, -1 <= x <= 1)
        c = JuMP.constraint_object(cref)
        @test c.func isa JuMP.GenericAffExpr
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, -1 <= x <= sum(0.5 for i in 1:2))
        c = JuMP.constraint_object(cref)
        @test c.func isa JuMP.GenericAffExpr
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(-1.0, 1.0)

        cref = @constraint(m, 1 >= x >= 0)
        c = JuMP.constraint_object(cref)
        @test c.func isa JuMP.GenericAffExpr
        @test JuMP.isequal_canonical(c.func, 1x)
        @test c.set == MOI.Interval(0.0, 1.0)

        @test_throws ErrorException @constraint(m, x <= t <= y)
        @test_throws ErrorException @constraint(m, 0 <= Dict() <= 1)
        @test_macro_throws ErrorException @constraint(1 <= x <= 2, foo = :bar)

        @test JuMP.isequal_canonical(
            @expression(m, 3x - y - 3.3(w + 2z) + 5),
            3 * x - y - 3.3 * w - 6.6 * z + 5,
        )
        @test JuMP.isequal_canonical(
            @expression(m, quad, (w + 3) * (2x + 1) + 10),
            2 * w * x + 6 * x + w + 13,
        )

        cref = @constraint(m, 3 + 5 * 7 <= 0)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, zero(AffExpr))
        @test c.set == MOI.LessThan(-38.0)
    end

    @testset "Unicode comparison operators" begin
        model = ModelType()
        @variable(model, x)
        @variable(model, y)
        t = 10.0

        cref = @constraint(model, (x + y) / 2 ≥ 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 0.5 * x + 0.5 * y)
        @test c.set == MOI.GreaterThan(1.0)

        cref = @constraint(model, (x + y) / 2 ≤ 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, 0.5 * x + 0.5 * y)
        @test c.set == MOI.LessThan(1.0)

        cref = @constraint(model, -1 ≤ x - y ≤ t)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x - y)
        @test c.set == MOI.Interval(-1.0, t)

        cref = @constraint(model, 1 ≥ x ≥ 0)
        c = JuMP.constraint_object(cref)
        @test c.func isa JuMP.GenericAffExpr
        @test JuMP.isequal_canonical(c.func, 1 * x)
        @test c.set == MOI.Interval(0.0, 1.0)
    end

    @testset "Constraint Naming" begin
        model = ModelType()
        @variable(model, x)

        cref = @constraint(model, x == 0)
        @test name(cref) == ""

        cref = @constraint(model, x == 0, base_name = "cat")
        @test name(cref) == "cat"

        cref = @constraint(model, c1, x == 0)
        @test name(cref) == "c1"

        cref = @constraint(model, c2, x == 0, base_name = "cat")
        @test name(cref) == "cat"

        crefs = @constraint(model, [1:2], x == 0, base_name = "cat")
        @test name.(crefs) == ["cat[1]", "cat[2]"]

        @test_macro_throws ErrorException @constraint(model, c3[1:2])
        @test_macro_throws ErrorException @constraint(model, "c"[1:2], x == 0)
    end

    @testset "@build_constraint (scalar inequality)" begin
        model = ModelType()
        @variable(model, x)
        con = @build_constraint(3x == 1)
        @test con isa JuMP.ScalarConstraint
        @test JuMP.isequal_canonical(con.func, 3x)
        @test con.set == MOI.EqualTo(1.0)
    end

    @testset "@build_constraint (function-in-set)" begin
        model = ModelType()
        @variable(model, x[1:2])
        con = @build_constraint(x in JuMP.SecondOrderCone())
        @test con isa JuMP.VectorConstraint
        @test con.func == x
        @test con.set == MOI.SecondOrderCone(2)
    end

    @testset "@build_constraint (SOS1)" begin
        model = ModelType()
        @variable(model, x[1:3])
        con = @build_constraint(x in JuMP.SOS1())
        con2 = @build_constraint(x in JuMP.SOS1([4.0, 6.0, 1.0]))
        @test con isa JuMP.VectorConstraint
        @test con.func == x
        @test con.set == MOI.SOS1([1.0, 2.0, 3.0])
        @test_throws(
            ErrorException("Weight vector in SOS1 is not of length 3."),
            @build_constraint(x in JuMP.SOS1([1.0]))
        )
        @test con2 isa JuMP.VectorConstraint
        @test con2.func == x
        @test con2.set == MOI.SOS1([4.0, 6.0, 1.0])
    end

    @testset "@build_constraint (SOS2)" begin
        model = ModelType()
        @variable(model, x[1:3])
        con = @build_constraint(x in JuMP.SOS2())
        con2 = @build_constraint(x in JuMP.SOS2([4.0, 6.0, 1.0]))
        @test con isa JuMP.VectorConstraint
        @test con.func == x
        @test con.set == MOI.SOS2([1.0, 2.0, 3.0])
        @test_throws(
            ErrorException("Weight vector in SOS2 is not of length 3."),
            @build_constraint(x in JuMP.SOS2([1.0]))
        )
        @test con2 isa JuMP.VectorConstraint
        @test con2.func == x
        @test con2.set == MOI.SOS2([4.0, 6.0, 1.0])
    end

    @testset "@build_constraint (broadcast)" begin
        model = ModelType()
        @variable(model, x[1:2])
        ub = [1.0, 2.0]
        con = @build_constraint(x .<= ub)
        @test con isa Vector{<:JuMP.ScalarConstraint}
        @test JuMP.isequal_canonical(con[1].func, 1.0x[1])
        @test JuMP.isequal_canonical(con[2].func, 1.0x[2])
        @test con[1].set == MOI.LessThan(1.0)
        @test con[2].set == MOI.LessThan(2.0)
    end

    @testset "@build_constraint error" begin
        model = ModelType()
        @variable(model, x)
        @test_macro_throws ErrorException @build_constraint(2x + 1)
    end

    @testset "Promotion of SOS sets" begin
        model = Model()
        @variable(model, x[1:3])
        c_sos1 = @build_constraint(x in MOI.SOS1([1, 2, 3]))
        @test c_sos1.set == MOI.SOS1([1.0, 2.0, 3.0])
        c_sos2 = @build_constraint(x in MOI.SOS2([6, 5, 4]))
        @test c_sos2.set == MOI.SOS2([6.0, 5.0, 4.0])
    end

    build_constraint_keyword_test(ModelType)
    custom_expression_test(ModelType)
    build_constraint_extra_arg_test(ModelType)
    return custom_function_test(ModelType)
end

@testset "Macros for JuMP.Model" begin
    macros_test(Model, VariableRef)

    @testset "Adding anonymous variable and specify required constraint on it" begin
        model = Model()
        @test_macro_throws(
            ErrorException(
                "In `@variable(m, Int)`: Ambiguous variable name Int detected." *
                " To specify an anonymous integer variable, use `@variable(model, integer = true)` instead.",
            ),
            @variable(m, Int)
        )
        v = @variable(model, integer = true)
        @test name(v) == ""
        @test is_integer(v)

        @test_macro_throws(
            ErrorException(
                "In `@variable(m, Bin)`: Ambiguous variable name Bin detected." *
                " To specify an anonymous binary variable, use `@variable(model, binary = true)` instead.",
            ),
            @variable(m, Bin)
        )
        v = @variable(model, binary = true)
        @test name(v) == ""
        @test is_binary(v)

        @test_macro_throws(
            ErrorException(
                "In `@variable(m, PSD)`: Size of anonymous square matrix of positive semidefinite anonymous variables is not specified." *
                " To specify size of square matrix use `@variable(model, [1:n, 1:n], PSD)` instead.",
            ),
            @variable(m, PSD)
        )
        v = @variable(model, [1:1, 1:1], PSD)
        @test name(v[1]) == ""
    end

    @testset "Nested tuple destructuring" begin
        m = Model()
        d = Dict((1, 2) => 3)
        ex = @expression(m, sum(i + j + k for ((i, j), k) in d))
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
        # Julia v1.0 -> v1.3
        # ERROR: function getindex does not accept keyword arguments
        # Julia v1.3 onwards
        # ERROR: MethodError: no method matching getindex(::VariableRef; i=1)
        @test_throws Union{ErrorException,MethodError} x[i = 1]
        err = ErrorException("Unexpected assignment in expression `x[i=1]`.")
        @test_macro_throws ErrorException @constraint(m, x[i = 1] <= 1)
    end

    @testset "Lookup in model scope: @variable" begin
        model = Model()
        @variable(model, x)
        @test x === model[:x]
    end

    @testset "Lookup in model scope: @constraint" begin
        model = Model()
        @variable(model, x)
        @constraint(model, con, x + 1 <= 2)
        @test con === model[:con]
    end

    @testset "Lookup in model scope: @expression" begin
        model = Model()
        @variable(model, x)
        @expression(model, expr, 2x)
        @test expr === model[:expr]
    end

    @testset "Lookup in model scope: @NLexpression" begin
        model = Model()
        @variable(model, x)
        @NLexpression(model, nl_expr, sin(x))
        @test nl_expr === model[:nl_expr]
    end

    @testset "Lookup in model scope: @NLconstraint" begin
        model = Model()
        @variable(model, x)
        @NLconstraint(model, nl_con, sin(x) == 1.0)
        @test nl_con === model[:nl_con]
    end

    @testset "Error on duplicate names in model scope" begin
        model = Model()
        y = @variable(model, x)
        @test_throws ErrorException @constraint(model, x, 2y <= 1)
    end

    @testset "Duplicate variables are never created" begin
        model = Model()
        @variable(model, x)
        @test_throws ErrorException @variable(model, x)
        @test 1 == @inferred num_variables(model)
    end

    @testset "Anonymous expressions aren't registered" begin
        model = Model()
        x = @variable(model)
        ex = @expression(model, x + 1)
        @test length(JuMP.object_dictionary(model)) == 0
    end

    @testset "Anonymous NLexpressions aren't registered" begin
        model = Model()
        x = @variable(model)
        ex = @NLexpression(model, x + 1)
        @test length(JuMP.object_dictionary(model)) == 0
    end

    @testset "Adjoints" begin
        model = Model()
        @variable(model, x[1:2])
        obj = @objective(model, Min, x' * ones(2, 2) * x)
        @test JuMP.isequal_canonical(obj, x[1]^2 + 2 * x[1] * x[2] + x[2]^2)
        cref = @constraint(model, x' * ones(2, 2) * x <= 1)
        c = JuMP.constraint_object(cref)
        @test JuMP.isequal_canonical(c.func, x[1]^2 + 2 * x[1] * x[2] + x[2]^2)
        @test c.set == MOI.LessThan(1.0)
        @test JuMP.isequal_canonical(
            MA.@rewrite(x' * ones(2, 2)),
            x' * ones(2, 2),
        )
    end

    @testset "Nonliteral exponents in @constraint" begin
        model = Model()
        @variable(model, x)
        foo() = 2
        con1 = @build_constraint(x^(foo()) + x^(foo() - 1) + x^(foo() - 2) == 0)
        con2 = @build_constraint(
            (x - 1)^(foo()) + (x - 1)^2 + (x - 1)^1 + (x - 1)^0 == 0
        )
        con3 = @build_constraint(sum(x for i in 1:3)^(foo()) == 0)
        con4 = @build_constraint(sum(x for i in 1:3)^(foo() - 1) == 0)
        @test con1.func == x^2 + x
        @test con2.func == 2 * x^2 - 3 * x
        @test con3.func == 9 * x^2
        @test con4.func == convert(QuadExpr, 3 * x)
    end

    @testset "AffExpr in macros" begin
        eq = JuMP._math_symbol(REPLMode, :eq)
        ge = JuMP._math_symbol(REPLMode, :geq)

        model = Model()
        @variable(model, x)
        @variable(model, y)
        temp = x + 2 * y + 1
        a = 1.0 * x
        con1 = @constraint(model, 3 * temp - x - 2 >= 0)
        con2 = @constraint(model, (2 + 2) * ((3 + 4) * (1 + a)) == 0)
        con3 = @constraint(model, 1 + 0 * temp == 0)
        @test string(con1) == "2 x + 6 y $ge -1.0"
        @test string(con2) == "28 x $eq -28.0"
        @test string(con3) == "0 $eq -1.0"
    end

    @testset "@constraints" begin
        eq = JuMP._math_symbol(REPLMode, :eq)
        ge = JuMP._math_symbol(REPLMode, :geq)

        model = Model()
        @variable(model, x)
        @variable(model, y[1:3])

        @constraints(model, begin
            x + y[1] == 1
            ref[i = 1:3], y[1] + y[i] >= i
        end)

        @test string(model) == """
        Feasibility
        Subject to
         x + y[1] $eq 1.0
         ref[1] : 2 y[1] $ge 1.0
         ref[2] : y[1] + y[2] $ge 2.0
         ref[3] : y[1] + y[3] $ge 3.0
        """
    end

    @testset "NLparameters" begin
        model = Model()
        @NLparameters(model, begin
            a == 1
            b[i = 1:2] == i
        end)
        @test value(a) == 1
        @test value.(b) == [1, 2]
    end

    @testset "Index variables don't leak out of macros" begin
        model = Model()
        i = 10
        j = 10
        @expression(model, ex[j = 2:3], sum(i for i in 1:j))
        @test ex[2] == 3
        @test ex[3] == 6
        @test i == 10
        @test j == 10
    end

    @testset "Plural failures" begin
        model = Model()
        @test_macro_throws MethodError @variables(model)
        @test_macro_throws ErrorException("Invalid syntax for @variables") @variables(
            model,
            x
        )
        @test_macro_throws ErrorException("Invalid syntax for @variables") @variables(
            model,
            x >= 0
        )
        @test_macro_throws MethodError @variables(model, x >= 0, Bin)
    end

    @testset "Empty summation in @constraints" begin
        model = Model()
        @variable(model, x)
        c = @constraint(model, x == sum(1.0 for i in 1:0))
        @test isa(
            constraint_object(c).func,
            GenericAffExpr{Float64,VariableRef},
        )
    end

    @testset "Empty summation in @NLconstraints" begin
        model = Model()
        @variable(model, x)
        c = @NLconstraint(model, x == sum(1.0 for i in 1:0))
        @test sprint(show, c) == "x - 0 = 0" || sprint(show, c) == "x - 0 == 0"
    end

    @testset "Splatting error" begin
        model = Model()
        A = [1 0; 0 1]
        @variable(model, x)

        @test_macro_throws(
            ErrorException(
                "In `@variable(model, y[axes(A)...])`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @variable(model, y[axes(A)...])
        )

        f(a, b) = [a, b]
        @variable(model, z[f((1, 2)...)])
        @test length(z) == 2

        @test_macro_throws(
            ErrorException(
                "In `@constraint(model, [axes(A)...], x >= 1)`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @constraint(model, [axes(A)...], x >= 1)
        )

        @test_macro_throws(
            ErrorException(
                "In `@NLconstraint(model, [axes(A)...], x >= 1)`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @NLconstraint(model, [axes(A)...], x >= 1)
        )

        @test_macro_throws(
            ErrorException(
                "In `@expression(model, [axes(A)...], x)`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @expression(model, [axes(A)...], x)
        )

        @test_macro_throws(
            ErrorException(
                "In `@NLexpression(model, [axes(A)...], x)`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @NLexpression(model, [axes(A)...], x)
        )

        @test_macro_throws(
            ErrorException(
                "In `@NLparameter(model, p[axes(A)...] == x)`: cannot use splatting operator `...` in the definition of an index set.",
            ),
            @NLparameter(model, p[axes(A)...] == x)
        )
    end

    @testset "NaN in constraints" begin
        model = Model()
        @variable(model, x >= 0)
        @test_throws(
            ErrorException(
                "Expression contains an invalid NaN constant. This could be produced by `Inf - Inf`.",
            ),
            @constraint(model, x >= NaN)
        )
        @test_throws ErrorException(
            "Expression contains an invalid NaN constant. This could be produced by `Inf - Inf`.",
        ) @constraint(model, 1 <= x + NaN <= 2)
        @test_throws ErrorException(
            "Expression contains an invalid NaN constant. This could be produced by `Inf - Inf`.",
        ) @constraint(model, 1 <= x + Inf <= 2)
        @test_throws_strip(
            ErrorException(
                "In `@constraint(model, 1 <= x <= NaN)`: Invalid bounds, cannot contain NaN: [1, NaN].",
            ),
            @constraint(model, 1 <= x <= NaN)
        )
    end

    @testset "unregister" begin
        model = Model()
        @variable(model, x)
        @test num_variables(model) == 1
        @test_throws ErrorException @variable(model, x)
        unregister(model, :x)
        @variable(model, x)
        @test num_variables(model) == 2
    end

    @testset "unrecognized_variable_type" begin
        model = Model()
        err = ErrorException(
            "In `@variable(model, x, 2, variable_type = 1)`: Unrecognized " *
            "arguments: 2, 1. (You may have passed these as positional " *
            "arguments, or as a keyword value to `variable_type`.)\n\nIf " *
            "you're trying to create a JuMP extension, you need to implement " *
            "`build_variable`.",
        )
        @test_throws_strip err @variable(model, x, 2, variable_type = 1)
    end

    @testset "kwargs" begin
        model = Model()
        x = @variable(model; integer = true)
        @test is_integer(x)
        x = @variable(model; integer = false)
        @test !is_integer(x)
        @variable(model, y; integer = true)
        @test is_integer(y)
    end

    @testset "Model as index" begin
        m = Model()
        @variable(m, x)
        index_set = VERSION < v"1.1" ? "m=1:2" : "m = 1:2"
        msg = "Index m is the same symbol as the model. Use a different name for the index."
        @test_macro_throws(
            ErrorException("In `@variable(m, y[$(index_set)] <= m)`: $(msg)"),
            @variable(m, y[m = 1:2] <= m),
        )
        @test_macro_throws(
            ErrorException("In `@constraint(m, [m = 1:2], x <= m)`: $(msg)"),
            @constraint(m, [m = 1:2], x <= m),
        )
        @test_macro_throws(
            ErrorException("In `@expression(m, [m = 1:2], m * x)`: $(msg)"),
            @expression(m, [m = 1:2], m * x),
        )
        @test_macro_throws(
            ErrorException(
                "In `@NLconstraint(m, [m = 1:2], sqrt(x) <= m)`: $(msg)",
            ),
            @NLconstraint(m, [m = 1:2], sqrt(x) <= m),
        )
        @test_macro_throws(
            ErrorException("In `@NLexpression(m, [m = 1:2], x)`: $(msg)"),
            @NLexpression(m, [m = 1:2], x),
        )
        @test_macro_throws(
            ErrorException(
                "In `@NLparameter(m, p[$(index_set)] == m)`: $(msg)",
            ),
            @NLparameter(m, p[m = 1:2] == m),
        )
    end
end

@testset "Macros for JuMPExtension.MyModel" begin
    macros_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
