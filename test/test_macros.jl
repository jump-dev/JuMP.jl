#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestMacros

using JuMP
using Test

import LinearAlgebra
import SparseArrays

include(joinpath(@__DIR__, "utilities.jl"))

struct NewVariable <: AbstractVariable
    info::VariableInfo
end

function JuMP.build_variable(
    ::Function,
    info::VariableInfo,
    ::Type{NewVariable},
)
    return NewVariable(info)
end

function JuMP.add_variable(
    model::GenericModel,
    v::NewVariable,
    name::String = "",
)
    return add_variable(model, ScalarVariable(v.info), name * "_normal_add")
end

function JuMP.add_variable(
    model::GenericModel,
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
    return add_variable(model, new_v, names)
end

mutable struct MyVariable
    test_kw::Int
    info::VariableInfo
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
    return build_constraint(_error, f, moi_set)
end

struct MyConstrType end

struct BadPosArg end

function JuMP.build_constraint(
    _error::Function,
    f::GenericAffExpr,
    set::MOI.EqualTo,
    ::Type{MyConstrType};
    d = 0,
)
    new_set = MOI.LessThan{Float64}(set.value + d)
    return build_constraint(_error, f, new_set)
end

struct CustomType end

function JuMP.parse_constraint_head(_error::Function, ::Val{:(:=)}, lhs, rhs)
    return false, :(), :(build_constraint($_error, $(esc(lhs)), $(esc(rhs))))
end

struct CustomSet <: MOI.AbstractScalarSet end

Base.copy(set::CustomSet) = set

function JuMP.build_constraint(_error::Function, func, ::CustomType)
    return build_constraint(_error, func, CustomSet())
end

function JuMP.parse_constraint_call(_error::Function, ::Bool, ::Val{:f}, x)
    return :(), :(build_constraint($_error, $(esc(x)), $(esc(CustomType()))))
end

const MyVariableTuple{S,T,U,V} = Tuple{VariableInfo{S,T,U,V},Int,Int}

function JuMP.add_variable(
    model::GenericModel,
    v::MyVariableTuple,
    name::String = "",
)
    model.ext[:names][v] = name
    return v
end

# Since `VariableInfo` is an immutable struct, two objects with the same
# fields have the same hash hence we need to add an id to distinguish
# variables in the `names` dictionary.
let id = 0
    function JuMP.build_variable(
        ::Function,
        info::VariableInfo,
        ::Type{MyVariableTuple};
        test_kw::Int = 0,
    )
        return (info, test_kw, id += 1)
    end
end

function test_Check_Julia_generator_expression_parsing()
    sumexpr = :(sum(x[i, j] * y[i, j] for i in 1:N, j in 1:M if i != j))
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
    return
end

function test_Check_Julia_condition_expression_parsing()
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
    return
end

function test_add_positional_args()
    call = :(f(1; a = 2))
    @test JuMP._add_positional_args(call, [:(MyObject)]) isa Nothing
    @test call == :(f(1, $(Expr(:escape, :MyObject)); a = 2))
    return
end

function test_MutableArithmetics_Zero_Issue_2187()
    model = Model()
    c = @constraint(model, sum(1 for _ in 1:0) == sum(1 for _ in 1:0))
    @test constraint_object(c).func == AffExpr(0.0)
    @test constraint_object(c).set == MOI.EqualTo(0.0)
    return
end

function test_MutableArithmetics_Zero_Issue_2087()
    model = Model()
    @objective(model, Min, sum(1 for _ in 1:0))
    @test objective_function(model) == AffExpr(0.0)
    c = @constraint(model, sum(1 for _ in 1:0) in MOI.EqualTo(0.0))
    @test constraint_object(c).func == AffExpr(0.0)
    @test constraint_object(c).set == MOI.EqualTo(0.0)
    return
end

function test_variables_constrained_on_creation_2594()
    model = Model()
    @variable(model, 0 <= x <= 1, NewVariable, Bin)
    @test lower_bound(x) == 0
    @test upper_bound(x) == 1
    @test is_binary(x)
    @test name(x) == "x_normal_add"
    @variable(model, y[1:3] in SecondOrderCone(), NewVariable)
    @test name.(y) == ["y[$i]_constr_add" for i in 1:3]
    @test num_constraints(model, Vector{VariableRef}, MOI.SecondOrderCone) == 1
    return
end

function test_of_variable_with_build_variable_1029()
    m = Model()
    names = Dict{MyVariableTuple,String}()
    m.ext[:names] = names
    @variable(
        m,
        1 <= x <= 2,
        MyVariableTuple,
        binary = true,
        test_kw = 1,
        start = 3
    )
    @test isa(x, MyVariableTuple)
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
    @variable(m, y[1:3] >= 0, MyVariableTuple, test_kw = 2)
    @test isa(y, Vector{<:MyVariableTuple})
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
    z = @variable(
        m,
        variable_type = MyVariableTuple,
        upper_bound = 3,
        test_kw = 5
    )
    info = z[1]
    test_kw = z[2]
    @test isa(z, MyVariableTuple)
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
    return
end

function test_extension_build_constraint_keyword_test(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    cref1 = @constraint(model, [1, x, x] in PowerCone(0.5))
    @test constraint_object(cref1).set isa MOI.PowerCone{Float64}
    cref2 = @constraint(model, [1, x, x] in PowerCone(0.5), dual = true)
    @test constraint_object(cref2).set isa MOI.DualPowerCone{Float64}
    return
end

function test_extension_build_constraint_extra_arg_test(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    cref = @constraint(model, x == 0, MyConstrType)
    @test constraint_object(cref).set isa MOI.LessThan{Float64}
    cref = @constraint(model, c1, x == 0, MyConstrType, d = 1)
    @test constraint_object(cref).set == MOI.LessThan{Float64}(1)
    @test_throws_strip ErrorException @constraint(model, x == 0, BadPosArg)
    @test_throws_strip(
        ErrorException,
        @constraint(model, x == 0, BadPosArg, d = 1),
    )
    @test_macro_throws(
        ErrorException,
        @constraint(model, x == 0, MyConstrType, BadPosArg),
    )
    return
end

function test_extension_custom_expression_test(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, con_ref, x := CustomType())
    con = constraint_object(con_ref)
    @test jump_function(con) == x
    @test moi_set(con) isa CustomSet
    return
end

function test_extension_custom_function_test(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, con_ref, f(x))
    con = constraint_object(con_ref)
    @test jump_function(con) == x
    @test moi_set(con) isa CustomSet
    @test_macro_throws ErrorException @constraint(model, g(x))
    return
end

function test_extension_build_constraint_on_variable(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test build_constraint(error, x, MOI.GreaterThan(0.0)) isa
          ScalarConstraint{VariableRefType,MOI.GreaterThan{Float64}}
    @test build_constraint(error, x, MOI.LessThan(0.0)) isa
          ScalarConstraint{VariableRefType,MOI.LessThan{Float64}}
    @test build_constraint(error, x, MOI.EqualTo(0)) isa
          ScalarConstraint{VariableRefType,MOI.EqualTo{Int}}
    return
end

function test_extension_check_constraint_basics(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    AffExprType = GenericAffExpr{T,VariableRefType}
    m = ModelType()
    @variable(m, w)
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    t = T(10)
    α = T(33) / T(10)
    cref = @constraint(m, 3x - y == α * (w + 2z) + 5)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 3 * x - y - α * w - 2α * z)
    @test c.set == MOI.EqualTo(T(5))
    cref = @constraint(m, 3x - y == (w + 2z) * α + 5)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 3 * x - y - α * w - 2α * z)
    @test c.set == MOI.EqualTo(T(5))
    cref = @constraint(m, (x + y) / 2 == 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, T(inv(2)) * x + T(inv(2)) * y)
    @test c.set == MOI.EqualTo(T(1))
    cref = @constraint(m, -1 <= x - y <= t)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x - y)
    @test c.set == MOI.Interval(-T(1), t)
    cref = @constraint(m, -1 <= x + 1 <= 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, 1x)
    @test c.set == MOI.Interval(-T(2), T(0))
    cref = @constraint(m, -1 <= x <= 1)
    c = constraint_object(cref)
    @test c.func isa AffExprType
    @test isequal_canonical(c.func, 1x)
    @test c.set == MOI.Interval(-T(1), T(1))
    cref = @constraint(m, -1 <= x <= sum(T(inv(2)) for i in 1:2))
    c = constraint_object(cref)
    @test c.func isa AffExprType
    @test isequal_canonical(c.func, 1x)
    @test c.set == MOI.Interval(-T(1), T(1))
    cref = @constraint(m, 1 >= x >= 0)
    c = constraint_object(cref)
    @test c.func isa AffExprType
    @test isequal_canonical(c.func, 1x)
    @test c.set == MOI.Interval(T(0), T(1))
    @test_throws ErrorException @constraint(m, x <= t <= y)
    @test_throws ErrorException @constraint(m, 0 <= Dict() <= 1)
    @test_macro_throws ErrorException @constraint(1 <= x <= 2, foo = :bar)
    @test isequal_canonical(
        @expression(m, 3x - y - α * (w + 2z) + 5),
        3 * x - y - α * w - T(66) / T(10) * z + 5,
    )
    @test isequal_canonical(
        @expression(m, quad, (w + 3) * (2x + 1) + 10),
        2 * w * x + 6 * x + w + 13,
    )
    cref = @constraint(m, 3 + 5 * 7 <= 0)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, zero(AffExpr))
    @test c.set == MOI.LessThan(-T(38))
    return
end

function test_extension_Unicode_comparison_operators(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    @variable(model, y)
    t = T(10)
    cref = @constraint(model, (x + y) / 2 ≥ 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, inv(T(2)) * x + inv(T(2)) * y)
    @test c.set == MOI.GreaterThan(T(1))
    cref = @constraint(model, (x + y) / 2 ≤ 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, inv(T(2)) * x + inv(T(2)) * y)
    @test c.set == MOI.LessThan(T(1))
    cref = @constraint(model, -1 ≤ x - y ≤ t)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x - y)
    @test c.set == MOI.Interval(-T(1), t)
    cref = @constraint(model, 1 ≥ x ≥ 0)
    c = constraint_object(cref)
    @test c.func isa GenericAffExpr
    @test isequal_canonical(c.func, 1 * x)
    @test c.set == MOI.Interval(T(0), T(1))
    return
end

function test_extension_constraint_naming(
    ModelType = Model,
    VariableRefType = VariableRef,
)
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
    return
end

function test_extension_build_constraint_scalar_inequality(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    con = @build_constraint(3x == 1)
    @test con isa ScalarConstraint
    @test isequal_canonical(con.func, 3x)
    @test con.set == MOI.EqualTo(T(1))
    return
end

function test_extension_build_constraint_function_in_set(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    con = @build_constraint(x in SecondOrderCone())
    @test con isa VectorConstraint
    @test con.func == x
    @test con.set == MOI.SecondOrderCone(2)
    return
end

function test_extension_build_constraint_SOS1(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:3])
    con = @build_constraint(x in SOS1())
    con2 = @build_constraint(x in SOS1(T[4, 6, 1]))
    @test con isa VectorConstraint
    @test con.func == x
    @test con.set == MOI.SOS1(T[1, 2, 3])
    @test_throws(
        ErrorException("Weight vector in SOS1 is not of length 3."),
        @build_constraint(x in SOS1(T[1]))
    )
    @test con2 isa VectorConstraint
    @test con2.func == x
    @test con2.set == MOI.SOS1(T[4, 6, 1])
    return
end

function test_extension_build_constraint_SOS2(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:3])
    con = @build_constraint(x in SOS2())
    con2 = @build_constraint(x in SOS2(T[4, 6, 1]))
    @test con isa VectorConstraint
    @test con.func == x
    @test con.set == MOI.SOS2(T[1, 2, 3])
    @test_throws(
        ErrorException("Weight vector in SOS2 is not of length 3."),
        @build_constraint(x in SOS2(T[1]))
    )
    @test con2 isa VectorConstraint
    @test con2.func == x
    @test con2.set == MOI.SOS2(T[4, 6, 1])
    return
end

function test_extension_build_constraint_broadcast(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x[1:2])
    ub = T[1, 2]
    con = @build_constraint(x .<= ub)
    @test con isa Vector{<:ScalarConstraint}
    @test isequal_canonical(con[1].func, 1.0x[1])
    @test isequal_canonical(con[2].func, 1.0x[2])
    @test con[1].set == MOI.LessThan(T(1))
    @test con[2].set == MOI.LessThan(T(2))
    return
end

function test_extension_extension_build_constraint_error(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test_macro_throws ErrorException @build_constraint(2x + 1)
    return
end

function test_Promotion_of_SOS_sets()
    model = Model()
    @variable(model, x[1:3])
    c_sos1 = @build_constraint(x in MOI.SOS1([1, 2, 3]))
    @test c_sos1.set == MOI.SOS1([1.0, 2.0, 3.0])
    c_sos2 = @build_constraint(x in MOI.SOS2([6, 5, 4]))
    @test c_sos2.set == MOI.SOS2([6.0, 5.0, 4.0])
    return
end

function test_Adding_anonymous_variable_and_specify_required_constraint_on_it()
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
    return
end

function test_Nested_tuple_destructuring()
    m = Model()
    d = Dict((1, 2) => 3)
    ex = @expression(m, sum(i + j + k for ((i, j), k) in d))
    @test ex == 6
    return
end

function test_Error_on_unexpected_comparison()
    m = Model()
    @variable(m, x)
    @test_macro_throws ErrorException @expression(m, x <= 1)
    return
end

function test_Lookup_in_model_scope_variable()
    model = Model()
    @variable(model, x)
    @test x === model[:x]
    return
end

function test_Lookup_in_model_scope_constraint()
    model = Model()
    @variable(model, x)
    @constraint(model, con, x + 1 <= 2)
    @test con === model[:con]
    return
end

function test_Lookup_in_model_scope_expression()
    model = Model()
    @variable(model, x)
    @expression(model, expr, 2x)
    @test expr === model[:expr]
    return
end

function test_Lookup_in_model_scope_NLexpression()
    model = Model()
    @variable(model, x)
    @NLexpression(model, nl_expr, sin(x))
    @test nl_expr === model[:nl_expr]
    return
end

function test_Lookup_in_model_scope_NLconstraint()
    model = Model()
    @variable(model, x)
    @NLconstraint(model, nl_con, sin(x) == 1.0)
    @test nl_con === model[:nl_con]
    return
end

function test_Error_on_duplicate_names_in_model_scope()
    model = Model()
    y = @variable(model, x)
    @test_throws ErrorException @constraint(model, x, 2y <= 1)
    return
end

function test_Duplicate_variables_are_never_created()
    model = Model()
    @variable(model, x)
    @test_throws ErrorException @variable(model, x)
    @test 1 == @inferred num_variables(model)
    return
end

function test_Anonymous_expressions_arent_registered()
    model = Model()
    x = @variable(model)
    @expression(model, x + 1)
    @test length(object_dictionary(model)) == 0
    return
end

function test_Anonymous_NLexpressions_arent_registered()
    model = Model()
    x = @variable(model)
    @NLexpression(model, x + 1)
    @test length(object_dictionary(model)) == 0
    return
end

function test_Adjoints()
    model = Model()
    @variable(model, x[1:2])
    obj = @objective(model, Min, x' * ones(2, 2) * x)
    @test isequal_canonical(obj, x[1]^2 + 2 * x[1] * x[2] + x[2]^2)
    cref = @constraint(model, x' * ones(2, 2) * x <= 1)
    c = constraint_object(cref)
    @test isequal_canonical(c.func, x[1]^2 + 2 * x[1] * x[2] + x[2]^2)
    @test c.set == MOI.LessThan(1.0)
    @test isequal_canonical(JuMP._MA.@rewrite(x' * ones(2, 2)), x' * ones(2, 2))
    return
end

function test_Nonliteral_exponents_in_constraint()
    model = Model()
    @variable(model, x)
    foo() = 2
    con1 = @build_constraint(x^(foo()) + x^(foo() - 1) + x^(foo() - 2) == 0)
    # TODO(odow): `con2` fails to build due to a bug in MutableArithmetics. To
    # fix, we need MutableArithmetics in the current scope.
    MutableArithmetics = JuMP._MA
    con2 = @build_constraint(
        (x - 1)^(foo()) + (x - 1)^2 + (x - 1)^1 + (x - 1)^0 == 0
    )
    con3 = @build_constraint(sum(x for i in 1:3)^(foo()) == 0)
    con4 = @build_constraint(sum(x for i in 1:3)^(foo() - 1) == 0)
    @test con1.func == x^2 + x
    @test con2.func == 2 * x^2 - 3 * x
    @test con3.func == 9 * x^2
    @test con4.func == convert(QuadExpr, 3 * x)
    return
end

function test_AffExpr_in_macros()
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)

    model = Model()
    @variable(model, x)
    @variable(model, y)
    temp = x + 2 * y + 1
    a = 1.0 * x
    con1 = @constraint(model, 3 * temp - x - 2 >= 0)
    con2 = @constraint(model, (2 + 2) * ((3 + 4) * (1 + a)) == 0)
    con3 = @constraint(model, 1 + 0 * temp == 0)
    @test string(con1) == "2 x + 6 y $ge -1"
    @test string(con2) == "28 x $eq -28"
    @test string(con3) == "0 $eq -1"
    return
end

function test_constraints()
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    model = Model()
    @variable(model, x)
    @variable(model, y[1:3])
    @constraints(model, begin
        x + y[1] == 1
        ref[i = 1:3], y[1] + y[i] >= i
    end)
    @test string(model) ==
          "Feasibility\n" *
          "Subject to\n" *
          " x + y[1] $eq 1\n" *
          " ref[1] : 2 y[1] $ge 1\n" *
          " ref[2] : y[1] + y[2] $ge 2\n" *
          " ref[3] : y[1] + y[3] $ge 3\n"
    return
end

function test_NLparameters()
    model = Model()
    @NLparameters(model, begin
        a == 1
        b[i = 1:2] == i
    end)
    @test value(a) == 1
    @test value.(b) == [1, 2]
    return
end

function test_Index_variables_dont_leak_out_of_macros()
    model = Model()
    i = 10
    j = 10
    @expression(model, ex[j = 2:3], sum(i for i in 1:j))
    @test ex[2] == 3
    @test ex[3] == 6
    @test i == 10
    @test j == 10
    return
end

function test_Plural_failures()
    model = Model()
    @test_macro_throws MethodError @variables(model)
    err = ErrorException(
        "Invalid syntax for @variables. The second argument must be a `begin end` " *
        "block. For example:\n" *
        "```julia\n@variables(model, begin\n    # ... lines here ...\nend)\n```.",
    )
    @test_macro_throws err @variables(model, x)
    @test_macro_throws err @variables(model, x >= 0)
    @test_macro_throws MethodError @variables(model, x >= 0, Bin)
    return
end

function test_Plural_returns()
    model = Model()
    vars = @variables(model, begin
        x
        y[1:2]
    end)
    @test vars == (x, y)
    eqs = @constraints(model, begin
        E_x, x == 0
        E_y[i = 1:2], y[i] == 0
    end)
    @test eqs == (E_x, E_y)
    return
end

function test_Empty_summation_in_constraints()
    model = Model()
    @variable(model, x)
    c = @constraint(model, x == sum(1.0 for i in 1:0))
    @test isa(constraint_object(c).func, GenericAffExpr{Float64,VariableRef})
    return
end

function test_Empty_summation_in_NLconstraints()
    model = Model()
    @variable(model, x)
    c = @NLconstraint(model, x == sum(1.0 for i in 1:0))
    @test sprint(show, c) == "x - 0.0 = 0" || sprint(show, c) == "x - 0.0 == 0"
    return
end

function test_Splatting_error()
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
    return
end

function test_NaN_in_constraints()
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
    return
end

function test_unregister()
    model = Model()
    @variable(model, x)
    @test num_variables(model) == 1
    @test_throws ErrorException @variable(model, x)
    unregister(model, :x)
    @variable(model, x)
    @test num_variables(model) == 2
    return
end

function test_unrecognized_variable_type()
    model = Model()
    err = ErrorException(
        "In `@variable(model, x, 2, variable_type = 1)`: " *
        "Unrecognized positional arguments: (2, 1). (You may have " *
        "passed it as a positional argument, or as a keyword value to " *
        "`variable_type`.)\n\nIf you're trying to create a JuMP " *
        "extension, you need to implement `build_variable`. Read the " *
        "docstring for more details.",
    )
    @test_throws_strip err @variable(model, x, 2, variable_type = 1)
    return
end

function test_unrecognized_kwarg()
    model = Model()
    err = ErrorException(
        "In `@variable(model, x, foo = 1)`: " *
        "Unrecognized keyword argument: foo.\n\nIf you're trying to " *
        "create a JuMP extension, you need to implement " *
        "`build_variable`. Read the docstring for more details.",
    )
    @test_throws_strip err @variable(model, x, foo = 1)
    return
end

function test_kwargs()
    model = Model()
    x = @variable(model; integer = true)
    @test is_integer(x)
    x = @variable(model; integer = false)
    @test !is_integer(x)
    @variable(model, y; integer = true)
    @test is_integer(y)
    return
end

function test_Model_as_index()
    m = Model()
    @variable(m, x)
    msg = "Index m is the same symbol as the model. Use a different name for the index."
    @test_macro_throws(
        ErrorException("In `@variable(m, y[m = 1:2] <= m)`: $(msg)"),
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
        ErrorException("In `@NLparameter(m, p[m = 1:2] == m)`: $(msg)"),
        @NLparameter(m, p[m = 1:2] == m),
    )
    return
end

function test_Broadcasting_sparse_arrays()
    model = Model()
    A = SparseArrays.sparse([1, 2], [1, 2], [1, 1])
    b = SparseArrays.sparsevec([1, 2], [1, 2])
    @variable(model, x[1:2])
    le_cons = @constraint(model, A * x .<= b)
    @test length(le_cons) ==
          num_constraints(model, AffExpr, MOI.LessThan{Float64})
    int_cons = @constraint(model, -b .<= A * x .<= b)
    @test length(int_cons) ==
          num_constraints(model, AffExpr, MOI.Interval{Float64})
    return
end

function test_singular_plural_error()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, begin\n    x\nend)`: " *
            "Invalid syntax. Did you mean to use `@variables`?",
        ),
        @variable(model, begin
            x
        end),
    )
    @variable(model, x)
    @test_macro_throws(
        ErrorException(
            "In `@constraint(model, begin\n    x >= 0\nend)`: " *
            "Invalid syntax. Did you mean to use `@constraints`?",
        ),
        @constraint(model, begin
            x >= 0
        end),
    )
    @test_macro_throws(
        ErrorException(
            "In `@expression(model, begin\n    x\nend)`: " *
            "Invalid syntax. Did you mean to use `@expressions`?",
        ),
        @expression(model, begin
            x
        end),
    )
    @test_macro_throws(
        ErrorException(
            "In `@NLconstraint(model, begin\n    x >= 0\nend)`: " *
            "Invalid syntax. Did you mean to use `@NLconstraints`?",
        ),
        @NLconstraint(model, begin
            x >= 0
        end),
    )
    @test_macro_throws(
        ErrorException(
            "In `@NLexpression(model, begin\n    x\nend)`: " *
            "Invalid syntax. Did you mean to use `@NLexpressions`?",
        ),
        @NLexpression(model, begin
            x
        end),
    )
    @test_macro_throws(
        ErrorException(
            "In `@NLparameter(model, begin\n    x == 1\nend)`: " *
            "Invalid syntax: did you mean to use `@NLparameters`?",
        ),
        @NLparameter(model, begin
            x == 1
        end),
    )
    return
end

function test_supports_shift_constant()
    model = Model()
    @variable(model, x)
    @constraint(model, c, x + 0.5 in MOI.Integer())
    obj = constraint_object(c)
    @test obj.func == x + 0.5
    @test obj.set == MOI.Integer()
    @test_throws ErrorException add_to_function_constant(c, 1.0)
    return
end

function test_Interval_errors()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "Interval constraint contains non-constant left- or right-hand " *
        "sides. Reformulate as two separate constraints, or move all " *
        "variables into the central term.",
    )
    @test_throws err @NLconstraint(model, x <= x <= 2x)
    return
end

function test_Interval_errors_lhs()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, x <= x <= 2)`: Interval constraint " *
        "contains non-constant left- or right-hand sides. Reformulate as " *
        "two separate constraints, or move all variables into the " *
        "central term.",
    )
    @test_throws_strip err @constraint(model, x <= x <= 2)
    return
end

function test_Interval_errors_rhs()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, 2 <= x <= x)`: Interval constraint " *
        "contains non-constant left- or right-hand sides. Reformulate as " *
        "two separate constraints, or move all variables into the " *
        "central term.",
    )
    @test_throws_strip err @constraint(model, 2 <= x <= x)
    return
end

function test_Interval_errors_lhs_and_rhs()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, x <= x <= x)`: Interval constraint " *
        "contains non-constant left- or right-hand sides. Reformulate as " *
        "two separate constraints, or move all variables into the " *
        "central term.",
    )
    @test_throws_strip err @constraint(model, x <= x <= x)
    return
end

function test_nlparameter_basic()
    model = Model()
    @NLparameter(model, p == 1)
    @test p isa NonlinearParameter
    return
end

function test_nlparameter_Vector()
    model = Model()
    @NLparameter(model, p[i = 1:3] == i)
    @test p isa Vector{NonlinearParameter}
    @test value(p[2]) == 2
    return
end

function test_nlparameter_DenseAxisArray()
    model = Model()
    @NLparameter(model, p[i = 2:3] == i)
    @test p isa Containers.DenseAxisArray{NonlinearParameter}
    @test value(p[2]) == 2
    return
end

function test_nlparameter_SparseAxisArray()
    model = Model()
    @NLparameter(model, p[i = 1:3; iseven(i)] == i)
    @test p isa Containers.SparseAxisArray{NonlinearParameter}
    @test value(p[2]) == 2
    return
end

function test_nlparameter_requested_container()
    model = Model()
    S = 1:3
    @NLparameter(model, p[i = S] == i, container = Array)
    @test p isa Vector{NonlinearParameter}
    @test value(p[2]) == 2
    return
end

function test_nlparameter_too_many_positional_args()
    msg = "Invalid syntax: too many positional arguments."
    model = Model()
    @test_macro_throws(
        ErrorException("In `@NLparameter(model, p == 1, Int)`: $(msg)"),
        @NLparameter(model, p == 1, Int),
    )
    return
end

function test_nlparameter_unsupported_keyword_args()
    msg = "Invalid syntax: unsupported keyword arguments."
    model = Model()
    @test_macro_throws(
        ErrorException("In `@NLparameter(model, p == 1, bad = false)`: $(msg)"),
        @NLparameter(model, p == 1, bad = false),
    )
    return
end

function test_nlparameter_invalid_syntax()
    msg = "Invalid syntax: expected syntax of form `param == value`."
    model = Model()
    @test_macro_throws(
        ErrorException("In `@NLparameter(model, p)`: $(msg)"),
        @NLparameter(model, p),
    )
    return
end

function test_nlparameter_anonymous()
    model = Model()
    p = @NLparameter(model, [i = 1:2] == i)
    @test p isa Vector{NonlinearParameter}
    @test length(p) == 2
    q = @NLparameter(model, value = 1)
    @test q isa NonlinearParameter
    return
end

function test_nlparameter_anonymous_error()
    msg = "Invalid syntax: no positional args allowed for anonymous parameters."
    model = Model()
    @test_macro_throws(
        ErrorException("In `@NLparameter(model, p, value = 1)`: $(msg)"),
        @NLparameter(model, p, value = 1),
    )
    return
end

function test_nlparameter_invalid_number()
    msg = "Parameter value is not a number."
    model = Model()
    @test_throws_strip(
        ErrorException("In `@NLparameter(model, p == :a)`: $(msg)"),
        @NLparameter(model, p == :a),
    )
    return
end

function test_nlparameter_register()
    model = Model()
    @NLparameter(model, p == 1)
    @test model[:p] === p
    @NLparameter(model, q[i = 1:3; isodd(i)] == i)
    @test model[:q] === q
    @test_throws ErrorException @NLparameter(model, p == 1)
    r = @NLparameter(model, value = 1)
    @test_throws KeyError model[:r]
    return
end

function test_variable_vector_lowerbound()
    msg = """
    Passing arrays as variable bounds without indexing them is not supported.

    Instead of:
    ```julia
    @variable(model, x[1:2] >= lb)
    ```
    use
    ```julia
    @variable(model, x[i=1:2] >= lb[i])
    ```
    or
    ```julia
    @variable(model, x[1:2])
    set_lower_bound.(x, lb)
    ```
    """
    model = Model()
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] >= [1, 2])`: $(msg)"),
        @variable(model, x[1:2] >= [1, 2]),
    )
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] .>= [1, 2])`: $(msg)"),
        @variable(model, x[1:2] .>= [1, 2]),
    )
    return
end

function test_variable_vector_upperbound()
    msg = """
    Passing arrays as variable bounds without indexing them is not supported.

    Instead of:
    ```julia
    @variable(model, x[1:2] <= ub)
    ```
    use
    ```julia
    @variable(model, x[i=1:2] <= ub[i])
    ```
    or
    ```julia
    @variable(model, x[1:2])
    set_upper_bound.(x, ub)
    ```
    """
    model = Model()
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] <= [1, 2])`: $(msg)"),
        @variable(model, x[1:2] <= [1, 2]),
    )
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] .<= [1, 2])`: $(msg)"),
        @variable(model, x[1:2] .<= [1, 2]),
    )
    return
end

function test_variable_vector_fixed()
    msg = """
    Passing arrays as variable bounds without indexing them is not supported.

    Instead of:
    ```julia
    @variable(model, x[1:2] == fx)
    ```
    use
    ```julia
    @variable(model, x[i=1:2] == fx[i])
    ```
    or
    ```julia
    @variable(model, x[1:2])
    fix.(x, fx)
    ```
    """
    model = Model()
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] == [1, 2])`: $(msg)"),
        @variable(model, x[1:2] == [1, 2]),
    )
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2] .== [1, 2])`: $(msg)"),
        @variable(model, x[1:2] .== [1, 2]),
    )
    return
end

function test_variable_vector_start()
    msg = """
    Passing arrays as variable starts without indexing them is not supported.

    Instead of:
    ```julia
    @variable(model, x[1:2], start = x0)
    ```
    use
    ```julia
    @variable(model, x[i=1:2], start = x0[i])
    ```
    or
    ```julia
    @variable(model, x[1:2])
    set_start_value.(x, x0)
    ```
    """
    model = Model()
    @test_throws_strip(
        ErrorException("In `@variable(model, x[1:2], start = [1, 2])`: $(msg)"),
        @variable(model, x[1:2], start = [1, 2]),
    )
    return
end

function test_variable_vector_interval()
    msg = """
    Passing arrays as variable bounds without indexing them is not supported.

    Instead of:
    ```julia
    @variable(model, x[1:2] <= ub)
    ```
    use
    ```julia
    @variable(model, x[i=1:2] <= ub[i])
    ```
    or
    ```julia
    @variable(model, x[1:2])
    set_upper_bound.(x, ub)
    ```
    """
    model = Model()
    @test_throws_strip(
        ErrorException(
            "In `@variable(model, 0 <= x[2:3, 3:4] <= rand(2, 2))`: $(msg)",
        ),
        @variable(model, 0 <= x[2:3, 3:4] <= rand(2, 2)),
    )
    @test_throws_strip(
        ErrorException(
            "In `@variable(model, 0 .<= x[2:3, 3:4] .<= rand(2, 2))`: $(msg)",
        ),
        @variable(model, 0 .<= x[2:3, 3:4] .<= rand(2, 2)),
    )
    return
end

function test_invalid_name_errors()
    model = Model()
    @test_macro_throws(
        ErrorException("Expression x.y cannot be used as a name."),
        @variable(model, x.y),
    )
    return
end

function test_invalid_name_errors_denseaxisarray()
    model = Model()
    @test_macro_throws(
        ErrorException("Expression x.y cannot be used as a name."),
        @variable(model, x.y[2:3, 1:2]),
    )
    return
end

function test_invalid_name_errors_sparseaxisarray()
    model = Model()
    @test_macro_throws(
        ErrorException("Expression x.y cannot be used as a name."),
        @variable(model, x.y[i = 1:3; isodd(i)]),
    )
    return
end

function test_invalid_variable_syntax()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, MyInfo(1))`: Invalid syntax: your syntax " *
            "is wrong, but we don't know why. Consult the documentation for " *
            "various ways to create variables in JuMP.",
        ),
        @variable(model, MyInfo(1)),
    )
    return
end

function test_broadcasting_variable_in_set()
    model = Model()
    @variable(model, z[1:2] in MOI.GreaterThan.([3.0, 2.0]))
    @test num_variables(model) == 2
    @test lower_bound.(z) == [3.0, 2.0]
    @test has_upper_bound.(z) == [false, false]
    @variable(model, w[1:2] in MOI.GreaterThan(3.0))
    @test num_variables(model) == 4
    @test lower_bound.(w) == [3.0, 3.0]
    @test has_upper_bound.(w) == [false, false]
    @variable(model, v[1:2] in [MOI.GreaterThan(3.0), MOI.LessThan(3.0)])
    @test num_variables(model) == 6
    @test lower_bound(v[1]) == 3.0
    @test has_upper_bound(v[1]) == false
    @test upper_bound(v[2]) == 3.0
    @test has_lower_bound(v[2]) == false
    @test_throws(
        ErrorException,
        @variable(model, k in MOI.GreaterThan.([3.0, 2.0])),
    )
    @test_throws(
        ErrorException,
        @variable(model, u[1:3] in MOI.GreaterThan.([3.0, 2.0])),
    )
    @test num_variables(model) == 6
    # SparseAxisArray
    @variable(model, b[i = 1:2, j = 1:2; i + j == 3] in MOI.GreaterThan(3.0))
    @test num_variables(model) == 8
    @variable(
        model,
        a[i = 1:2, j = 1:2; i + j == 3] in Containers.SparseAxisArray(
            Dict(
                (1, 2) => MOI.GreaterThan(3.0),
                (2, 1) => MOI.GreaterThan(3.0),
            ),
        )
    )
    @test num_variables(model) == 10
    # DenseAxisArray
    @variable(model, dense_a[i = ["x", "xx"]] in MOI.GreaterThan(3.0))
    @test num_variables(model) == 12
    @variable(model, dense_b[1:2, 2:4] in MOI.GreaterThan(3.0))
    @test num_variables(model) == 18
    return
end

function test_parse_constraint_head_error()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, {x == 0})`: " *
        "Unsupported constraint expression: we don't know how to parse " *
        "constraints containing expressions of type :braces.\n\nIf you are " *
        "writing a JuMP extension, implement " *
        "`parse_constraint_head(::Function, ::Val{:braces}, args...)",
    )
    @test_macro_throws(err, @constraint(model, {x == 0}))
    return
end

function test_parse_constraint_head_inconsistent_vectorize()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, 1 .<= [x, x] <= 2)`: " *
        "Operators are inconsistently vectorized.",
    )
    @test_macro_throws(err, @constraint(model, 1 .<= [x, x] <= 2))
    return
end

function test_parse_constraint_head_inconsistent_signs()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, 1 >= x <= 2)`: " *
        "unsupported mix of comparison operators `1 >= ... <= 2`.\n\n" *
        "Two-sided rows must of the form `1 <= ... <= 2` or `2 >= ... >= 1`.",
    )
    @test_macro_throws(err, @constraint(model, 1 >= x <= 2))
    return
end

function test_MA_Zero_objective()
    model = Model()
    @test @objective(model, Min, sum(i for i in 1:0)) === 0.0
    return
end

function test_MA_Zero_expression()
    model = Model()
    @test @expression(model, sum(i for i in 1:0)) === 0.0
    @expression(model, expr[j = 1:2], sum(i for i in j:0))
    @test expr == [0.0, 0.0]
    @test expr isa Vector{Float64}
    return
end

function test_reorder_keyword_arguments()
    model = Model()
    @variable(model, integer = true, x)
    @test is_integer(x)
    c = @constraint(model, base_name = "my_c", x <= 1)
    @test name(c) == "my_c"
    return
end

function test_vectorized_constraint_name()
    model = Model()
    @variable(model, x)
    c = @constraint(model, [i = 1:2], [x, x] .>= [i, i + 1], base_name = "my_c")
    @test name(c[1][1]) == "my_c[1]"
    @test name(c[1][2]) == "my_c[1]"
    @test name(c[2][1]) == "my_c[2]"
    @test name(c[2][2]) == "my_c[2]"
    return
end

function test_expr_index_vars()
    model = Model()
    S = [(1, 2), (3, 4)]
    @variable(model, x[(i, j) in S])
    @test x[(1, 2)] isa VariableRef
    @test x[(3, 4)] isa VariableRef
    @test length(x) == 2
    return
end

function test_broadcasted_SparseAxisArray_constraint()
    model = Model()
    @variable(model, u[i = 1:2, i:3])
    Containers.@container(c[i = 1:2, i:3], rand())
    cons = @constraint(model, u .<= c)
    @test cons isa Containers.SparseAxisArray
    @test length(cons) == 5
    return
end

function test_broadcasted_DenseAxisArray_constraint()
    model = Model()
    S = 1:2
    @variable(model, u[S])
    Containers.@container(c[S], rand())
    cons = @constraint(model, u .<= c)
    @test cons isa Containers.DenseAxisArray
    @test length(cons) == 2
    return
end

function test_set_string_name_variables()
    model = Model()
    @variable(model, w[i = 1:2], set_string_name = isodd(i))
    @test !isempty(name(w[1]))
    @test isempty(name(w[2]))
    @test model[:w] == w
    @variable(model, x, set_string_name = false)
    @test isempty(name(x))
    @test model[:x] == x
    @variable(model, y[1:2], set_string_name = false)
    @test all(isempty.(name.(y)))
    @test model[:y] == y
    flag = false
    @variable(model, z, set_string_name = flag)
    @test isempty(name(z))
    return
end

function test_set_string_name_constraints()
    model = Model()
    @variable(model, x)
    @constraint(model, c1[i = 1:2], x <= i, set_string_name = isodd(i))
    @test !isempty(name(c1[1]))
    @test isempty(name(c1[2]))
    @test model[:c1] == c1
    @constraint(model, c2, x <= 1, set_string_name = false)
    @test isempty(name(c2))
    @test model[:c2] == c2
    @constraint(model, c3[i = 1:2], x <= i, set_string_name = false)
    @test all(isempty.(name.(c3)))
    @test model[:c3] == c3
    flag = false
    @constraint(model, c4, x <= 1, set_string_name = flag)
    @test isempty(name(c4))
    return
end

function test_set_string_name_model()
    model = Model()
    set_string_names_on_creation(model, false)
    @variable(model, w[i = 1:2], set_string_name = isodd(i))
    @test !isempty(name(w[1]))
    @test isempty(name(w[2]))
    @test model[:w] == w
    @variable(model, x)
    @test isempty(name(x))
    @test model[:x] == x
    @variable(model, y[1:2])
    @test all(isempty.(name.(y)))
    @test model[:y] == y
    flag = false
    @variable(model, z, set_string_name = flag)
    @test isempty(name(z))
    return
end

function test_set_string_name_Symmetric()
    model = Model()
    x = @variable(model, [1:2, 1:2], Symmetric, set_string_name = false)
    @test name.(x) == fill("", 2, 2)
    @test sprint.(show, x) == ["_[1]" "_[2]"; "_[2]" "_[3]"]
    return
end

function test_set_string_name_PSDCone()
    model = Model()
    set_string_names_on_creation(model, false)
    x = @variable(model, [1:2, 1:2] in PSDCone())
    @test name.(x) == fill("", 2, 2)
    @test sprint.(show, x) == ["_[1]" "_[2]"; "_[2]" "_[3]"]
    return
end

function test_nonlinear_unicode_operators()
    model = Model()
    @variable(model, x)
    c1 = @NLconstraint(model, x^2 ≤ 1)
    @test c1 isa ConstraintRef
    c2 = @NLconstraint(model, x^2 ≥ 1)
    @test c2 isa ConstraintRef
    c3 = @NLconstraint(model, 1 ≤ x^2 ≤ 2)
    @test c3 isa ConstraintRef
    return
end

function test_bad_model_type()
    model = "Not a model"
    @test_throws(
        ErrorException(
            "Expected model to be a JuMP model, but it has type String",
        ),
        @variable(model, x),
    )
    return
end

function test_constraint_not_enough_arguments()
    model = Model()
    @test_macro_throws(
        ErrorException("In `@constraint(model)`: Not enough arguments"),
        @constraint(model),
    )
    return
end

function test_constraint_no_constraint_expression_detected()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "In `@constraint(model, x = 2)`: No constraint expression detected. " *
        "If you are trying to construct an equality constraint, use `==` " *
        "instead of `=`.",
    )
    @test_macro_throws(err, @constraint(model, x = 2))
    return
end

function test_objective_not_enough_arguments()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@objective(model, Min)`: needs three arguments: model, " *
            "objective sense (Max or Min) and expression.",
        ),
        @objective(model, Min),
    )
    return
end

function test_expression_not_enough_arguments()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@expression(model)`: needs at least two arguments.",
        ),
        @expression(model),
    )
    return
end

function test_build_constraint_invalid()
    model = Model()
    @variable(model, x)
    @test_macro_throws(
        ErrorException(
            "In `@build_constraint(x)`: Incomplete constraint specification " *
            "x. Are you missing a comparison (<=, >=, or ==)?",
        ),
        @build_constraint(x),
    )
    return
end

function test_variable_reverse_sense()
    model = Model()
    @variable(model, 1 <= a)
    @variable(model, 1 ≤ b)
    @variable(model, 1 >= c)
    @variable(model, 1 ≥ d)
    @variable(model, 1 == e)
    @test lower_bound(a) == lower_bound(b) == 1
    @test upper_bound(c) == upper_bound(d) == 1
    @test fix_value(e) == 1
    return
end

function test_variable_unsupported_operator()
    model = Model()
    @test_macro_throws(
        ErrorException("In `@variable(model, x ⊕ 1)`: unsupported operator ⊕"),
        @variable(model, x ⊕ 1),
    )
    return
end

function test_constraint_unsupported_operator()
    model = Model()
    @variable(model, x)
    @test_macro_throws(
        ErrorException(
            "In `@constraint(model, x ⊕ 1)`: unsupported operator ⊕",
        ),
        @constraint(model, x ⊕ 1),
    )
    return
end

function test_variable_anon_bounds()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, [1:2] >= 0)`: Cannot use explicit bounds " *
            "via >=, <= with an anonymous variable",
        ),
        @variable(model, [1:2] >= 0),
    )
    return
end

function test_nonlinear_flatten_expressions()
    model = Model()
    @variable(model, x[1:2, 1:2])
    a = @NLexpression(model, sum(x[i, j] for i in 1:2, j in 1:2))
    @test string(a) == "subexpression[1]: x[1,1] + x[1,2] + x[2,1] + x[2,2]"
    b = @NLexpression(model, sum(x[i, j] for i in 1:2 for j in 1:2))
    @test string(b) == "subexpression[2]: x[1,1] + x[1,2] + x[2,1] + x[2,2]"
    c = @NLexpression(model, sum(x[i, j] for i in 1:2 for j in i:2))
    @test string(c) == "subexpression[3]: x[1,1] + x[1,2] + x[2,2]"
    d = @NLexpression(
        model,
        sum(k * x[i, j] for i in 1:2 for j in i:2 for k in 2:3),
    )
    @test string(d) ==
          "subexpression[4]: 2.0 * x[1,1] + 3.0 * x[1,1] + 2.0 * x[1,2] + 3.0 * x[1,2] + 2.0 * x[2,2] + 3.0 * x[2,2]"
    return
end

function test_variable_Bool_argument()
    model = Model()
    err = ErrorException(
        "In `@variable(model, x, Bool)`: " *
        "Unsupported positional argument `Bool`. If you intended to create a " *
        "`{0, 1}` decision variable, use `Bin` instead. For example, " *
        "`@variable(model, x, Bin)` or `@variable(model, x, binary = true)`.",
    )
    @test_throws_strip(err, @variable(model, x, Bool))
    err = ErrorException(
        "In `@variable(model, x, Bool = true)`: " *
        "Unsupported keyword argument: Bool.\n\nIf you intended to " *
        "create a `{0, 1}` decision variable, use the `binary` keyword " *
        "argument instead: `@variable(model, x, binary = true)`.",
    )
    @test_throws_strip(err, @variable(model, x, Bool = true))
    return
end

function test_nonlinear_generator_init_sum()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, sum(x for i in 1:1; init = 0))
    @test string(a) == "subexpression[1]: +x"
    a = @NLexpression(model, sum(x for i in 1:0; init = 0))
    @test string(a) == "subexpression[2]: 0.0"
    a = @NLexpression(model, sum(x for i in 1:0; init = 1))
    @test string(a) == "subexpression[3]: 1.0"
    y = 3
    a = @NLexpression(model, sum(x for i in 1:0; init = y))
    @test string(a) == "subexpression[4]: 3.0"
    a = @NLexpression(model, sum(x for i in 1:0, j in 1:0; init = 1))
    @test string(a) == "subexpression[5]: 1.0"
    a = @NLexpression(model, sum(x for i in 1:2, j in 1:2; init = 3))
    @test string(a) == "subexpression[6]: 3.0 + x + x + x + x"
    a = @NLexpression(model, sum(x for i in 1:2, j in 1:0; init = 3))
    @test string(a) == "subexpression[7]: 3.0"
    a = @NLexpression(model, sum(x for i in 1:1; init = 2.5))
    @test string(a) == "subexpression[8]: 2.5 + x"
    return
end

function test_nonlinear_generator_init_prod()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, prod(x for i in 1:1; init = 0))
    @test string(a) == "subexpression[1]: 0.0 * x"
    a = @NLexpression(model, prod(x for i in 1:0; init = 1))
    @test string(a) == "subexpression[2]: 1.0"
    y = 3
    a = @NLexpression(model, prod(x for i in 1:0; init = y))
    @test string(a) == "subexpression[3]: 3.0"
    a = @NLexpression(model, prod(x for i in 1:0; init = 2))
    @test string(a) == "subexpression[4]: 2.0"
    a = @NLexpression(model, prod(x for i in 1:1; init = 2.5))
    @test string(a) == "subexpression[5]: 2.5 * x"
    return
end

function test_nonlinear_generator_init_min()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, minimum(x for i in 1:1; init = 0))
    @test string(a) == "subexpression[1]: min(0.0, x)"
    a = @NLexpression(model, minimum(x for i in 1:0; init = 1))
    @test string(a) == "subexpression[2]: 1.0"
    y = 3
    a = @NLexpression(model, minimum(x for i in 1:0; init = y))
    @test string(a) == "subexpression[3]: 3.0"
    a = @NLexpression(model, minimum(x for i in 1:0; init = 2))
    @test string(a) == "subexpression[4]: 2.0"
    a = @NLexpression(model, minimum(x for i in 1:1; init = -2))
    @test string(a) == "subexpression[5]: min(-2.0, x)"
    return
end

function test_nonlinear_generator_bad_init()
    model = Model()
    @variable(model, x)
    expr = :(sum((x for i in 1:1); bad_init = 3))
    @test_macro_throws(
        ErrorException("Unsupported nonlinear expression: $expr"),
        @NLexpression(model, sum(x for i in 1:1; bad_init = 3))
    )
    return
end

#!format: off

function test_nonlinear_generator_pos_init_sum()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, sum(x for i in 1:1, init = 1))
    @test string(a) == "subexpression[1]: 1.0 + x"
    a = @NLexpression(model, sum(x for i in 1:0, init = 0))
    @test string(a) == "subexpression[2]: 0.0"
    a = @NLexpression(model, sum(x for i in 1:0, init = 1))
    @test string(a) == "subexpression[3]: 1.0"
    y = 3
    a = @NLexpression(model, sum(x for i in 1:0, init = y))
    @test string(a) == "subexpression[4]: 3.0"
    return
end

function test_nonlinear_generator_pos_init_prod()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, prod(x for i in 1:1, init = 0))
    @test string(a) == "subexpression[1]: 0.0 * x"
    a = @NLexpression(model, prod(x for i in 1:0, init = 1))
    @test string(a) == "subexpression[2]: 1.0"
    y = 3
    a = @NLexpression(model, prod(x for i in 1:0, init = y))
    @test string(a) == "subexpression[3]: 3.0"
    a = @NLexpression(model, prod(x for i in 1:0, init = 2))
    @test string(a) == "subexpression[4]: 2.0"
    return
end

function test_nonlinear_generator_pos_init_min()
    model = Model()
    @variable(model, x)
    a = @NLexpression(model, minimum(x for i in 1:1, init = 0))
    @test string(a) == "subexpression[1]: min(0.0, x)"
    a = @NLexpression(model, minimum(x for i in 1:0, init = 1))
    @test string(a) == "subexpression[2]: 1.0"
    y = 3
    a = @NLexpression(model, minimum(x for i in 1:0, init = y))
    @test string(a) == "subexpression[3]: 3.0"
    a = @NLexpression(model, minimum(x for i in 1:0, init = 2))
    @test string(a) == "subexpression[4]: 2.0"
    return
end

#!format: on

function test_matrix_in_vector_set()
    model = Model()
    @variable(model, X[1:2, 1:2])
    A = [1 2; 3 4]
    @test_throws_strip(
        ErrorException(
            "In `@constraint(model, X >= A)`: " *
            "Unsupported matrix in vector-valued set. Did you mean to use the " *
            "broadcasting syntax `.>=` instead? Alternatively, perhaps you are " *
            "missing a set argument like `@constraint(model, X >= 0, PSDCone())` " *
            "or `@constraint(model, X >= 0, HermmitianPSDCone())`.",
        ),
        @constraint(model, X >= A),
    )
    @test_throws_strip(
        ErrorException(
            "In `@constraint(model, X <= A)`: " *
            "Unsupported matrix in vector-valued set. Did you mean to use the " *
            "broadcasting syntax `.<=` instead? Alternatively, perhaps you are " *
            "missing a set argument like `@constraint(model, X <= 0, PSDCone())` " *
            "or `@constraint(model, X <= 0, HermmitianPSDCone())`.",
        ),
        @constraint(model, X <= A),
    )
    @test_throws_strip(
        ErrorException(
            "In `@constraint(model, X == A)`: " *
            "Unsupported matrix in vector-valued set. Did you mean to use the " *
            "broadcasting syntax `.==` for element-wise equality? Alternatively, " *
            "this syntax is supported in the special case that the matrices are " *
            "`LinearAlgebra.Symmetric` or `LinearAlgebra.Hermitian`.",
        ),
        @constraint(model, X == A),
    )
    return
end

function test_hermitian_variable_tag()
    model = Model()
    @variable(model, x[1:3, 1:3], Hermitian)
    @test x isa LinearAlgebra.Hermitian
    @test num_variables(model) == 9
    return
end

function test_unsupported_operator_errors()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x > 0)`: " *
            "unsupported operator `>`.\n\n" *
            "JuMP does not support strict inequalities, use `>=` instead.\n\n" *
            "If you require a strict inequality, you will need to use a " *
            "tolerance. For example, instead of `x > 1`, do `x >= 1 + 1e-4`. " *
            "If the variable must take integer values, use a tolerance of " *
            "`1.0`. If the variable may take continuous values, note that this " *
            "work-around can cause numerical issues, and your bound may not " *
            "hold exactly.",
        ),
        @variable(model, x > 0),
    )
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, x < 0)`: " *
            "unsupported operator `<`.\n\n" *
            "JuMP does not support strict inequalities, use `<=` instead.\n\n" *
            "If you require a strict inequality, you will need to use a " *
            "tolerance. For example, instead of `x < 1`, do `x <= 1 - 1e-4`. " *
            "If the variable must take integer values, use a tolerance of " *
            "`1.0`. If the variable may take continuous values, note that this " *
            "work-around can cause numerical issues, and your bound may not " *
            "hold exactly.",
        ),
        @variable(model, x < 0),
    )
    @variable(model, x)
    @test_macro_throws(
        ErrorException(
            "In `@constraint(model, x > 0)`: " *
            "unsupported operator `>`.\n\n" *
            "JuMP does not support strict inequalities, use `>=` instead.\n\n" *
            "If you require a strict inequality, you will need to use a " *
            "tolerance. For example, instead of `x > 1`, do `x >= 1 + 1e-4`. " *
            "If the constraint must take integer values, use a tolerance of " *
            "`1.0`. If the constraint may take continuous values, note that this " *
            "work-around can cause numerical issues, and your constraint may not " *
            "hold exactly.",
        ),
        @constraint(model, x > 0),
    )
    @test_macro_throws(
        ErrorException(
            "In `@constraint(model, x < 0)`: " *
            "unsupported operator `<`.\n\n" *
            "JuMP does not support strict inequalities, use `<=` instead.\n\n" *
            "If you require a strict inequality, you will need to use a " *
            "tolerance. For example, instead of `x < 1`, do `x <= 1 - 1e-4`. " *
            "If the constraint must take integer values, use a tolerance of " *
            "`1.0`. If the constraint may take continuous values, note that this " *
            "work-around can cause numerical issues, and your constraint may not " *
            "hold exactly.",
        ),
        @constraint(model, x < 0),
    )
    return
end

function test_unsupported_ternary_operator()
    model = Model()
    @test_macro_throws(
        ErrorException(
            "In `@variable(model, 1 < x < 2)`: " *
            "unsupported mix of comparison operators `1 < ... < 2`.\n\n" *
            "Two-sided variable bounds must of the form `1 <= ... <= 2` or " *
            "`2 >= ... >= 1`.",
        ),
        @variable(model, 1 < x < 2),
    )
    return
end

end  # module
