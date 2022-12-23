#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/print.jl
# Testing for all pretty-printing-related functionality
#############################################################################

module TestPrint

using JuMP
using LinearAlgebra
using Test

# Helper function to test IO methods work correctly
function _io_test_show(::MIME"text/plain", obj, exp_str)
    @test sprint(show, obj) == exp_str
    return
end

function _io_test_show(::MIME"text/latex", obj, exp_str)
    @test sprint(show, "text/latex", obj) == string("\$\$ ", exp_str, " \$\$")
    return
end

function _io_test_print(::MIME"text/plain", obj, exp_str)
    @test sprint(print, obj) == exp_str
    return
end

function _io_test_print(::MIME"text/latex", obj, exp_str)
    @test sprint(show, "text/latex", obj) == string("\$\$ ", exp_str, " \$\$")
    return
end

function _io_test_print(::MIME"text/plain", obj::AbstractModel, exp_str)
    @test sprint(print, obj) == exp_str
    @test model_string(MIME("text/plain"), obj) == exp_str
    return
end

function _io_test_print(::MIME"text/latex", obj::AbstractModel, exp_str)
    model = latex_formulation(obj)
    @test sprint(io -> show(io, MIME("text/latex"), model)) ==
          string("\$\$ ", exp_str, " \$\$")
    @test model_string(MIME("text/latex"), obj) ==
          string("\$\$ ", exp_str, " \$\$")
    # TODO(odow): I don't know how to test an IJulia display without adding
    # IJulia as a test-dependency, so just print and check it doesn't error.
    print(obj)
    return
end

function io_test(mode, obj, exp_str; repl = :both)
    if repl == :show || repl == :both
        _io_test_show(mode, obj, exp_str)
    end
    if repl == :print || repl == :both
        _io_test_print(mode, obj, exp_str)
    end
    return
end

# Used to test that JuMP printing works correctly for types for which
# oneunit is not convertible to Float64
struct UnitNumber <: Number
    α::Float64
end
Base.zero(::Union{UnitNumber,Type{UnitNumber}}) = UnitNumber(0.0)
Base.oneunit(::Union{UnitNumber,Type{UnitNumber}}) = UnitNumber(1.0)
Base.:(+)(u::UnitNumber, v::UnitNumber) = UnitNumber(u.α + v.α)
Base.:(-)(u::UnitNumber, v::UnitNumber) = UnitNumber(u.α - v.α)
Base.:(*)(α::Float64, u::UnitNumber) = UnitNumber(α * u.α)
Base.abs(u::UnitNumber) = UnitNumber(abs(u.α))
Base.isless(u::UnitNumber, v::UnitNumber) = isless(u.α, v.α)

# Used to test extensibility of JuMP printing for `AbstractConstraint`
struct CustomConstraint{S<:AbstractShape} <: AbstractConstraint
    function_str::String
    in_set_str::String
    shape::S
end

function JuMP.function_string(mode, constraint::CustomConstraint)
    return constraint.function_str
end

function JuMP.in_set_string(mode, constraint::CustomConstraint)
    return constraint.in_set_str
end
struct CustomIndex
    value::Int
end

function JuMP.add_constraint(
    model::GenericModel,
    constraint::CustomConstraint,
    name::String,
)
    if !haskey(model.ext, :custom)
        model.ext[:custom_constraints] = CustomConstraint[]
        model.ext[:custom_names] = String[]
    end
    constraints = model.ext[:custom_constraints]
    push!(constraints, constraint)
    push!(model.ext[:custom_names], name)
    return ConstraintRef(
        model,
        CustomIndex(length(constraints)),
        constraint.shape,
    )
end

function JuMP.constraint_object(cref::ConstraintRef{Model,CustomIndex})
    return cref.model.ext[:custom_constraints][cref.index.value]
end

function JuMP.name(cref::ConstraintRef{Model,CustomIndex})
    return cref.model.ext[:custom_names][cref.index.value]
end

function test_printing_expressions()
    # Most of the expression logic is well covered by test/operator.jl
    # This is really just to check IJulia printing for expressions
    mod = Model()
    @variable(mod, x[1:5])
    @variable(mod, y[i = 2:4, j = i:5])
    @variable(mod, z)

    ex = @expression(mod, x[1] + 2 * y[2, 3])
    io_test(MIME("text/plain"), ex, "x[1] + 2 y[2,3]")
    io_test(MIME("text/latex"), ex, "x_{1} + 2 y_{2,3}")

    ex = @expression(mod, x[1] + 2 * y[2, 3] + x[1])
    io_test(MIME("text/plain"), ex, "2 x[1] + 2 y[2,3]")
    io_test(MIME("text/latex"), ex, "2 x_{1} + 2 y_{2,3}")

    # TODO: These tests shouldn't depend on order of the two variables in
    # quadratic terms, i.e., x*y vs y*x.
    ex = @expression(mod, (x[1] + x[2]) * (y[2, 2] + 3.0))
    io_test(
        MIME("text/plain"),
        ex,
        "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2]",
    )
    io_test(
        MIME("text/latex"),
        ex,
        "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + 3 x_{1} + 3 x_{2}",
    )

    ex = @expression(mod, (y[2, 2] + 3.0) * (x[1] + x[2]))
    io_test(
        MIME("text/plain"),
        ex,
        "y[2,2]*x[1] + y[2,2]*x[2] + 3 x[1] + 3 x[2]",
    )
    io_test(
        MIME("text/latex"),
        ex,
        "y_{2,2}\\times x_{1} + y_{2,2}\\times x_{2} + 3 x_{1} + 3 x_{2}",
    )

    ex = @expression(mod, (x[1] + x[2]) * (y[2, 2] + 3.0) + z^2 - 1)
    repl_sq = JuMP._math_symbol(MIME("text/plain"), :sq)
    io_test(
        MIME("text/plain"),
        ex,
        "x[1]*y[2,2] + x[2]*y[2,2] + z$repl_sq + 3 x[1] + 3 x[2] - 1",
    )
    ijulia_sq = JuMP._math_symbol(MIME("text/latex"), :sq)
    io_test(
        MIME("text/latex"),
        ex,
        "x_{1}\\times y_{2,2} + x_{2}\\times y_{2,2} + z$ijulia_sq + 3 x_{1} + 3 x_{2} - 1",
    )

    ex = @expression(mod, -z * x[1] - x[1] * z + x[1] * x[2] + 0 * z^2)
    io_test(MIME("text/plain"), ex, "-2 z*x[1] + x[1]*x[2]")
    io_test(MIME("text/latex"), ex, "-2 z\\times x_{1} + x_{1}\\times x_{2}")

    ex = z^2 + x[1] - z^2 - x[1]
    io_test(MIME("text/plain"), ex, "0 z² + 0 x[1]")
    io_test(MIME("text/latex"), ex, "0 z$ijulia_sq + 0 x_{1}")
    return
end

# See https://github.com/jump-dev/JuMP.jl/pull/1352
function test_printing_expression_unit_coefficient_type()
    m = Model()
    @variable m x
    @variable m y
    u = UnitNumber(2.0)
    z = UnitNumber(0.0)
    aff = GenericAffExpr(zero(u), x => u, y => zero(u))
    io_test(MIME("text/plain"), aff, "$u x + $z y")
    io_test(MIME("text/latex"), aff, "$u x + $z y")
    drop_zeros!(aff)
    io_test(MIME("text/plain"), aff, "$u x")
    io_test(MIME("text/latex"), aff, "$u x")
    quad = aff * x
    io_test(MIME("text/plain"), quad, "$u x² + $z")
    io_test(MIME("text/latex"), quad, "$u x^2 + $z")
    return
end

function test_nonlinear_expressions()
    model = Model()
    @variable(model, x)
    expr = @NLexpression(model, x + 1)
    io_test(MIME("text/plain"), expr, "subexpression[1]: x + 1.0")
    return
end

function test_nonlinear_parameters()
    model = Model()
    param = @NLparameter(model, value = 1.0)
    io_test(MIME("text/plain"), param, "parameter[1] == 1.0")
    return
end

function test_registered_nonlinear_parameters()
    model = Model()
    @NLparameter(model, param == 1.0)
    io_test(MIME("text/plain"), param, "param == 1.0")
    return
end

function test_printing_NLPEvaluator()
    model = Model()
    evaluator = NLPEvaluator(model)
    io_test(
        MIME("text/plain"),
        evaluator,
        "Nonlinear.Evaluator with available features:\n" *
        "  * :Grad\n" *
        "  * :Jac\n" *
        "  * :JacVec\n" *
        "  * :Hess\n" *
        "  * :HessVec\n" *
        "  * :ExprGraph",
    )
    return
end

function test_nonlinear_constraints()
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)

    model = Model()
    @variable(model, x)
    constr_le = @NLconstraint(model, sin(x) <= 1)
    constr_ge = @NLconstraint(model, sin(x) >= 1)
    constr_eq = @NLconstraint(model, sin(x) == 1)
    constr_range = @NLconstraint(model, 0 <= sin(x) <= 1)
    constr_exponent_1 = @NLconstraint(model, x^4 <= 1)
    constr_exponent_2 = @NLconstraint(model, x^40.23 <= 1)
    constr_exponent_3 = @NLconstraint(model, x^4 + x^3 + x^2 <= 1)
    constr_exponent_4 = @NLconstraint(model, x^-4 + x^-3 + x^-2 <= 1)
    constr_exponent_5 = @NLconstraint(model, x^(x * x) <= 1)
    constr_exponent_6 = @NLconstraint(model, x^(2x) <= 1)
    constr_exponent_7 = @NLconstraint(model, x^(2 * 5) <= 1)
    constr_exponent_8 = @NLconstraint(model, x^(x^2) <= 1)

    io_test(MIME("text/plain"), constr_le, "sin(x) - 1.0 $le 0")
    io_test(MIME("text/plain"), constr_ge, "sin(x) - 1.0 $ge 0")
    io_test(MIME("text/plain"), constr_eq, "sin(x) - 1.0 $eq 0")
    # Note: This is inconsistent with the "x in [-1, 1]" printing for
    # regular constraints.
    io_test(MIME("text/plain"), constr_range, "0 $le sin(x) $le 1")
    io_test(MIME("text/plain"), constr_exponent_1, "x ^ 4.0 - 1.0 $le 0")
    io_test(MIME("text/plain"), constr_exponent_2, "x ^ 40.23 - 1.0 $le 0")
    io_test(
        MIME("text/plain"),
        constr_exponent_3,
        "(x ^ 4.0 + x ^ 3.0 + x ^ 2.0) - 1.0 $le 0",
    )
    io_test(
        MIME("text/plain"),
        constr_exponent_4,
        "(x ^ -4.0 + x ^ -3.0 + x ^ -2.0) - 1.0 $le 0",
    )
    io_test(MIME("text/plain"), constr_exponent_5, "x ^ (x * x) - 1.0 $le 0")
    io_test(MIME("text/plain"), constr_exponent_6, "x ^ (2.0 * x) - 1.0 $le 0")
    io_test(
        MIME("text/plain"),
        constr_exponent_7,
        "x ^ (2.0 * 5.0) - 1.0 $le 0",
    )
    io_test(MIME("text/plain"), constr_exponent_8, "x ^ (x ^ 2.0) - 1.0 $le 0")

    io_test(MIME("text/latex"), constr_le, "sin(x) - 1.0 \\leq 0")
    io_test(MIME("text/latex"), constr_ge, "sin(x) - 1.0 \\geq 0")
    io_test(MIME("text/latex"), constr_eq, "sin(x) - 1.0 = 0")
    io_test(MIME("text/latex"), constr_range, "0 \\leq sin(x) \\leq 1")
    io_test(MIME("text/latex"), constr_exponent_1, "x ^ {4.0} - 1.0 \\leq 0")
    io_test(MIME("text/latex"), constr_exponent_2, "x ^ {40.23} - 1.0 \\leq 0")
    io_test(
        MIME("text/latex"),
        constr_exponent_3,
        "(x ^ {4.0} + x ^ {3.0} + x ^ {2.0}) - 1.0 \\leq 0",
    )
    io_test(
        MIME("text/latex"),
        constr_exponent_4,
        "(x ^ {-4.0} + x ^ {-3.0} + x ^ {-2.0}) - 1.0 \\leq 0",
    )
    io_test(MIME("text/latex"), constr_exponent_5, "x ^ {x * x} - 1.0 \\leq 0")
    io_test(
        MIME("text/latex"),
        constr_exponent_6,
        "x ^ {2.0 * x} - 1.0 \\leq 0",
    )
    io_test(
        MIME("text/latex"),
        constr_exponent_7,
        "x ^ {2.0 * 5.0} - 1.0 \\leq 0",
    )
    io_test(
        MIME("text/latex"),
        constr_exponent_8,
        "x ^ {x ^ {2.0}} - 1.0 \\leq 0",
    )
    return
end

function test_nonlinear_constraint_with_anon_param_and_subexpression()
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    model = Model()
    @variable(model, x)
    expr = @NLexpression(model, x + 1)
    param = @NLparameter(model, value = 1.0)
    constr = @NLconstraint(model, expr - param <= 0)
    io_test(
        MIME("text/plain"),
        constr,
        "(subexpression[1] - parameter[1]) - 0.0 $le 0",
    )
    io_test(
        MIME("text/latex"),
        constr,
        "(subexpression_{1} - parameter_{1}) - 0.0 \\leq 0",
    )
    return
end

function test_nonlinear_constraint_with_named_param_and_subexpression()
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    model = Model()
    @variable(model, x)
    expr = @NLexpression(model, x + 1)
    model[:param] = @NLparameter(model, param == 1.0)
    constr = @NLconstraint(model, expr - param <= 0)
    io_test(
        MIME("text/plain"),
        constr,
        "(subexpression[1] - param) - 0.0 $le 0",
    )
    io_test(
        MIME("text/latex"),
        constr,
        "(subexpression_{1} - param) - 0.0 \\leq 0",
    )
    return
end

function test_printing_custom_constraint()
    model = Model()
    function test_constraint(function_str, in_set_str, name)
        constraint = CustomConstraint(function_str, in_set_str, ScalarShape())
        cref = add_constraint(model, constraint, name)
        @test string(cref) == "$name : $function_str $in_set_str"
    end
    test_constraint("fun", "set", "name")
    test_constraint("a", "b", "c")
    return
end

function test_extension_printing_variable_ref(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, 0 <= x <= 2)
    @test name(x) == "x"
    io_test(MIME("text/plain"), x, "x")
    io_test(MIME("text/latex"), x, "x")

    set_name(x, "x2")
    @test name(x) == "x2"
    io_test(MIME("text/plain"), x, "x2")
    io_test(MIME("text/latex"), x, "x2")
    set_name(x, "")
    @test name(x) == ""
    if x isa GenericVariableRef
        io_test(MIME("text/plain"), x, "_[1]")
        io_test(MIME("text/latex"), x, "{\\_}_{1}")
    else
        io_test(MIME("text/plain"), x, "anon")
        io_test(MIME("text/latex"), x, "anon")
    end
    @variable(m, z[1:2, 3:5])
    @test name(z[1, 3]) == "z[1,3]"
    io_test(MIME("text/plain"), z[1, 3], "z[1,3]")
    io_test(MIME("text/latex"), z[1, 3], "z_{1,3}")
    @test name(z[2, 4]) == "z[2,4]"
    io_test(MIME("text/plain"), z[2, 4], "z[2,4]")
    io_test(MIME("text/latex"), z[2, 4], "z_{2,4}")
    @test name(z[2, 5]) == "z[2,5]"
    io_test(MIME("text/plain"), z[2, 5], "z[2,5]")
    io_test(MIME("text/latex"), z[2, 5], "z_{2,5}")
    @variable(m, w[3:9, ["red", "blue", "green"]])
    @test name(w[7, "green"]) == "w[7,green]"
    io_test(MIME("text/plain"), w[7, "green"], "w[7,green]")
    io_test(MIME("text/latex"), w[7, "green"], "w_{7,green}")
    rng = 2:5
    @variable(m, v[rng, rng, rng, rng, rng, rng, rng])
    a_v = v[4, 5, 2, 3, 2, 2, 4]
    @test name(a_v) == "v[4,5,2,3,2,2,4]"
    io_test(MIME("text/plain"), a_v, "v[4,5,2,3,2,2,4]")
    io_test(MIME("text/latex"), a_v, "v_{4,5,2,3,2,2,4}")
    return
end

function test_extension_printing_base_name_keyword_argument(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x, base_name = "foo")
    @variable(m, y[1:3], base_name = "bar")
    num = 123
    @variable(m, z[[:red, :blue]], base_name = "color_$num")
    @variable(m, v[1:2, 1:2], PSD, base_name = string("i", "$num", num))
    @variable(m, w[1:3, 1:3], Symmetric, base_name = "symm")
    io_test(MIME("text/plain"), x, "foo")
    io_test(MIME("text/latex"), x, "foo")
    io_test(MIME("text/plain"), y[2], "bar[2]")
    io_test(MIME("text/latex"), y[2], "bar_{2}")
    io_test(MIME("text/plain"), z[:red], "color_123[red]")
    io_test(MIME("text/latex"), z[:red], "color\\_123_{red}")
    io_test(MIME("text/plain"), v[2, 1], "i123123[1,2]")
    io_test(MIME("text/latex"), v[2, 1], "i123123_{1,2}")
    io_test(MIME("text/plain"), w[1, 3], "symm[1,3]")
    io_test(MIME("text/latex"), w[1, 3], "symm_{1,3}")
    return
end

function test_extension_printing_vectorofvariables(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    model = ModelType()
    @variable(model, x)
    @variable(model, y)
    zero_constr = @constraint(model, [x, y] in MOI.Zeros(2))
    io_test(
        MIME("text/plain"),
        zero_constr,
        "[x, y] $in_sym MathOptInterface.Zeros(2)",
    )
    # TODO: Test in IJulia mode and do nice printing for Zeros().
    return
end

function test_extension_printing_scalaraffinefunction_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    model = ModelType()
    @variable(model, x)
    @constraint(model, linear_le, x + 0 <= 1)
    @constraint(model, linear_ge, x + 0 >= 1)
    @constraint(model, linear_eq, x + 0 == 1)
    @constraint(model, linear_range, -1 <= x + 0 <= 1)
    linear_noname = @constraint(model, x + 0 <= 1)
    io_test(MIME("text/plain"), linear_le, "linear_le : x $le 1")
    io_test(MIME("text/plain"), linear_eq, "linear_eq : x $eq 1")
    io_test(
        MIME("text/plain"),
        linear_range,
        "linear_range : x $in_sym [-1, 1]",
    )
    io_test(MIME("text/plain"), linear_noname, "x $le 1")
    # io_test doesn't work here because constraints print with a mix of math
    # and non-math.
    @test sprint(show, "text/latex", linear_le) == "linear_le : \$ x \\leq 1 \$"
    @test sprint(show, "text/latex", linear_ge) == "linear_ge : \$ x \\geq 1 \$"
    @test sprint(show, "text/latex", linear_eq) == "linear_eq : \$ x = 1 \$"
    @test sprint(show, "text/latex", linear_range) ==
          "linear_range : \$ x \\in \\[-1, 1\\] \$"
    @test sprint(show, "text/latex", linear_noname) == "\$\$ x \\leq 1 \$\$"
    return
end

function test_extension_printing_vectoraffinefunction_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    model = ModelType()
    @variable(model, x)
    @constraint(model, soc_constr, [x - 1, x + 1] in SecondOrderCone())
    io_test(
        MIME("text/plain"),
        soc_constr,
        "soc_constr : " *
        "[x - 1, x + 1] $in_sym MathOptInterface.SecondOrderCone(2)",
    )
    # TODO: Test in IJulia mode.
    return
end

function test_extension_printing_scalarquadraticfunction_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    le = JuMP._math_symbol(MIME("text/plain"), :leq)
    sq = JuMP._math_symbol(MIME("text/plain"), :sq)

    model = ModelType()
    @variable(model, x)
    quad_constr = @constraint(model, 2x^2 <= 1)

    io_test(MIME("text/plain"), quad_constr, "2 x$sq $le 1")
    # TODO: Test in IJulia mode.
    return
end

function test_extension_printing_indicator_constraints(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    le = JuMP._math_symbol(MIME("text/plain"), :leq)

    model = ModelType()
    @variable(model, x, Bin)
    @variable(model, y)
    ind_constr = @constraint(model, !x => {y <= 1})

    io_test(MIME("text/plain"), ind_constr, "!x => {y $le 1}")
    # TODO: Test in IJulia mode.
    return
end

# Test printing of models of type `ModelType` for which the model is stored in
# an MOI backend
function test_model_printing()
    ModelType = Model
    repl(s) = JuMP._math_symbol(MIME("text/plain"), s)
    le, ge, eq, inset = repl(:leq), repl(:geq), repl(:eq), repl(:in)
    #------------------------------------------------------------------

    model_1 = ModelType()
    @variable(model_1, a >= 1)
    @variable(model_1, b <= 1)
    @variable(model_1, -1 <= c <= 1)
    @variable(model_1, a1 >= 1, Int)
    @variable(model_1, b1 <= 1, Int)
    @variable(model_1, -1 <= c1 <= 1, Int)
    @variable(model_1, x, Bin)
    @variable(model_1, y)
    @variable(model_1, z, Int)
    @variable(model_1, u[1:3], Bin)
    @variable(model_1, fi == 9)
    @objective(model_1, Max, a - b + 2a1 - 10x)
    @constraint(model_1, con, a + b - 10c + c1 - 2x <= 1)
    @constraint(model_1, a * b <= 2)
    @constraint(model_1, soc, [1 - a; u] in SecondOrderCone())
    @constraint(model_1, [a b; c x] in PSDCone())
    @constraint(model_1, Symmetric([a b; b x]) in PSDCone())
    @constraint(model_1, [a, b, c] in MOI.PositiveSemidefiniteConeTriangle(2))
    @constraint(model_1, [a, b, c, x] in MOI.PositiveSemidefiniteConeSquare(2))

    VariableType = typeof(a)

    io_test(
        MIME("text/plain"),
        model_1,
        """
Max a - b + 2 a1 - 10 x
Subject to
 con : a + b - 10 c + c1 - 2 x $le 1
 a*b $le 2
 [a  b;
  b  x] $inset $(PSDCone())
 [a, b, c] $inset $(MOI.PositiveSemidefiniteConeTriangle(2))
 [a  b;
  c  x] $inset $(PSDCone())
 [a, b, c, x] $inset $(MOI.PositiveSemidefiniteConeSquare(2))
 soc : [-a + 1, u[1], u[2], u[3]] $inset $(MOI.SecondOrderCone(4))
 fi $eq 9
 a $ge 1
 c $ge -1
 a1 $ge 1
 c1 $ge -1
 b $le 1
 c $le 1
 b1 $le 1
 c1 $le 1
 a1 integer
 b1 integer
 c1 integer
 z integer
 x binary
 u[1] binary
 u[2] binary
 u[3] binary
""";
        repl = :print,
    )

    io_test(
        MIME("text/plain"),
        model_1,
        """
A JuMP Model
Maximization problem with:
Variables: 13
Objective function type: $(GenericAffExpr{Float64,VariableType})
`$(GenericAffExpr{Float64,VariableType})`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
`$(GenericQuadExpr{Float64,VariableType})`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
`$(Array{VariableType,1})`-in-`MathOptInterface.PositiveSemidefiniteConeTriangle`: 2 constraints
`$(Array{VariableType,1})`-in-`MathOptInterface.PositiveSemidefiniteConeSquare`: 2 constraints
`$(Array{GenericAffExpr{Float64,VariableType},1})`-in-`MathOptInterface.SecondOrderCone`: 1 constraint
`$VariableType`-in-`MathOptInterface.EqualTo{Float64}`: 1 constraint
`$VariableType`-in-`MathOptInterface.GreaterThan{Float64}`: 4 constraints
`$VariableType`-in-`MathOptInterface.LessThan{Float64}`: 4 constraints
`$VariableType`-in-`MathOptInterface.Integer`: 4 constraints
`$VariableType`-in-`MathOptInterface.ZeroOne`: 4 constraints
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: a, a1, b, b1, c, c1, con, fi, soc, u, x, y, z""";
        repl = :show,
    )

    io_test(
        MIME("text/latex"),
        model_1,
        "\\begin{aligned}\n" *
        "\\max\\quad & a - b + 2 a1 - 10 x\\\\\n" *
        "\\text{Subject to} \\quad & a + b - 10 c + c1 - 2 x \\leq 1\\\\\n" *
        " & a\\times b \\leq 2\\\\\n" *
        " & \\begin{bmatrix}\n" *
        "a & b\\\\\n" *
        "\\cdot & x\\\\\n" *
        "\\end{bmatrix} \\in \\text{$(PSDCone())}\\\\\n" *
        " & [a, b, c] \\in \\text{MathOptInterface.PositiveSemidefiniteConeTriangle(2)}\\\\\n" *
        " & \\begin{bmatrix}\n" *
        "a & b\\\\\n" *
        "c & x\\\\\n" *
        "\\end{bmatrix} \\in \\text{$(PSDCone())}\\\\\n" *
        " & [a, b, c, x] \\in \\text{MathOptInterface.PositiveSemidefiniteConeSquare(2)}\\\\\n" *
        " & [-a + 1, u_{1}, u_{2}, u_{3}] \\in \\text{MathOptInterface.SecondOrderCone(4)}\\\\\n" *
        " & fi = 9\\\\\n" *
        " & a \\geq 1\\\\\n" *
        " & c \\geq -1\\\\\n" *
        " & a1 \\geq 1\\\\\n" *
        " & c1 \\geq -1\\\\\n" *
        " & b \\leq 1\\\\\n" *
        " & c \\leq 1\\\\\n" *
        " & b1 \\leq 1\\\\\n" *
        " & c1 \\leq 1\\\\\n" *
        " & a1 \\in \\mathbb{Z}\\\\\n" *
        " & b1 \\in \\mathbb{Z}\\\\\n" *
        " & c1 \\in \\mathbb{Z}\\\\\n" *
        " & z \\in \\mathbb{Z}\\\\\n" *
        " & x \\in \\{0, 1\\}\\\\\n" *
        " & u_{1} \\in \\{0, 1\\}\\\\\n" *
        " & u_{2} \\in \\{0, 1\\}\\\\\n" *
        " & u_{3} \\in \\{0, 1\\}\\\\\n" *
        "\\end{aligned}";
        repl = :print,
    )

    #------------------------------------------------------------------

    model_2 = ModelType()
    @variable(model_2, x, Bin)
    @variable(model_2, y, Int)
    @constraint(model_2, x * y <= 1)
    io_test(
        MIME("text/plain"),
        model_2,
        """
A JuMP Model
Feasibility problem with:
Variables: 2
`$(GenericQuadExpr{Float64,VariableType})`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
`$VariableType`-in-`MathOptInterface.Integer`: 1 constraint
`$VariableType`-in-`MathOptInterface.ZeroOne`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x, y""";
        repl = :show,
    )

    model_3 = ModelType()
    @variable(model_3, x)
    @constraint(model_3, x <= 3)

    io_test(
        MIME("text/plain"),
        model_3,
        """
A JuMP Model
Feasibility problem with:
Variable: 1
`$(GenericAffExpr{Float64,VariableType})`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x""";
        repl = :show,
    )
    return
end

function test_printing_model_with_nonlinear()
    eq = JuMP._math_symbol(MIME("text/plain"), :eq)
    model = Model()
    @variable(model, x)
    @NLobjective(model, Max, sin(x))
    c = @NLexpression(model, cos(x))
    @NLconstraint(model, c == 0)

    io_test(
        MIME("text/plain"),
        model,
        """
A JuMP Model
Maximization problem with:
Variable: 1
Objective function type: Nonlinear
Nonlinear: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x""";
        repl = :show,
    )

    io_test(
        MIME("text/plain"),
        model,
        """
Max sin(x)
Subject to
 subexpression[1] - 0.0 $eq 0
With NL expressions
 subexpression[1]: cos(x)
""";
        repl = :print,
    )

    io_test(
        MIME("text/latex"),
        model,
        "\\begin{aligned}\n" *
        "\\max\\quad & sin(x)\\\\\n" *
        "\\text{Subject to} \\quad & subexpression_{1} - 0.0 = 0\\\\\n" *
        "\\text{With NL expressions} \\quad & subexpression_{1}: cos(x)\\\\\n" *
        "\\end{aligned}";
        repl = :print,
    )
    return
end

function test_SingleVariable_constraints()
    ge = JuMP._math_symbol(MIME("text/plain"), :geq)
    model = Model()
    @variable(model, x >= 10)
    zero_one = @constraint(model, x in MOI.ZeroOne())
    io_test(MIME("text/plain"), LowerBoundRef(x), "x $ge 10")
    io_test(MIME("text/plain"), zero_one, "x binary")
    # TODO: Test in IJulia mode
    return
end

function test_latex_feasibility()
    io_test(
        MIME("text/latex"),
        Model(),
        "\\begin{aligned}\n\\text{feasibility}\\\\\n\\end{aligned}";
        repl = :print,
    )
    return
end

function test_latex_min()
    model = Model()
    @variable(model, x)
    @objective(model, Min, x)
    io_test(
        MIME("text/latex"),
        model,
        "\\begin{aligned}\n\\min\\quad & x\\\\\n\\end{aligned}";
        repl = :print,
    )
    return
end

function test_print_IJulia_with_math_operators()
    model = Model()
    @variable(model, x_ab)
    io_test(MIME("text/latex"), x_ab, "x\\_ab")
    @variable(model, y[1:2], base_name = "y[:a]")
    io_test(MIME("text/latex"), y[1], "y[:a]_{1}")
    io_test(MIME("text/latex"), y[2], "y[:a]_{2}")
    @variable(model, z, base_name = "z^1")
    io_test(MIME("text/latex"), z, "z\\^1")
    return
end

struct _UnsupportedNameOptimizer <: MOI.AbstractOptimizer end
MOI.is_empty(::_UnsupportedNameOptimizer) = true

function test_Name_direct_mode()
    model = direct_model(_UnsupportedNameOptimizer())
    @test name(model) == "A JuMP Model"
    return
end

function test_print_summary_min_sense()
    model = Model()
    @variable(model, x)
    @objective(model, Min, x)
    @test occursin("Minimization problem with:", sprint(show, model))
end

function test_show_latex_parameter()
    model = Model()
    @NLparameter(model, p == 1)
    @test sprint((io, p) -> show(io, MIME("text/latex"), p), p) == "p == 1.0"
    return
end

function test_complex_expr()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    f = 1.0im * x + 1.0im
    @test sprint(show, im * f) == "-x - 1"
    @test sprint(show, -f) == "-x im - im"
    @test sprint(show, f * x) == "x² im + x im"
    @test sprint(show, f * y) == "x*y im + y im"
    @test sprint(show, f + y) == "x im + y + im"
    @test sprint(show, 2f) == "2im x + 2im"
    return
end

function test_print_hermitian_psd_cone()
    model = Model()
    @variable(model, x[1:2])
    in_sym = JuMP._math_symbol(MIME("text/plain"), :in)
    H = Hermitian([x[1] 1im; -1im x[2]])
    c = @constraint(model, H in HermitianPSDCone())
    @test sprint(io -> show(io, MIME("text/plain"), c)) ==
          "[x[1]  im;\n -im   x[2]] $in_sym $(HermitianPSDCone())"
    @test sprint(io -> show(io, MIME("text/latex"), c)) ==
          "\$\$ \\begin{bmatrix}\nx_{1} & im\\\\\n-im & x_{2}\\\\\n\\end{bmatrix} \\in \\text{$(HermitianPSDCone())} \$\$"
    return
end

end
