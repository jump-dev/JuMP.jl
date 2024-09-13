#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestNLPExpr

using JuMP
using Test

import LinearAlgebra
import MutableArithmetics as MA

function test_extension_univariate_operators(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
        if f in (:+, :-, :abs2)
            op = getfield(Base, f)
            @test op(sin(x)) isa GenericNonlinearExpr{VariableRefType}
        elseif isdefined(Base, f)
            op = getfield(Base, f)
            @test op(x) isa GenericNonlinearExpr{VariableRefType}
        elseif isdefined(MOI.Nonlinear.SpecialFunctions, f)
            op = getfield(MOI.Nonlinear.SpecialFunctions, f)
            @test op(x) isa GenericNonlinearExpr{VariableRefType}
        end
    end
    return
end

function test_extension_binary_operators(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    num, aff, quad, nlp = 1.0, 1.0 + x, x^2, sin(x)
    for op in (+, -, *, /), a in (num, x, aff, quad, nlp)
        @test op(a, nlp) isa GenericNonlinearExpr{VariableRefType}
        @test op(nlp, a) isa GenericNonlinearExpr{VariableRefType}
    end
    for op in (*, /), a in (x, aff)
        @test op(a, quad) isa GenericNonlinearExpr{VariableRefType}
        @test op(quad, a) isa GenericNonlinearExpr{VariableRefType}
    end
    for a in (num, x, aff, quad), b in (x, aff, quad)
        @test /(a, b) isa GenericNonlinearExpr{VariableRefType}
    end
    return
end

function test_extension_objective(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @objective(model, Min, 2.0 * sin(x)^2 + cos(x) / x)
    @test objective_function(model) isa GenericNonlinearExpr{VariableRefType}
    return
end

function test_extension_expression(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @variable(model, y[1:3])
    @test string(@expression(model, *(y...))) == "(y[1]*y[2]) * y[3]"
    @test string(@expression(model, sin(x))) == "sin(x)"
    @test string(@expression(model, ifelse(x >= 0, x, 0))) ==
          "ifelse(x >= 0, x, 0)"
    @test string(@expression(model, 2^x)) == "2.0 ^ x"
    @test string(@expression(model, x^x)) == "x ^ x"
    @test string(@expression(model, sin(x)^2)) == "sin(x) ^ 2.0"
    @test string(@expression(model, sin(x)^2.0)) == "sin(x) ^ 2.0"
    @test string(@expression(model, 2 * sin(x)^2.0)) == "2.0 * (sin(x) ^ 2.0)"
    @test string(@expression(model, 1 + sin(x))) == "1.0 + sin(x)"
    @test string(@expression(model, 1 + 2 * sin(x))) == "1.0 + (2.0 * sin(x))"
    @test string(@expression(model, 2.0 * sin(x)^2 + cos(x) / x)) ==
          "(2.0 * (sin(x) ^ 2.0)) + (cos(x) / x)"
    @test string(@expression(model, 2.0 * sin(x)^2 - cos(x) / x)) ==
          "(2.0 * (sin(x) ^ 2.0)) - (cos(x) / x)"
    return
end

function test_extension_flatten_nary(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    expr_plus = GenericNonlinearExpr{VariableRefType}(:+, Any[x])
    expr_mult = GenericNonlinearExpr{VariableRefType}(:*, Any[x])
    expr_sin = GenericNonlinearExpr{VariableRefType}(:sin, Any[x])
    to_string(x) = string(flatten!(x))
    @test to_string(+(expr_plus, 1)) == "x + 1.0"
    @test to_string(+(1, expr_plus)) == "1.0 + x"
    @test to_string(+(expr_plus, x)) == "x + x"
    @test to_string(+(expr_sin, x)) == "sin(x) + x"
    @test to_string(+(x, expr_plus)) == "x + x"
    @test to_string(+(x, expr_sin)) == "x + sin(x)"
    @test to_string(+(expr_plus, expr_plus)) == "x + x"
    @test to_string(+(expr_plus, expr_sin)) == "x + sin(x)"
    @test to_string(+(expr_sin, expr_plus)) == "sin(x) + x"
    @test to_string(+(expr_sin, expr_sin)) == "sin(x) + sin(x)"
    @test to_string(*(expr_mult, 2)) == "x * 2.0"
    @test to_string(*(2, expr_mult)) == "2.0 * x"
    @test to_string(*(expr_mult, x)) == "x * x"
    @test to_string(*(expr_sin, x)) == "sin(x) * x"
    @test to_string(*(x, expr_mult)) == "x * x"
    @test to_string(*(x, expr_sin)) == "x * sin(x)"
    @test to_string(*(expr_mult, expr_mult)) == "x * x"
    @test to_string(*(expr_mult, expr_sin)) == "x * sin(x)"
    @test to_string(*(expr_sin, expr_mult)) == "sin(x) * x"
    @test to_string(*(expr_sin, expr_sin)) == "sin(x) * sin(x)"
    @test to_string(sin(+(expr_plus, 1))) == "sin(x + 1.0)"
    @test to_string(sin(*(expr_mult, expr_mult))) == "sin(x * x)"
    return
end

function test_extension_zero_one(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    @test string(zero(GenericNonlinearExpr{VariableRefType})) == "+(0.0)"
    @test string(one(GenericNonlinearExpr{VariableRefType})) == "+(1.0)"
    return
end

function test_extension_latex(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x)
    @test function_string(MIME("text/latex"), sin(x)) ==
          raw"\textsf{sin}\left({x}\right)"
    @test function_string(MIME("text/latex"), sin(x)^x) ==
          raw"{\textsf{sin}\left({x}\right)} ^ {x}"
    @test function_string(MIME("text/latex"), sin(x)^(x + 1)) ==
          raw"{\textsf{sin}\left({x}\right)} ^ {\left(x + 1\right)}"
    @test function_string(MIME("text/latex"), (x + 1)^x) ==
          raw"{\left(x + 1\right)} ^ {x}"
    @test function_string(MIME("text/latex"), (x + 1)^(x + 1)) ==
          raw"{\left(x + 1\right)} ^ {\left(x + 1\right)}"
    @test function_string(MIME("text/latex"), (x + 1)^sin(x)) ==
          raw"{\left(x + 1\right)} ^ {\textsf{sin}\left({x}\right)}"

    @expression(model, g, ifelse(x > 0, sin(x), x + cos(x)^2))
    @test function_string(MIME("text/latex"), g) ==
          raw"\textsf{ifelse}\left({{x} > {0}}, {\textsf{sin}\left({x}\right)}, {{x} + {\left({\textsf{cos}\left({x}\right)} ^ {2.0}\right)}}\right)"
    return
end

function test_extension_latex2(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x[1:2])
    @test function_string(MIME("text/latex"), sin(x[1])) ==
          raw"\textsf{sin}\left({x_{1}}\right)"
    @test function_string(MIME("text/latex"), sin(x[1])^x[2]) ==
          raw"{\textsf{sin}\left({x_{1}}\right)} ^ {x_{2}}"
    @expression(
        model,
        expr,
        (x[1] <= x[2]) || (x[1] >= x[2]) && (x[1] == x[1]),
    )
    @test function_string(MIME("text/latex"), expr) ==
          raw"{\left({x_{1}} \le {x_{2}}\right)} \vee {\left({\left({x_{1}} \ge {x_{2}}\right)} \wedge {\left({x_{1}} = {x_{1}}\right)}\right)}"
    return
end

function test_extension_expression_addmul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test string(@expression(model, x + 3 * sin(x))) == "x + (3.0 * sin(x))"
    @test string(@expression(model, 2 * x + 3 * sin(x))) ==
          "(2 x) + (3.0 * sin(x))"
    @test string(@expression(model, x^2 + 3 * sin(x))) ==
          "($(x^2)) + (3.0 * sin(x))"
    @test string(@expression(model, sin(x) + 3 * sin(x))) ==
          "sin(x) + (3.0 * sin(x))"
    @test string(@expression(model, sin(x) + 3 * x)) == "sin(x) + (3 x)"
    @test string(@expression(model, sin(x) + 3 * x * x)) ==
          "sin(x) + (3 $(x^2))"
    return
end

function test_extension_expression_explicit_add_mul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    f = sin(x)
    @test string(MA.operate!!(MA.add_mul, 1, 2, f)) == "1.0 + (2.0 * $f)"
    @test string(MA.operate!!(MA.add_mul, 1, f, 2)) == "1.0 + ($f * 2.0)"
    @test string(MA.operate!!(MA.add_mul, 1, f, f)) == "1.0 + ($f * $f)"
    @test string(MA.operate!!(MA.add_mul, f, 2, f)) == "$f + (2.0 * $f)"
    @test string(MA.operate!!(MA.add_mul, f, f, 2)) == "$f + ($f * 2.0)"
    @test string(MA.operate!!(MA.add_mul, f, f, f)) == "$f + ($f * $f)"
    return
end

function test_extension_expression_submul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test string(@expression(model, x - 3 * sin(x))) == "x - (3.0 * sin(x))"
    @test string(@expression(model, 2 * x - 3 * sin(x))) ==
          "(2 x) - (3.0 * sin(x))"
    @test string(@expression(model, x^2 - 3 * sin(x))) ==
          "($(x^2)) - (3.0 * sin(x))"
    @test string(@expression(model, sin(x) - 3 * sin(x))) ==
          "sin(x) - (3.0 * sin(x))"
    @test string(@expression(model, sin(x) - 3 * x)) == "sin(x) - (3 x)"
    @test string(@expression(model, sin(x) - 3 * x * x)) ==
          "sin(x) - (3 $(x^2))"
    return
end

function test_extension_aff_expr_convert(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    _to_string(x) = string(convert(GenericNonlinearExpr{VariableRefType}, x))
    @test _to_string(AffExpr(0.0)) == "+(0.0)"
    @test _to_string(AffExpr(1.0)) == "+(1.0)"
    @test _to_string(x + 1) == "x + 1.0"
    @test _to_string(2x + 1) == "(2.0 * x) + 1.0"
    @test _to_string(2x) == "2.0 * x"
    return
end

function test_extension_quad_expr_convert(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    _to_string(x) = string(convert(GenericNonlinearExpr{VariableRefType}, x))
    @test _to_string(QuadExpr(AffExpr(0.0))) == "+(0.0)"
    @test _to_string(QuadExpr(AffExpr(1.0))) == "+(1.0)"
    @test _to_string(x^2 + 1) == "(x * x) + 1.0"
    @test _to_string(2x^2 + 1) == "(2.0 * x * x) + 1.0"
    @test _to_string(2x^2) == "2.0 * x * x"
    @test _to_string(x^2 + x + 1) == "x + (x * x) + 1.0"
    @test _to_string(2x^2 + x + 1) == "x + (2.0 * x * x) + 1.0"
    @test _to_string(2x^2 + x) == "x + (2.0 * x * x)"
    @test _to_string(x^2 + 2x + 1) == "(2.0 * x) + (x * x) + 1.0"
    @test _to_string(2x^2 + 2x + 1) == "(2.0 * x) + (2.0 * x * x) + 1.0"
    @test _to_string(2x^2 + 2x) == "(2.0 * x) + (2.0 * x * x)"
    return
end

function test_extension_constraint_name(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, c, sin(x) <= 1)
    @test name(c) == "c"
    set_name(c, "d")
    @test name(c) == "d"
    @test startswith(string(c), "d : ")
    return
end

function test_extension_constraint_lessthan(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x <= 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.LessThan(zero(value_type(ModelType)))
    return
end

function test_extension_constraint_greaterthan(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x >= 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.GreaterThan(zero(value_type(ModelType)))
    return
end

function test_extension_constraint_equalto(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x == 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.EqualTo(zero(value_type(ModelType)))
    return
end

function test_extension_constraint_interval(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @constraint(model, c, 0 <= 2.0 * sin(x)^2 + cos(x) / x <= 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x)
    T = value_type(ModelType)
    @test obj.set == MOI.Interval(zero(T), one(T))
    return
end

function test_user_defined_function_overload()
    model = Model()
    @variable(model, x)
    f(x::Real) = x^2
    f(x::AbstractJuMPScalar) = NonlinearExpr(:f, x)
    register(model, :f, 1, f; autodiff = true)
    @test string(@expression(model, f(x))) == "f(x)"
    @test string(f(x) + f(x)) == "f(x) + f(x)"
    @test string(1 / (f(x) + f(x))) == "1.0 / (f(x) + f(x))"
    return
end

function test_extension_nonlinear_matrix_algebra(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, X[1:3, 1:3], Symmetric)
    @objective(model, Max, sum(X^4 .- X^3))
    @test objective_function(model) isa GenericNonlinearExpr{VariableRefType}
    return
end

"""
This test checks that we can work with expressions of arbitrary depth. Don't use
recursion!
"""
function test_extension_recursion_stackoverflow(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    expr = sin(x)
    for _ in 1:20_000
        expr = sin(expr)
    end
    @test @objective(model, Min, expr) isa GenericNonlinearExpr{VariableRefType}
    @test string(expr) isa String
    return
end

function test_nlobjective_with_nlexpr()
    model = Model()
    @variable(model, x)
    y = sin(x)
    @NLobjective(model, Min, y^2)
    nlp = nonlinear_model(model)
    @test isequal_canonical(jump_function(model, nlp.objective), sin(x)^2)
    return
end

function test_nlconstraint_with_nlexpr()
    model = Model()
    @variable(model, x)
    y = sin(x)
    @NLconstraint(model, c, y^2 <= 1)
    nlp = nonlinear_model(model)
    @test isequal_canonical(
        jump_function(model, nlp.constraints[index(c)].expression),
        sin(x)^2 - 1,
    )
    return
end

function test_constraint_object()
    model = Model()
    @variable(model, x)
    y = sin(x)
    @NLconstraint(model, c, y^2 <= 1)
    con = constraint_object(c)
    @test isequal_canonical(con.func, sin(x)^2 - 1.0)
    @test con.set == MOI.LessThan(0.0)
    return
end

function test_extension_expr_mle(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    data = [1.0, 2.0, 4.0, 8.0]
    n = length(data)
    @variable(model, x)
    @variable(model, y)
    obj = @expression(
        model,
        n / 2 * log(1 / (2 * y^2)) -
        sum((data[i] - x)^2 for i in 1:n) / (2 * y^2)
    )
    @test string(obj) ==
          "(2.0 * log(1.0 / (2 $(y^2)))) - ((4 $(x^2) - 30 x + 85) / (2 $(y^2)))"
    return
end

function test_extension_nl_macro(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test isequal_canonical(
        @expression(model, ifelse(x, 1, 2)),
        GenericNonlinearExpr{VariableRefType}(:ifelse, Any[x, 1, 2]),
    )
    @test isequal_canonical(
        @expression(model, x || 1),
        GenericNonlinearExpr{VariableRefType}(:||, Any[x, 1]),
    )
    @test isequal_canonical(
        @expression(model, x && 1),
        GenericNonlinearExpr{VariableRefType}(:&&, Any[x, 1]),
    )
    @test isequal_canonical(
        @expression(model, x < 0),
        GenericNonlinearExpr{VariableRefType}(:<, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x > 0),
        GenericNonlinearExpr{VariableRefType}(:>, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x <= 0),
        GenericNonlinearExpr{VariableRefType}(:<=, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x >= 0),
        GenericNonlinearExpr{VariableRefType}(:>=, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x == 0),
        GenericNonlinearExpr{VariableRefType}(:(==), Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, 0 < x <= 1),
        GenericNonlinearExpr{VariableRefType}(
            :&&,
            Any[@expression(model, 0 < x), @expression(model, x <= 1)],
        ),
    )
    @test isequal_canonical(
        @expression(model, ifelse(x > 0, x^2, sin(x))),
        GenericNonlinearExpr{VariableRefType}(
            :ifelse,
            Any[@expression(model, x > 0), x^2, sin(x)],
        ),
    )
    return
end

function test_register_univariate()
    model = Model()
    @variable(model, x)
    @operator(model, f, 1, x -> x^2)
    @test f isa NonlinearOperator
    @test sprint(show, f) == "NonlinearOperator($(f.func), :f)"
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_eval_non_jump()
    model = Model()
    @variable(model, x)
    @operator(model, f, 1, x -> x^2)
    @test f(2.0) == 4.0
    @operator(model, g, 2, (x, y) -> x^2 - sin(y))
    @test g(2.0, 3.0) == 4.0 - sin(3.0)
    return
end

function test_register_univariate_gradient()
    model = Model()
    @variable(model, x)
    @operator(model, f, 1, x -> x^2, x -> 2 * x)
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_univariate_gradient_hessian()
    model = Model()
    @variable(model, x)
    @operator(model, f, 1, x -> x^2, x -> 2 * x, x -> 2.0)
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_multivariate()
    model = Model()
    @variable(model, x[1:2])
    f = (x...) -> sum(x .^ 2)
    @operator(model, foo, 2, f)
    @test isequal_canonical(@expression(model, foo(x...)), foo(x...))
    @test isequal_canonical(foo(x...), NonlinearExpr(:foo, Any[x...]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:foo, 2) in attrs
    return
end

function test_register_multivariate_gradient()
    model = Model()
    @variable(model, x[1:2])
    f = (x...) -> sum(x .^ 2)
    ∇f = (g, x...) -> (g .= 2 .* x)
    @operator(model, foo, 2, f, ∇f)
    @test isequal_canonical(@expression(model, foo(x...)), foo(x...))
    @test isequal_canonical(foo(x...), NonlinearExpr(:foo, Any[x...]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:foo, 2) in attrs
    return
end

function test_register_multivariate_gradient_hessian()
    model = Model()
    @variable(model, x[1:2])
    f = (x...) -> sum(x .^ 2)
    ∇f = (g, x...) -> (g .= 2 .* x)
    ∇²f = (H, x...) -> begin
        for i in 1:2
            H[i, i] = 2.0
        end
    end
    @operator(model, foo, 2, f, ∇f, ∇²f)
    @test isequal_canonical(@expression(model, foo(x...)), foo(x...))
    @test isequal_canonical(foo(x...), NonlinearExpr(:foo, Any[x...]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:foo, 2) in attrs
    return
end

function test_register_multivariate_many_args()
    model = Model()
    @variable(model, x[1:10])
    f = (x...) -> sum(x .^ 2)
    @operator(model, foo, 10, f)
    @test isequal_canonical(foo(x...), NonlinearExpr(:foo, Any[x...]))
    @test foo((1:10)...) == 385
    return
end

function test_register_errors()
    model = Model()
    f = x -> x^2
    @test_throws(
        ErrorException(
            "Unable to add operator foo: invalid number of " *
            "functions provided. Got 4, but expected 1 (if function only), " *
            "2 (if function and gradient), or 3 (if function, gradient, and " *
            "hesssian provided)",
        ),
        @operator(model, foo, 2, f, f, f, f),
    )
    return
end

function test_value_expression()
    model = Model()
    @variable(model, x)
    f = x -> 1.1
    @test value(f, sin(x)) ≈ sin(1.1)
    @test value(f, sin(x) + cos(x)) ≈ sin(1.1) + cos(1.1)
    @test value(f, x^1.3 / x) ≈ 1.1^1.3 / 1.1
    @test value(f, @expression(model, ifelse(x > 1, 1, 2))) ≈ 1
    @test value(f, @expression(model, ifelse(x < 1, 1, 2))) ≈ 2
    @test value(f, @expression(model, ifelse(x < 1 || x > 2, 1, 2))) ≈ 2
    @test value(f, @expression(model, ifelse(x < 1 && x > 2, 1, 2))) ≈ 2
    @test value(f, sin(x + 1)) ≈ sin(1.1 + 1)
    @test value(f, sin(x^2 + x + 1)) ≈ sin(1.1^2 + 1.1 + 1)
    foo(x) = (x - 1)^2
    bar(x, y) = sqrt(x - y)
    @operator(model, my_foo, 1, foo)
    @operator(model, my_bar, 2, bar)
    @test value(f, my_foo(x)) ≈ (1.1 - 1)^2
    @test value(f, my_foo(x + 1)) ≈ (1.1 + 1 - 1)^2
    @test value(f, my_foo(x^2 + 1)) ≈ (1.1^2 + 1 - 1)^2
    @test value(f, my_foo(x^2 + x + 1)) ≈ (1.1^2 + 1.1 + 1 - 1)^2
    y = QuadExpr(x + 1)
    @test value(f, my_foo(y)) ≈ (value(f, y) - 1)^2
    @test value(f, my_bar(2.2, x)) ≈ sqrt(2.2 - 1.1)
    bad_udf = NonlinearOperator(f, :bad_udf)
    @test_throws(
        ErrorException(
            "Unable to evaluate nonlinear operator bad_udf because it was " *
            "not added as an operator.",
        ),
        value(f, bad_udf(x)),
    )
    return
end

function test_value_result()
    model = Model() do
        return MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    end
    @variable(model, x)
    optimize!(model)
    mock = unsafe_backend(model)
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.ResultCount(), 2)
    MOI.set(mock, MOI.PrimalStatus(1), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.PrimalStatus(2), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.VariablePrimal(1), optimizer_index(x), 1.1)
    MOI.set(mock, MOI.VariablePrimal(2), optimizer_index(x), 2.2)
    f = sin(x)
    @test value(f; result = 1) ≈ sin(1.1)
    @test value(f; result = 2) ≈ sin(2.2)
    return
end

function test_nonlinear_expr_owner_model()
    model = Model()
    @variable(model, x)
    f = NonlinearExpr(:sin, Any[x])
    # This shouldn't happen in regular code, but let's test against it to check
    # we get something similar to AffExpr and QuadExpr.
    empty!(f.args)
    @test owner_model(f) === nothing
    return
end

function test_operate_shortcut_ma_operate!!_add_mul()
    model = Model()
    @variable(model, x)
    @test isequal_canonical(
        @expression(model, sum(sin(x) for i in 1:3)),
        NonlinearExpr(:+, Any[sin(x), sin(x), sin(x)]),
    )
    @test isequal_canonical(JuMP._MA.add_mul(sin(x), 2, x), sin(x) + 2x)
    @test isequal_canonical(JuMP._MA.add_mul(sin(x), 2, x, 2), sin(x) + 4x)
    @test isequal_canonical(JuMP._MA.sub_mul(sin(x), 2, x), sin(x) - 2x)
    @test isequal_canonical(JuMP._MA.sub_mul(sin(x), 2, x, 2), sin(x) - 4x)
    return
end

function test_show_nonlinear_model()
    model = Model()
    @variable(model, x >= -1)
    @objective(model, Min, exp(x))
    @constraint(model, sin(x) <= 0)
    str = sprint(show, model)
    @test occursin("NonlinearExpr", str)
    return
end

function test_error_both_nl_interfaces_constraint()
    model = Model()
    @variable(model, x)
    @constraint(model, log(x) <= 1)
    @NLconstraint(model, log(x) <= 1)
    @test_throws(
        ErrorException(
            "Cannot optimize a model which contains the features from " *
            "both the legacy (macros beginning with `@NL`) and new " *
            "(`NonlinearExpr`) nonlinear interfaces. You must use one or " *
            "the other.",
        ),
        optimize!(model),
    )
    return
end

function test_error_both_nl_interfaces_objective()
    model = Model()
    @variable(model, x)
    @objective(model, Max, log(x))
    @NLconstraint(model, log(x) <= 1)
    @test_throws(
        ErrorException(
            "Cannot optimize a model which contains the features from " *
            "both the legacy (macros beginning with `@NL`) and new " *
            "(`NonlinearExpr`) nonlinear interfaces. You must use one or " *
            "the other.",
        ),
        optimize!(model),
    )
    return
end

function test_VectorNonlinearFunction_moi_function()
    model = Model()
    @variable(model, x)
    F = [sin(x)]
    @test moi_function_type(typeof(F)) == MOI.VectorNonlinearFunction
    @test isapprox(
        moi_function(F),
        MOI.VectorNonlinearFunction([
            MOI.ScalarNonlinearFunction(:sin, Any[index(x)]),
        ]),
    )
    @test MOI.VectorNonlinearFunction(F) ≈ moi_function(F)
    @test jump_function_type(model, MOI.VectorNonlinearFunction) ==
          Vector{NonlinearExpr}
    @test isequal_canonical(jump_function(model, moi_function(F)), F)
    return
end

function test_VectorNonlinearFunction_moi_function_conversion()
    model = Model()
    @variable(model, x)
    F = [sin(x), x, x + 1, x^2]
    @test F isa Vector{NonlinearExpr}
    @test moi_function_type(typeof(F)) == MOI.VectorNonlinearFunction
    @test isapprox(
        moi_function(F),
        MOI.VectorNonlinearFunction([
            MOI.ScalarNonlinearFunction(:sin, Any[index(x)]),
            MOI.ScalarNonlinearFunction(:+, Any[index(x)]),
            MOI.ScalarNonlinearFunction(:+, Any[index(x), 1.0]),
            MOI.ScalarNonlinearFunction(:*, Any[index(x), index(x)]),
        ]),
    )
    @test MOI.VectorNonlinearFunction(F) ≈ moi_function(F)
    @test jump_function_type(model, MOI.VectorNonlinearFunction) ==
          Vector{NonlinearExpr}
    @test isequal_canonical(jump_function(model, moi_function(F)), F)
    return
end

function test_VectorNonlinearFunction_moi_function_conversion_variable()
    model = Model()
    @variable(model, x)
    F = [sin(x), x]
    @test F isa Vector{NonlinearExpr}
    @test moi_function_type(typeof(F)) == MOI.VectorNonlinearFunction
    @test isapprox(
        moi_function(F),
        MOI.VectorNonlinearFunction([
            MOI.ScalarNonlinearFunction(:sin, Any[index(x)]),
            MOI.ScalarNonlinearFunction(:+, Any[index(x)]),
        ]),
    )
    @test MOI.VectorNonlinearFunction(F) ≈ moi_function(F)
    @test jump_function_type(model, MOI.VectorNonlinearFunction) ==
          Vector{NonlinearExpr}
    @test isequal_canonical(jump_function(model, moi_function(F)), F)
    return
end

function test_VectorNonlinearFunction_objective()
    model = Model()
    @variable(model, x)
    F = [sin(x), sqrt(x)]
    @objective(model, Min, F)
    @test objective_function_type(model) == Vector{NonlinearExpr}
    @test isequal_canonical(objective_function(model), F)
    return
end

function test_operator_overload_complex_error()
    model = Model()
    @variable(model, x)
    f = (1 + 2im) * x
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($f)::$(typeof(f))`",
        ),
        sin(f),
    )
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($(1 + 2im))::$(typeof(1 + 2im))`",
        ),
        +(sin(x), 1 + 2im),
    )
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($(1 + 2im))::$(typeof(1 + 2im))`",
        ),
        +(1 + 2im, sin(x)),
    )
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($f)::$(typeof(f))`",
        ),
        +(f, sin(x)),
    )
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($f)::$(typeof(f))`",
        ),
        +(sin(x), f),
    )
    return
end

function test_redefinition_of_function()
    model = Model()
    f(x) = x^2
    err = try
        JuMP._catch_redefinition_constant_error(:f, f)
    catch err
        err
    end
    @test_throws(err, @operator(model, f, 1, f))
    return
end

function test_moi_function_abstract_jump_scalar()
    model = Model()
    @variable(model, x)
    y = AbstractJuMPScalar[x, sin(x)]
    @test_throws ErrorException moi_function(y)
    return
end

function test_linear_algebra_errors()
    model = Model()
    @variable(model, x[1:2, 1:2])
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.det(x)
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.logdet(x)
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.norm(x)
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.nullspace(x)
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.qr(x)
    y = 2.0 .* x[:, 2] .+ 1.0
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.norm(y)
    @test_throws MOI.UnsupportedNonlinearOperator LinearAlgebra.nullspace(y)
    return
end

function test_ma_zero_in_operate!!()
    model = Model()
    @variable(model, x)
    y = @expression(model, sum(sin(x) for i in 1:2) + sum(1 for i in 1:0))
    @test isequal_canonical(y, sin(x) + sin(x))
    a = sin(x) + sin(x)
    y = MA.operate!!(MA.add_mul, a, MA.Zero())
    @test y === a
    @test isequal_canonical(y, sin(x) + sin(x))
    return
end

function test_ma_operate!!_nested_sum()
    model = Model()
    @variable(model, x)
    y = NonlinearExpr(:+, Any[x])
    z = MA.operate!!(MA.add_mul, y, y)
    @test isequal_canonical(y, @force_nonlinear(+x))
    @test isequal_canonical(z, @force_nonlinear(+(+x, +x)))
    return
end

function test_nonlinear_operator_inferred()
    model = Model()
    @variable(model, x)
    @inferred op_less_than_or_equal_to(x, 1)
    @test @inferred(op_less_than_or_equal_to(1, 2)) == true
    return
end

function test_generic_nonlinear_expr_infer_variable_type()
    model = Model()
    @variable(model, x)
    @inferred GenericNonlinearExpr(:sin, x)
    @inferred GenericNonlinearExpr GenericNonlinearExpr(:sin, Any[x])
    f = sin(x)
    @test isequal_canonical(GenericNonlinearExpr(:sin, x), f)
    @test isequal_canonical(GenericNonlinearExpr(:sin, Any[x]), f)
    g = @expression(model, 1 <= x)
    @inferred GenericNonlinearExpr(:<=, 1, x)
    @inferred GenericNonlinearExpr GenericNonlinearExpr(:<=, Any[1, x])
    @test isequal_canonical(GenericNonlinearExpr(:<=, 1, x), g)
    @test isequal_canonical(GenericNonlinearExpr(:<=, Any[1, x]), g)
    @test_throws(
        ErrorException(
            "Unable to create a nonlinear expression because it did not " *
            "contain any JuMP scalars. head = `:sin`, args = `(1,)`.",
        ),
        GenericNonlinearExpr(:sin, 1),
    )
    @test_throws(
        ErrorException(
            "Unable to create a nonlinear expression because it did not " *
            "contain any JuMP scalars. head = `:sin`, args = `Any[1]`.",
        ),
        GenericNonlinearExpr(:sin, Any[1]),
    )
    return
end

function test_add_to_expression!()
    model = Model()
    @variable(model, x)
    y = zero(NonlinearExpr)
    @test_throws(
        ErrorException(
            """
            `add_to_expression!` is not supported for expressions of type
            `$(typeof(y))` because they cannot be modified in-place.
            Instead of `add_to_expression!(expr, args..)`, use one of the following:
            ```julia
            expr += *(args...)
            # or
            import MutableArithmetics as MA
            expr = MA.add_mul!!(expr, args...)
            ```
            """,
        ),
        add_to_expression!(y, 2.0, sin(x)),
    )
    return
end

function test_operator_min()
    model = Model()
    @variable(model, x)
    @test isequal_canonical(min(x, 1), NonlinearExpr(:min, Any[x, 1.0]))
    @test isequal_canonical(min(1, x, x^2), min(min(1.0, x), x^2))
    for f in (x, 1.0 * x + 2.0, x^2, sin(x))
        @test min(f) === f
    end
    return
end

function test_operator_max()
    model = Model()
    @variable(model, x)
    @test isequal_canonical(max(x, 1), NonlinearExpr(:max, Any[x, 1.0]))
    @test isequal_canonical(max(1, x, x^2), max(max(1.0, x), x^2))
    for f in (x, 1.0 * x + 2.0, x^2, sin(x))
        @test max(f) === f
    end
    return
end

function test_variable_ref_type()
    for V in (GenericVariableRef{Int}, VariableRef)
        @test variable_ref_type(GenericNonlinearExpr{V}) == V
    end
    return
end

function test_printing_truncation()
    model = Model()
    @variable(model, x[1:100])
    y = @expression(model, sum(sin.(x) .* 2))
    @test occursin(
        "(sin(x[72]) * 2.0) + [[...41 terms omitted...]] + (sin(x[30]) * 2.0)",
        function_string(MIME("text/plain"), y),
    )
    @test occursin(
        "{\\left({\\textsf{sin}\\left({x_{72}}\\right)} * {2.0}\\right) + {[[\\ldots\\text{41 terms omitted}\\ldots]]} + {\\left({\\textsf{sin}\\left({x_{30}}\\right)} * {2.0}\\right)}",
        function_string(MIME("text/latex"), y),
    )
    return
end

function test_convert_vector_aff_expr()
    model = Model()
    @variable(model, x)
    @test [sin(x), x] isa Vector{NonlinearExpr}
    @test [sin(x), x + 1] isa Vector{NonlinearExpr}
    @test [sin(x), convert(AffExpr, x)] isa Vector{NonlinearExpr}
    return
end

function test_convert_float_nonlinear_expr()
    model = Model()
    @variable(model, x)
    @test [0.0, x, sin(x)] isa Vector{NonlinearExpr}
    @test [0.0, sin(x), x] isa Vector{NonlinearExpr}
    @test [x, 0.0, sin(x)] isa Vector{NonlinearExpr}
    @test [x, sin(x), 0.0] isa Vector{NonlinearExpr}
    @test [sin(x), 0.0, x] isa Vector{NonlinearExpr}
    @test [sin(x), x, 0.0] isa Vector{NonlinearExpr}
    @test isequal_canonical(
        convert(NonlinearExpr, 2),
        NonlinearExpr(:+, Any[2]),
    )
    return
end

function test_nlp_adjoint()
    model = Model()
    @variable(model, x)
    y = sin(x)
    @test y' === y
    @test conj(y) === y
    @test real(y) === y
    @test isequal_canonical(imag(y), zero(y))
    @test isequal_canonical(abs2(y), y^2)
    @test isreal(y)
    return
end

function test_nlp_matrix_adjoint()
    model = Model()
    @variable(model, x[1:2])
    y = sin.(x)
    @test isequal_canonical(
        @expression(model, expr, y' * y),
        NonlinearExpr(:+, Any[0.0, y[2]*y[2], y[1]*y[1]]),
    )
    return
end

function test_error_complex_literal_pow()
    model = Model()
    @variable(model, x in ComplexPlane())
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($x)::$(typeof(x))`",
        ),
        x^3,
    )
    return
end

function test_error_nonlinear_expr_complex_constructor()
    model = Model()
    @variable(model, x in ComplexPlane())
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($x)::$(typeof(x))`",
        ),
        NonlinearExpr(:^, Any[x, 3]),
    )
    @test_throws(
        ErrorException(
            "Cannot build `GenericNonlinearExpr` because a term is complex-" *
            "valued: `($x)::$(typeof(x))`",
        ),
        NonlinearExpr(:^, x, 3),
    )
    return
end

function test_error_legacy_expression_constructor()
    model = Model()
    @variable(model, x)
    @NLexpression(model, arg, x^3)
    err = ErrorException(
        """
        Cannot mix a legacy NonlinearExpression with the new nonlinear API.

        Got: $arg

        To update, replace all calls to `@NLexpression` with `@expression`.
        """,
    )
    @test_throws err @objective(model, Min, arg)
    @test_throws err @constraint(model, arg <= 0)
    @test_throws err @constraint(model, arg in MOI.LessThan(0.0))
    @test_throws err moi_function(arg)
    @test_throws err x * arg
    @test_throws err arg * x
    return
end

function test_error_legacy_parameter_constructor()
    model = Model()
    @variable(model, x)
    @NLparameter(model, p == 1)
    err = ErrorException(
        """
        Cannot mix a legacy NonlinearParameter with the new nonlinear API.

        Got: $p

        To update, replace calls to:
        ```julia
        @NLparameter(model, p == 1)
        ```
        with
        ```julia
        @variable(model, p in Parameter(1))
        ```
        """,
    )
    @test_throws err @objective(model, Min, p)
    @test_throws err @constraint(model, p <= 0)
    @test_throws err @constraint(model, p in MOI.LessThan(0.0))
    @test_throws err moi_function(p)
    @test_throws err x * p
    @test_throws err p * x
    return
end

end  # module
