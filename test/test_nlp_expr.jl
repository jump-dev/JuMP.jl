#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestNLPExpr

using JuMP
using Test

function test_extension_univariate_operators(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
        if f in (:+, :-, :abs2)
            op = getfield(Base, f)
            @test op(sin(x)) isa NonlinearExpr{VariableRefType}
        elseif isdefined(Base, f)
            op = getfield(Base, f)
            @test op(x) isa NonlinearExpr{VariableRefType}
        elseif isdefined(MOI.Nonlinear.SpecialFunctions, f)
            op = getfield(MOI.Nonlinear.SpecialFunctions, f)
            @test op(x) isa NonlinearExpr{VariableRefType}
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
        @test op(a, nlp) isa NonlinearExpr{VariableRefType}
        @test op(nlp, a) isa NonlinearExpr{VariableRefType}
    end
    for op in (*, /), a in (x, aff)
        @test op(a, quad) isa NonlinearExpr{VariableRefType}
        @test op(quad, a) isa NonlinearExpr{VariableRefType}
    end
    for a in (num, x, aff, quad), b in (x, aff, quad)
        @test /(a, b) isa NonlinearExpr{VariableRefType}
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
    @test objective_function(model) isa NonlinearExpr{VariableRefType}
    return
end

function test_extension_expression(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @variable(model, y[1:3])
    @test string(@expression(model, *(y...))) == "(y[1]*y[2] * y[3])"
    @test string(@expression(model, sin(x))) == "sin(x)"
    @test string(@expression(model, 2^x)) == "(2.0 ^ x)"
    @test string(@expression(model, x^x)) == "(x ^ x)"
    @test string(@expression(model, sin(x)^2)) == "(sin(x) ^ 2.0)"
    @test string(@expression(model, sin(x)^2.0)) == "(sin(x) ^ 2.0)"
    @test string(@expression(model, 2 * sin(x)^2.0)) == "(2.0 * (sin(x) ^ 2.0))"
    @test string(@expression(model, 1 + sin(x))) == "(1.0 + sin(x))"
    @test string(@expression(model, 1 + 2 * sin(x))) == "(1.0 + (2.0 * sin(x)))"
    @test string(@expression(model, 2.0 * sin(x)^2 + cos(x) / x)) ==
          "((2.0 * (sin(x) ^ 2.0)) + (cos(x) / x))"
    @test string(@expression(model, 2.0 * sin(x)^2 - cos(x) / x)) ==
          "((2.0 * (sin(x) ^ 2.0)) - (cos(x) / x))"
    return
end

function test_extension_flatten_nary(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    expr_plus = NonlinearExpr{VariableRefType}(:+, Any[x])
    expr_mult = NonlinearExpr{VariableRefType}(:*, Any[x])
    expr_sin = NonlinearExpr{VariableRefType}(:sin, Any[x])
    @test string(+(expr_plus, 1)) == "(x + 1.0)"
    @test string(+(1, expr_plus)) == "(1.0 + x)"
    @test string(+(expr_plus, x)) == "(x + x)"
    @test string(+(expr_sin, x)) == "(sin(x) + x)"
    @test string(+(x, expr_plus)) == "(x + x)"
    @test string(+(x, expr_sin)) == "(x + sin(x))"
    @test string(+(expr_plus, expr_plus)) == "(x + x)"
    @test string(+(expr_plus, expr_sin)) == "(x + sin(x))"
    @test string(+(expr_sin, expr_plus)) == "(sin(x) + x)"
    @test string(+(expr_sin, expr_sin)) == "(sin(x) + sin(x))"
    @test string(*(expr_mult, 2)) == "(x * 2.0)"
    @test string(*(2, expr_mult)) == "(2.0 * x)"
    @test string(*(expr_mult, x)) == "(x * x)"
    @test string(*(expr_sin, x)) == "(sin(x) * x)"
    @test string(*(x, expr_mult)) == "(x * x)"
    @test string(*(x, expr_sin)) == "(x * sin(x))"
    @test string(*(expr_mult, expr_mult)) == "(x * x)"
    @test string(*(expr_mult, expr_sin)) == "(x * sin(x))"
    @test string(*(expr_sin, expr_mult)) == "(sin(x) * x)"
    @test string(*(expr_sin, expr_sin)) == "(sin(x) * sin(x))"
    return
end

function test_extension_error_associative(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    if VariableRefType == VariableRef
        return
    end
    model1 = ModelType()
    @variable(model1, x1)
    model2 = Model()
    @variable(model2, x2)
    @test_throws(
        ErrorException(
            "Unable to call + with nonlinear expressions of different " *
            "variable type",
        ),
        +(sin(x1), sin(x2)),
    )
    return
end

function test_extension_zero_one(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    @test string(zero(NonlinearExpr{VariableRefType})) == "+(0.0)"
    @test string(one(NonlinearExpr{VariableRefType})) == "+(1.0)"
    return
end

function test_extension_latex(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x)
    @test function_string(MIME("text/latex"), sin(x)) ==
          raw"\textsf{sin}\left({x}\right)"
    @test function_string(MIME("text/plain"), sin(x)) == "sin(x)"
    @expression(model, g, ifelse(x > 0, sin(x), x + cos(x)^2))
    @test function_string(MIME("text/latex"), g) ==
          raw"\textsf{ifelse}\left({\left({x} > {0}\right)}, {\textsf{sin}\left({x}\right)}, {\left({x} + {\left({\textsf{cos}\left({x}\right)} ^ {2.0}\right)}\right)}\right)"
    return
end

function test_extension_expression_addmul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test string(@expression(model, x + 3 * sin(x))) == "(x + (3.0 * sin(x)))"
    @test string(@expression(model, 2 * x + 3 * sin(x))) ==
          "(2 x + (3.0 * sin(x)))"
    @test string(@expression(model, x^2 + 3 * sin(x))) ==
          "($(x^2) + (3.0 * sin(x)))"
    @test string(@expression(model, sin(x) + 3 * sin(x))) ==
          "(sin(x) + (3.0 * sin(x)))"
    @test string(@expression(model, sin(x) + 3 * x)) == "(sin(x) + 3 x)"
    @test string(@expression(model, sin(x) + 3 * x * x)) ==
          "(sin(x) + 3 $(x^2))"
    return
end

function test_extension_expression_submul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test string(@expression(model, x - 3 * sin(x))) == "(x - (3.0 * sin(x)))"
    @test string(@expression(model, 2 * x - 3 * sin(x))) ==
          "(2 x - (3.0 * sin(x)))"
    @test string(@expression(model, x^2 - 3 * sin(x))) ==
          "($(x^2) - (3.0 * sin(x)))"
    @test string(@expression(model, sin(x) - 3 * sin(x))) ==
          "(sin(x) - (3.0 * sin(x)))"
    @test string(@expression(model, sin(x) - 3 * x)) == "(sin(x) - 3 x)"
    @test string(@expression(model, sin(x) - 3 * x * x)) ==
          "(sin(x) - 3 $(x^2))"
    return
end

function test_extension_aff_expr_convert(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    _to_string(x) = string(convert(NonlinearExpr{VariableRefType}, x))
    @test _to_string(AffExpr(0.0)) == "0.0"
    @test _to_string(AffExpr(1.0)) == "1.0"
    @test _to_string(x + 1) == "(x + 1.0)"
    @test _to_string(2x + 1) == "((2.0 * x) + 1.0)"
    @test _to_string(2x) == "(2.0 * x)"
    return
end

function test_extension_quad_expr_convert(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    _to_string(x) = string(convert(NonlinearExpr{VariableRefType}, x))
    @test _to_string(QuadExpr(AffExpr(0.0))) == "0.0"
    @test _to_string(QuadExpr(AffExpr(1.0))) == "1.0"
    @test _to_string(x^2 + 1) == "((x * x) + 1.0)"
    @test _to_string(2x^2 + 1) == "((2.0 * x * x) + 1.0)"
    @test _to_string(2x^2) == "(2.0 * x * x)"
    @test _to_string(x^2 + x + 1) == "(x + (x * x) + 1.0)"
    @test _to_string(2x^2 + x + 1) == "(x + (2.0 * x * x) + 1.0)"
    @test _to_string(2x^2 + x) == "(x + (2.0 * x * x))"
    @test _to_string(x^2 + 2x + 1) == "((2.0 * x) + (x * x) + 1.0)"
    @test _to_string(2x^2 + 2x + 1) == "((2.0 * x) + (2.0 * x * x) + 1.0)"
    @test _to_string(2x^2 + 2x) == "((2.0 * x) + (2.0 * x * x))"
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
    @test obj.set == MOI.LessThan(0.0)
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
    @test obj.set == MOI.GreaterThan(0.0)
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
    @test obj.set == MOI.EqualTo(0.0)
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
    @test obj.set == MOI.Interval(0.0, 1.0)
    return
end

function test_user_defined_function_overload()
    model = Model()
    @variable(model, x)
    f(x::Real) = x^2
    f(x::AbstractJuMPScalar) = NonlinearExpr{VariableRef}(:f, x)
    register(model, :f, 1, f; autodiff = true)
    @test string(@expression(model, f(x))) == "f(x)"
    @test string(f(x) + f(x)) == "(f(x) + f(x))"
    return
end

function test_extension_nonlinear_matrix_algebra(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, X[1:3, 1:3], Symmetric)
    @objective(model, Max, sum(X^4 .- X^3))
    @test objective_function(model) isa NonlinearExpr{VariableRefType}
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
    @test @objective(model, Min, expr) isa NonlinearExpr{VariableRefType}
    @test string(expr) isa String
    return
end

function test_nlparameter_interaction()
    model = Model()
    @variable(model, x)
    @NLparameter(model, p == 1)
    e = x + p
    @test e isa NonlinearExpr
    @test string(e) == "(x + $p)"
    return
end

function test_nlexpression_interaction()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sin(x))
    e = x + expr
    @test e isa NonlinearExpr
    @test string(e) == "(x + $expr)"
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

function test_jump_function_nonlinearexpr()
    model = Model()
    @variable(model, x)
    @NLparameter(model, p == 1)
    @NLexpression(model, expr1, sin(p + x))
    @NLexpression(model, expr2, sin(expr1))
    nlp = nonlinear_model(model)
    @test string(jump_function(model, nlp[index(expr1)])) == "sin(($p + $x))"
    @test string(jump_function(model, nlp[index(expr2)])) == "sin($expr1)"
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
          "((2.0 * log((1.0 / 2 $(y^2)))) - (4 $(x^2) - 30 x + 85 / 2 $(y^2)))"
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
        NonlinearExpr(:ifelse, Any[x, 1, 2]),
    )
    @test isequal_canonical(
        @expression(model, x || 1),
        NonlinearExpr(:||, Any[x, 1]),
    )
    @test isequal_canonical(
        @expression(model, x && 1),
        NonlinearExpr(:&&, Any[x, 1]),
    )
    @test isequal_canonical(
        @expression(model, x < 0),
        NonlinearExpr(:<, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x > 0),
        NonlinearExpr(:>, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x <= 0),
        NonlinearExpr(:<=, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x >= 0),
        NonlinearExpr(:>=, Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, x == 0),
        NonlinearExpr(:(==), Any[x, 0]),
    )
    @test isequal_canonical(
        @expression(model, 0 < x <= 1),
        NonlinearExpr(
            :&&,
            Any[@expression(model, 0 < x), @expression(model, x <= 1)],
        ),
    )
    @test isequal_canonical(
        @expression(model, ifelse(x > 0, x^2, sin(x))),
        NonlinearExpr(:ifelse, Any[@expression(model, x > 0), x^2, sin(x)]),
    )
    return
end

function test_register_univariate()
    model = Model()
    @variable(model, x)
    @register(model, f, 1, x -> x^2)
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_univariate_gradient()
    model = Model()
    @variable(model, x)
    @register(model, f, 1, x -> x^2, x -> 2 * x)
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_univariate_gradient_hessian()
    model = Model()
    @variable(model, x)
    @register(model, f, 1, x -> x^2, x -> 2 * x, x -> 2.0)
    @test isequal_canonical(@expression(model, f(x)), f(x))
    @test isequal_canonical(f(x), NonlinearExpr(:f, Any[x]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:f, 1) in attrs
    return
end

function test_register_multivariate_()
    model = Model()
    @variable(model, x[1:2])
    f = (x...) -> sum(x .^ 2)
    @register(model, foo, 2, f)
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
    @register(model, foo, 2, f, ∇f)
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
    @register(model, foo, 2, f, ∇f, ∇²f)
    @test isequal_canonical(@expression(model, foo(x...)), foo(x...))
    @test isequal_canonical(foo(x...), NonlinearExpr(:foo, Any[x...]))
    attrs = MOI.get(model, MOI.ListOfModelAttributesSet())
    @test MOI.UserDefinedFunction(:foo, 2) in attrs
    return
end

function test_register_errors()
    model = Model()
    @test_throws(
        ErrorException(
            "Unable to register user-defined function foo: invalid number of " *
            "functions provided. Got 0, but expected 1 (if function only), " *
            "2 (if function and gradient), or 3 (if function, gradient, and " *
            "hesssian provided)",
        ),
        @register(model, foo, 2),
    )
    return
end

function test_expression_no_variable()
    head, args = :sin, Any[1]
    @test_throws(
        ErrorException(
            "Unable to create a nonlinear expression because it did not " *
            "contain any JuMP scalars. head = $head, args = $args.",
        ),
        NonlinearExpr(head, args),
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
    @register(model, my_foo, 1, foo)
    @register(model, my_bar, 2, bar)
    @test value(f, my_foo(x)) ≈ (1.1 - 1)^2
    @test value(f, my_foo(x + 1)) ≈ (1.1 + 1 - 1)^2
    @test value(f, my_foo(x^2 + 1)) ≈ (1.1^2 + 1 - 1)^2
    @test value(f, my_foo(x^2 + x + 1)) ≈ (1.1^2 + 1.1 + 1 - 1)^2
    y = QuadExpr(x + 1)
    @test value(f, my_foo(y)) ≈ (value(f, y) - 1)^2
    @test value(f, my_bar(2.2, x)) ≈ sqrt(2.2 - 1.1)
    bad_udf = UserDefinedFunction(:bad_udf)
    @test_throws(
        ErrorException(
            "Unable to evaluate nonlinear operator bad_udf because it is not " *
            "registered",
        ),
        value(f, bad_udf(x)),
    )
    return
end

end  # module
