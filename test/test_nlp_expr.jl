#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module TestNLPExpr

using JuMP
using Test

function test_univariate_operators()
    model = Model()
    @variable(model, x)
    for f in MOI.Nonlinear.DEFAULT_UNIVARIATE_OPERATORS
        if f in (:+, :-, :abs2)
            op = getfield(Base, f)
            @test op(sin(x)) isa NonlinearExpr
        elseif isdefined(Base, f)
            op = getfield(Base, f)
            @test op(x) isa NonlinearExpr
        elseif isdefined(MOI.Nonlinear.SpecialFunctions, f)
            op = getfield(MOI.Nonlinear.SpecialFunctions, f)
            @test op(x) isa NonlinearExpr
        end
    end
    return
end

function test_binary_operators()
    model = Model()
    @variable(model, x)
    num, aff, quad, nlp = 1.0, 1.0 + x, x^2, sin(x)
    for op in (+, -, *, /), a in (num, x, aff, quad, nlp)
        @test op(a, nlp) isa NonlinearExpr
        @test op(nlp, a) isa NonlinearExpr
    end
    for op in (*, /), a in (x, aff)
        @test op(a, quad) isa NonlinearExpr
        @test op(quad, a) isa NonlinearExpr
    end
    for a in (num, x, aff, quad), b in (x, aff, quad)
        @test /(a, b) isa NonlinearExpr
    end
    return
end

function test_objective()
    model = Model()
    @variable(model, x)
    @objective(model, Min, 2.0 * sin(x)^2 + cos(x) / x)
    @test objective_function(model) isa NonlinearExpr
    return
end

function test_expression()
    model = Model()
    @variable(model, x)
    @variable(model, y[1:3])
    @test string(@expression(model, *(y...))) == "*(y[1]*y[2], y[3])"
    @test string(@expression(model, sin(x))) == "sin(x)"
    @test string(@expression(model, 2^x)) == "^(2.0, x)"
    @test string(@expression(model, x^x)) == "^(x, x)"
    @test string(@expression(model, sin(x)^2)) == "^(sin(x), 2.0)"
    @test string(@expression(model, sin(x)^2.0)) == "^(sin(x), 2.0)"
    @test string(@expression(model, 2 * sin(x)^2.0)) == "*(2.0, ^(sin(x), 2.0))"
    @test string(@expression(model, 1 + sin(x))) == "+(1.0, sin(x))"
    @test string(@expression(model, 1 + 2 * sin(x))) == "+(1.0, *(2.0, sin(x)))"
    @test string(@expression(model, 2.0 * sin(x)^2 + cos(x) / x)) ==
          "+(*(2.0, ^(sin(x), 2.0)), /(cos(x), x))"
    @test string(@expression(model, 2.0 * sin(x)^2 - cos(x) / x)) ==
          "-(*(2.0, ^(sin(x), 2.0)), /(cos(x), x))"
    return
end

function test_flatten_nary()
    model = Model()
    @variable(model, x)
    @test string(zero(NonlinearExpr) + 1) == "+(+(0.0), 1.0)"
    @test string(zero(NonlinearExpr) + x) == "+(+(0.0), x)"
    @test string(sin(x) + sin(x) + 1) == "+(+(sin(x), sin(x)), 1.0)"
    @test string(sin(x) + sin(x) + x) == "+(+(sin(x), sin(x)), x)"
    @test string(sin(x) * sin(x) * 1) == "*(*(sin(x), sin(x)), 1.0)"
    @test string(sin(x) * sin(x) * x) == "*(*(sin(x), sin(x)), x)"
    return
end

function test_zero_one()
    @test string(zero(NonlinearExpr)) == "+(0.0)"
    @test string(one(NonlinearExpr)) == "+(1.0)"
    return
end

function test_latex()
    model = Model()
    @variable(model, x)
    @test function_string(MIME("text/latex"), sin(x)) == "\\textsf{sin(x)}"
    @test function_string(MIME("text/plain"), sin(x)) == "sin(x)"
    return
end

function test_expression_addmul()
    model = Model()
    @variable(model, x)
    @test string(@expression(model, x + 3 * sin(x))) == "+(x, *(3.0, sin(x)))"
    @test string(@expression(model, 2 * x + 3 * sin(x))) ==
          "+(2 x, *(3.0, sin(x)))"
    @test string(@expression(model, x^2 + 3 * sin(x))) ==
          "+($(x^2), *(3.0, sin(x)))"
    @test string(@expression(model, sin(x) + 3 * sin(x))) ==
          "+(sin(x), *(3.0, sin(x)))"
    @test string(@expression(model, sin(x) + 3 * x)) == "+(sin(x), 3 x)"
    @test string(@expression(model, sin(x) + 3 * x * x)) ==
          "+(sin(x), 3 $(x^2))"
    return
end

function test_expression_submul()
    model = Model()
    @variable(model, x)
    @test string(@expression(model, x - 3 * sin(x))) == "-(x, *(3.0, sin(x)))"
    @test string(@expression(model, 2 * x - 3 * sin(x))) ==
          "-(2 x, *(3.0, sin(x)))"
    @test string(@expression(model, x^2 - 3 * sin(x))) ==
          "-($(x^2), *(3.0, sin(x)))"
    @test string(@expression(model, sin(x) - 3 * sin(x))) ==
          "-(sin(x), *(3.0, sin(x)))"
    @test string(@expression(model, sin(x) - 3 * x)) == "-(sin(x), 3 x)"
    @test string(@expression(model, sin(x) - 3 * x * x)) ==
          "-(sin(x), 3 $(x^2))"
    return
end

function test_aff_expr_convert()
    model = Model()
    @variable(model, x)
    _to_string(x) = string(convert(NonlinearExpr, x))
    @test _to_string(AffExpr(0.0)) == "0.0"
    @test _to_string(AffExpr(1.0)) == "1.0"
    @test _to_string(x + 1) == "+(x, 1.0)"
    @test _to_string(2x + 1) == "+(*(2.0, x), 1.0)"
    @test _to_string(2x) == "*(2.0, x)"
    return
end

function test_quad_expr_convert()
    model = Model()
    @variable(model, x)
    _to_string(x) = string(convert(NonlinearExpr, x))
    @test _to_string(QuadExpr(AffExpr(0.0))) == "0.0"
    @test _to_string(QuadExpr(AffExpr(1.0))) == "1.0"
    @test _to_string(x^2 + 1) == "+(*(x, x), 1.0)"
    @test _to_string(2x^2 + 1) == "+(*(2.0, x, x), 1.0)"
    @test _to_string(2x^2) == "*(2.0, x, x)"
    @test _to_string(x^2 + x + 1) == "+(x, *(x, x), 1.0)"
    @test _to_string(2x^2 + x + 1) == "+(x, *(2.0, x, x), 1.0)"
    @test _to_string(2x^2 + x) == "+(x, *(2.0, x, x))"
    @test _to_string(x^2 + 2x + 1) == "+(*(2.0, x), *(x, x), 1.0)"
    @test _to_string(2x^2 + 2x + 1) == "+(*(2.0, x), *(2.0, x, x), 1.0)"
    @test _to_string(2x^2 + 2x) == "+(*(2.0, x), *(2.0, x, x))"
    return
end

function test_constraint_name()
    model = Model()
    @variable(model, x)
    @constraint(model, c, sin(x) <= 1)
    @test name(c) == "c"
    set_name(c, "d")
    @test name(c) == "d"
    @test startswith(string(c), "d : ")
    return
end

function test_constraint_lessthan()
    model = Model()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x <= 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.LessThan(0.0)
    return
end

function test_constraint_greaterthan()
    model = Model()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x >= 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.GreaterThan(0.0)
    return
end

function test_constraint_equalto()
    model = Model()
    @variable(model, x)
    @constraint(model, c, 2.0 * sin(x)^2 + cos(x) / x == 1)
    obj = constraint_object(c)
    @test isequal_canonical(obj.func, 2.0 * sin(x)^2 + cos(x) / x - 1)
    @test obj.set == MOI.EqualTo(0.0)
    return
end

function test_constraint_interval()
    model = Model()
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
    f(x::AbstractJuMPScalar) = NonlinearExpr(:f, x)
    register(model, :f, 1, f; autodiff = true)
    @test string(@expression(model, f(x))) == "f(x)"
    @test string(f(x) + f(x)) == "+(f(x), f(x))"
    return
end

function test_nonlinear_matrix_algebra()
    model = Model()
    @variable(model, X[1:3, 1:3], Symmetric)
    @objective(model, Max, sum(X^4 .- X^3))
    @test objective_function(model) isa NonlinearExpr
    return
end

"""
This test checks that we can work with expressions of arbitrary depth. Don't use
recursion!
"""
function test_recursion_stackoverflow()
    model = Model()
    @variable(model, x)
    expr = sin(x)
    for _ in 1:20_000
        expr = sin(expr)
    end
    @test @objective(model, Min, expr) isa NonlinearExpr
    @test string(expr) isa String
    return
end

function test_nlparameter_interaction()
    model = Model()
    @variable(model, x)
    @NLparameter(model, p == 1)
    e = x + p
    @test e isa NonlinearExpr
    @test string(e) == "+(x, $p)"
    return
end

function test_nlexpression_interaction()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sin(x))
    e = x + expr
    @test e isa NonlinearExpr
    @test string(e) == "+(x, $expr)"
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
    @test string(jump_function(model, nlp[index(expr1)])) == "sin(+($p, $x))"
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

function test_expr_mle()
    data = [1.0, 2.0, 4.0, 8.0]
    n = length(data)
    model = Model()
    @variable(model, x)
    @variable(model, y)
    obj = @expression(
        model,
        n / 2 * log(1 / (2 * y^2)) -
        sum((data[i] - x)^2 for i in 1:n) / (2 * y^2)
    )
    @test string(obj) ==
          "-(*(2.0, log(/(1.0, 2 $(y^2)))), /(4 $(x^2) - 30 x + 85, 2 $(y^2)))"
    return
end

function test_nl_macro()
    model = Model()
    @variable(model, x)
    @test isequal_canonical(@NL(x < 0), NonlinearExpr(:<, Any[x, 0]))
    @test isequal_canonical(@NL(0 >= x), NonlinearExpr(:>=, Any[0, x]))
    @test isequal_canonical(
        @NL(0 < x <= 1),
        NonlinearExpr(:&&, Any[@NL(0 < x), @NL(x <= 1)]),
    )
    @test isequal_canonical(@NL(x || 1), NonlinearExpr(:||, Any[x, 1]))
    @test isequal_canonical(@NL(x && 1), NonlinearExpr(:&&, Any[x, 1]))
    @test isequal_canonical(
        @NL(ifelse(x, 1, 2)),
        NonlinearExpr(:ifelse, Any[x, 1, 2]),
    )
    expr = :(x[])
    @test_throws(
        ErrorException(
            "Unable to convert expression to `NonlinearExpr`: $expr\n\nThe " *
            "`@NL` macro must be used on a single function call. If the " *
            "expression above is not what you intended, try using " *
            "parentheses to disambiguate the parsing of the expression.",
        ),
        JuMP._to_nl(expr)
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

end  # module
