module TestNLP

using JuMP
using Test

import LinearAlgebra
import SparseArrays

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

include(joinpath(@__DIR__, "utilities.jl"))
include(joinpath(@__DIR__, "JuMPExtension.jl"))

function test_univariate_error()
    model = Model()
    @variable(model, x >= 0)
    @test_throws ErrorException @NLobjective(model, Min, g_doesnotexist(x))
end

function test_univariate_error_existing()
    model = Model()
    @variable(model, x >= 0)
    @NLexpression(model, ex, x^2)
    @test_throws ErrorException @NLobjective(model, Min, g_doestnotexist(ex))
end

function test_univariate()
    model = Model()
    @variable(model, x >= 0)
    g(x) = x^2
    @test_logs (:warn,) @NLobjective(model, Min, g(x))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 4.0
end

function test_univariate_register_twice()
    model = Model()
    @variable(model, x >= 0)
    g(x) = x^2
    @test_logs (:warn,) @NLobjective(model, Min, g(x))
    @test_logs @NLconstraint(model, g(x) <= 1)
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    y = [NaN]
    MOI.eval_constraint(d, y, x)
    @test y == [3.0]
end

function test_univariate_register_twice_error()
    model = Model()
    @variable(model, x >= 0)
    g(x) = x^2
    g(x, y) = x^2 + x^2
    @test_logs (:warn,) @NLobjective(model, Min, g(x))
    @test_throws ErrorException @NLconstraint(model, g(x, x) <= 1)
end

function test_univariate_existing_nlpdata()
    model = Model()
    @variable(model, x >= 0)
    @NLexpression(model, ex, x^2)
    g(x) = x^2
    @test_logs (:warn,) @NLobjective(model, Min, g(ex))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 16.0
end

function test_univariate_redefine()
    model = Model()
    @variable(model, x >= 0)
    g = (x) -> x^2
    @test_logs (:warn,) @NLobjective(model, Min, g(x))

    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 4.0

    g = (x) -> 2x^2
    @test MOI.eval_objective(d, x) == 4.0
end

function test_multivariate_error()
    model = Model()
    @variable(model, x >= 0)
    @test_throws ErrorException @NLobjective(model, Min, g_doesnotexist(x, x))
end

function test_multivariate_error_existing()
    model = Model()
    @variable(model, x >= 0)
    @NLexpression(model, ex, x^2)
    @test_throws ErrorException @NLobjective(model, Min, g_doestnotexist(ex, x))
end

function test_multivariate()
    model = Model()
    @variable(model, x >= 0)
    g(x, y) = x^2 + y^2
    @test_logs (:warn,) @NLobjective(model, Min, g(x, x))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 8.0
end

function test_multivariate_register_warn()
    model = Model()
    g(x, y) = x^2 + y^2
    function ∇g(g::Vector{T}, x::T, y::T) where {T<:Real}
        g[1] = y
        g[2] = x
        return
    end
    @test_logs (:warn,) register(model, :g, 2, g, ∇g; autodiff = true)
end

function test_multivariate_existing_nlpdata()
    model = Model()
    @variable(model, x >= 0)
    @NLexpression(model, ex, x^2)
    g(x, y) = x^2 + y^2
    @test_logs (:warn,) @NLobjective(model, Min, g(ex, x))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 20.0
end

function test_multivariate_redefine()
    model = Model()
    @variable(model, x >= 0)
    @NLexpression(model, ex, x^2)
    g = (x, y) -> x^2 + y^2
    @test_logs (:warn,) @NLobjective(model, Min, g(ex, x))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, Symbol[])
    x = [2.0]
    @test MOI.eval_objective(d, x) == 20.0

    g = (x, y) -> x^2 + y
    @test MOI.eval_objective(d, x) == 20.0
end

function test_multivariate_register_splat()
    model = Model()
    @variable(model, x[1:2])
    f(x, y) = x + y
    err = ErrorException("""
        Unrecognized function "f" used in nonlinear expression.

        You must register it as a user-defined function before building
        the model. For example, replacing `N` with the appropriate number
        of arguments, do:
        ```julia
        model = Model()
        register(model, :f, N, f, autodiff=true)
        # ... variables and constraints ...
        ```
        """)
    @test_throws err @NLexpression(model, ex, f(x...))
end

function test_multivariate_register_splat_existing()
    model = Model()
    @variable(model, x[1:2])
    f(x, y) = x + y
    @NLconstraint(model, x[1]^2 <= 1)
    err = ErrorException("""
        Unrecognized function "f" used in nonlinear expression.

        You must register it as a user-defined function before building
        the model. For example, replacing `N` with the appropriate number
        of arguments, do:
        ```julia
        model = Model()
        register(model, :f, N, f, autodiff=true)
        # ... variables and constraints ...
        ```
        """)
    @test_throws err @NLexpression(model, ex, f(x...))
end

function test_multivariate_max()
    m = Model()
    @variable(m, x)
    @NLobjective(m, Min, max(0, x))
    nlp = NLPEvaluator(m)
    MOI.initialize(nlp, [:ExprGraph])
    @test MOI.eval_objective(nlp, [-1.0]) == 0.0
    @test MOI.eval_objective(nlp, [1.0]) == 1.0
    @test MOI.objective_expr(nlp) == :(max(0.0, x[$(index(x))]))
    return
end

function test_multivariate_min()
    m = Model()
    @variable(m, x)
    @NLobjective(m, Max, min(0, x))
    nlp = NLPEvaluator(m)
    MOI.initialize(nlp, [:ExprGraph])
    @test MOI.eval_objective(nlp, [-1.0]) == -1.0
    @test MOI.eval_objective(nlp, [1.0]) == 0.0
    @test MOI.objective_expr(nlp) == :(min(0.0, x[$(index(x))]))
    return
end

function test_register_check_forwarddiff_univariate_f()
    model = Model()
    f(x::Float64) = log(x)
    @test_throws(ErrorException, register(model, :f, 1, f; autodiff = true))
    return
end

function test_register_check_forwarddiff_univariate_gradf()
    model = Model()
    f(x) = log(x)
    # This is a common case, where user's type their arguments
    ∇f(x::Float64) = 1 / x
    @test_throws(ErrorException, register(model, :f, 1, f, ∇f; autodiff = true))
    return
end

function test_register_check_forwarddiff_multivariate()
    model = Model()
    function f(x...)
        # This is a common case, where user's preallocate a Float64 storage.
        y = zeros(length(x))
        for i in 1:length(x)
            y[i] = log(x[i])
        end
        return sum(y)
    end
    @test_throws(ErrorException, register(model, :f, 3, f; autodiff = true))
    return
end

"""
    test_register_check_forwarddiff_multivariate_gradf()

Because we disable Hessians, the functions in the multivariate case do not need
to be differentiable.
"""
function test_register_check_forwarddiff_multivariate_gradf()
    model = Model()
    function f(x...)
        # This is a common case, where user's preallocate a Float64 storage.
        y = zeros(length(x))
        for i in 1:length(x)
            y[i] = log(x[i])
        end
        return sum(y)
    end
    function ∇f(x...)
        # This is a common case, where user's preallocate a Float64 storage.
        y = zeros(length(x))
        for i in 1:length(x)
            y[i] = 1 / x[i]
        end
        return sum(y)
    end
    register(model, :f, 3, f, ∇f)
    return
end

function test_all_nonlinear_constraints()
    model = Model()
    @variable(model, x)
    @NLconstraint(model, c1, x^2 <= 1)
    c2 = @NLconstraint(model, [i = 1:2], x^i >= -1)
    @test all_nonlinear_constraints(model) == [c1; c2]
    return
end

function test_parse_plus_binary()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    ex = @NLexpression(m, x + y)
    @test sprint(show, ex) == "subexpression[1]: x + y"
    return
end

function test_parse_plus_ternary()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    ex = @NLexpression(m, x + y + z)
    @test sprint(show, ex) == "subexpression[1]: x + y + z"
    return
end

function test_parse_mult_binary()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    ex = @NLexpression(m, x * y)
    @test sprint(show, ex) == "subexpression[1]: x * y"
    return
end

function test_parse_mult_ternary()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    ex = @NLexpression(m, x * y * z)
    @test sprint(show, ex) == "subexpression[1]: x * y * z"
    return
end

function test_parse_exp_binary()
    m = Model()
    @variable(m, x)
    ex = @NLexpression(m, x^3)
    @test sprint(show, ex) == "subexpression[1]: x ^ 3.0"
    return
end

function test_parse_sin()
    m = Model()
    @variable(m, x)
    ex = @NLexpression(m, sin(x))
    @test sprint(show, ex) == "subexpression[1]: sin(x)"
    return
end

function test_parse_ifelse()
    m = Model()
    @variable(m, x)
    ex = @NLexpression(m, ifelse(1 == 2 || 3 == 4 && 5 == 6, x, 0.0))
    @test sprint(show, ex) ==
          "subexpression[1]: ifelse(1.0 == 2.0 || 3.0 == 4.0 && 5.0 == 6.0, x, 0.0)"
    return
end

function test_parse_ifelse_comparison()
    m = Model()
    @variable(m, x)
    ex = @NLexpression(m, ifelse(1 <= 2 <= 3, x, 0.0))
    @test sprint(show, ex) ==
          "subexpression[1]: ifelse(1.0 <= 2.0 <= 3.0, x, 0.0)"
    return
end

function test_parse_sum()
    m = Model()
    @variable(m, x[1:2])
    ex = @NLexpression(m, sum(x[i] for i in 1:2))
    @test sprint(show, ex) == "subexpression[1]: x[1] + x[2]"
    return
end

function test_parse_prod()
    m = Model()
    @variable(m, x[1:2])
    ex = @NLexpression(m, prod(x[i] for i in 1:2))
    @test sprint(show, ex) == "subexpression[1]: x[1] * x[2]"
    return
end

function test_parse_subexpressions()
    m = Model()
    @variable(m, x)
    @NLexpression(m, ex, x^2)
    ex2 = @NLexpression(m, ex + 1)
    @test sprint(show, ex2) == "subexpression[2]: subexpression[1] + 1.0"
    return
end

function test_parse_parameters()
    m = Model()
    @NLparameter(m, param == 10)
    ex = @NLexpression(m, param + 1)
    @test sprint(show, ex) == "subexpression[1]: param + 1.0"
    return
end

function test_parse_user_defined_function_univariate()
    model = Model()
    @variable(model, x)
    user_function = x -> x
    JuMP.register(model, :f, 1, user_function, autodiff = true)
    ex = @NLexpression(model, f(x))
    @test sprint(show, ex) == "subexpression[1]: f(x)"
    return
end

function test_parse_user_defined_function_multivariate()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    user_function = (x, y) -> x
    JuMP.register(model, :f, 2, user_function, autodiff = true)
    ex = @NLexpression(model, f(x, y))
    @test sprint(show, ex) == "subexpression[1]: f(x, y)"
    return
end

function test_parse_splatting()
    model = Model()
    @variable(model, x[1:2])
    user_function = (x, y) -> x
    JuMP.register(model, :f, 2, user_function, autodiff = true)
    ex = @NLexpression(model, f(x...))
    @test sprint(show, ex) == "subexpression[1]: f(x[1], x[2])"
    return
end

function test_parse_mixed_splatting()
    model = Model()
    @variable(model, x[1:2])
    @variable(model, y)
    @variable(model, z[1:1])
    ex = @NLexpression(model, *(x..., y, z...))
    @test sprint(show, ex) == "subexpression[1]: x[1] * x[2] * y * z[1]"
    return
end

function test_error_splatting_non_symbols()
    model = Model()
    @variable(model, x[1:2])
    @test_throws ErrorException @NLexpression(model, (*)((x / 2)...))
    return
end

function test_error_on_begin_end()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "`begin...end` blocks are not supported in nonlinear macros. The " *
        "nonlinear expression must be a single statement.",
    )
    @test_macro_throws(err, @NLobjective(model, Max, begin
        sin(x) + 1
    end))
    return
end

function test_error_on_unexpected_splatting()
    model = Model()
    @variable(model, x[1:2])
    @test_throws(ErrorException, @NLexpression(model, x...))
    return
end

function test_error_on_getfield_in_expression()
    model = Model()
    @variable(model, x[1:2])
    @test_macro_throws(
        ErrorException,
        @NLexpression(model, sum(foo.bar(i) * x[i] for i in 1:2))
    )
    return
end

function test_error_on_sum()
    m = Model()
    x = [1, 2, 3]
    @test_throws ErrorException @NLexpression(m, sum(x))
    return
end

function test_error_on_non_scalar_expression()
    m = Model()
    x = [1, 2, 3]
    @test_throws ErrorException @NLexpression(m, x + 1)
    return
end

# Converts the lower-triangular sparse Hessian in MOI format into a dense
# matrix.
function _dense_hessian(hessian_sparsity, V, n)
    I = [i for (i, j) in hessian_sparsity]
    J = [j for (i, j) in hessian_sparsity]
    raw = SparseArrays.sparse(I, J, V, n, n)
    return Matrix(
        raw + raw' -
        SparseArrays.sparse(LinearAlgebra.diagm(0 => LinearAlgebra.diag(raw))),
    )
end

# Converts the sparse Jacobian in MOI format into a dense matrix.
function _dense_jacobian(jacobian_sparsity, V, m, n)
    I = [i for (i, j) in jacobian_sparsity]
    J = [j for (i, j) in jacobian_sparsity]
    raw = SparseArrays.sparse(I, J, V, m, n)
    return Matrix(raw)
end

function test_hessian_evaluation_issue_435()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)
    @NLexpression(m, foo, a * b + c^2)
    @NLobjective(m, Min, foo)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    values = [1.0, 2.0, 3.0] # Values for a, b, and c, respectively.
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 3) ≈
          [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0]
    # make sure we don't get NaNs in this case
    @NLobjective(m, Min, a * b + 3 * c^2)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess])
    values = [1.0, 2.0, -1.0]
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 3) ≈
          [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0]
    # Initialize again
    MOI.initialize(d, [:Hess])
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 3) ≈
          [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0]
    return
end

function test_NaN_corner_case_issue_695()
    m = Model()
    x0 = 0.0
    y0 = 0.0
    @variable(m, x)
    @variable(m, y)
    @NLobjective(m, Min, (x - x0) / (sqrt(y0) + sqrt(y)))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:HessVec])
    h = ones(2)
    v = [2.4, 3.5]
    values = [1.0, 2.0] # For x and y.
    MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, Float64[])
    correct = [0.0 -1/(2*2^(3/2)); -1/(2*2^(3/2)) 3/(4*2^(5/2))] * v
    @test h ≈ correct
    return
end

function test_NaN_corner_case_issue_1205()
    m = Model()
    @variable(m, x)
    @NLobjective(m, Min, x^1.0)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    values = zeros(1)
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 1) ≈ [0.0]
    return
end

function test_NaN_corner_case_ifelse_issue_1205()
    m = Model()
    @variable(m, x)
    @NLobjective(m, Min, ifelse(true, x, x^1.0))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    values = zeros(1)
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 1) ≈ [0.0]
    return
end

function test_prod_corner_case_issue_1181()
    model = Model()
    @variable(model, x[1:2])
    @NLobjective(model, Min, x[1] * x[2])
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, V, [0.659, 0.702], 1.0, Float64[])
    @test V == [0.0, 0.0, 1.0]
    return
end

function test_constant_ifelse_issue_2115()
    model = Model()
    @variable(model, x)
    @NLobjective(model, Min, ifelse(x >= 1, 1, 0))
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    @test length(hessian_sparsity) == 0
    return
end

function test_hessians_and_hessvec()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)
    @NLobjective(m, Min, a * b + c^2)
    @NLconstraint(m, c * b <= 1)
    @NLconstraint(m, a^2 / 2 <= 1)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:HessVec, :Hess])
    values = [1.0, 2.0, 3.0] # For a, b, c.
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, [2.0, 3.0])
    correct_hessian = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]
    @test _dense_hessian(hessian_sparsity, V, 3) ≈ correct_hessian
    h = ones(3) # The input values should be overwritten.
    v = [2.4, 3.5, 1.2]
    MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0, 3.0])
    @test h ≈ correct_hessian * v
    return
end

function test_hessians_and_hessvec_with_subexpressions()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)
    @NLexpression(m, ab, a * b)
    @NLobjective(m, Min, ab + c^2)
    @NLconstraint(m, c * b <= 1)
    @NLconstraint(m, a^2 / 2 <= 1)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:HessVec, :Hess])
    values = [1.0, 2.0, 3.0] # For a, b, c.
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, [2.0, 3.0])
    correct_hessian = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]
    @test _dense_hessian(hessian_sparsity, V, 3) ≈ correct_hessian
    h = ones(3) # The input values should be overwritten.
    v = [2.4, 3.5, 1.2]
    MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0, 3.0])
    @test h ≈ correct_hessian * v
    return
end

function test_jacobians_and_jacvec()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)
    @NLobjective(m, Min, a * b + c^2)
    @NLconstraint(m, c * b <= 1)
    @NLconstraint(m, a^2 / 2 <= 1)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:JacVec, :Jac])
    values = [1.0, 2.0, 3.0] # For a, b, c.
    jacobian_sparsity = MOI.jacobian_structure(d)
    V = zeros(length(jacobian_sparsity))
    MOI.eval_constraint_jacobian(d, V, values)
    correct_jacobian = [0.0 3.0 2.0; 1.0 0.0 0.0]
    @test _dense_jacobian(jacobian_sparsity, V, 2, 3) ≈ correct_jacobian
    v = [2.4, 3.5, 1.2]
    product_storage = zeros(2)
    MOI.eval_constraint_jacobian_product(d, product_storage, values, v)
    @test product_storage ≈ correct_jacobian * v
    w = [0.6, 4.3]
    product_storage = zeros(3)
    MOI.eval_constraint_jacobian_transpose_product(
        d,
        product_storage,
        values,
        w,
    )
    @test product_storage ≈ correct_jacobian' * w
    return
end

function test_jacobians_and_jacvec_with_subexpressions()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)
    @NLexpression(m, bc, b * c)
    @NLobjective(m, Min, a * b + c^2)
    @NLconstraint(m, bc <= 1)
    @NLconstraint(m, a^2 / 2 <= 1)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:JacVec, :Jac])
    values = [1.0, 2.0, 3.0] # For a, b, c.
    jacobian_sparsity = MOI.jacobian_structure(d)
    V = zeros(length(jacobian_sparsity))
    MOI.eval_constraint_jacobian(d, V, values)
    correct_jacobian = [0.0 3.0 2.0; 1.0 0.0 0.0]
    @test _dense_jacobian(jacobian_sparsity, V, 2, 3) ≈ correct_jacobian
    v = [2.4, 3.5, 1.2]
    product_storage = zeros(2)
    MOI.eval_constraint_jacobian_product(d, product_storage, values, v)
    @test product_storage ≈ correct_jacobian * v
    w = [0.6, 4.3]
    product_storage = zeros(3)
    MOI.eval_constraint_jacobian_transpose_product(
        d,
        product_storage,
        values,
        w,
    )
    @test product_storage ≈ correct_jacobian' * w
    return
end

function test_expression_graphs()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    @NLobjective(m, Min, x^2 + y^2)
    @NLexpression(m, ex, exp(x))
    @NLconstraint(m, ex - y == 0)
    @NLconstraint(m, ex + 1 == 0)

    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:ExprGraph])
    xidx = x.index
    yidx = y.index
    @test MOI.objective_expr(d) == :(x[$xidx]^2.0 + x[$yidx]^2.0)
    @test MOI.constraint_expr(d, 1) ==
          :((exp(x[$xidx]) - x[$yidx]) - 0.0 == 0.0)
    @test MOI.constraint_expr(d, 2) == :((exp(x[$xidx]) + 1) - 0.0 == 0.0)
    return
end

function test_more_expression_graphs()
    m = Model()
    @variable(m, x)
    @variable(m, y)
    ψ(x) = 1
    t(x, y) = 2
    JuMP.register(m, :ψ, 1, ψ, autodiff = true)
    JuMP.register(m, :t, 2, t, autodiff = true)
    @NLobjective(m, Min, x^y)
    @NLconstraint(m, sin(x) * cos(y) == 5)
    @NLconstraint(m, nlconstr[i = 1:2], i * x^2 == i)
    @NLconstraint(m, -0.5 <= sin(x) <= 0.5)
    @NLconstraint(m, ψ(x) + t(x, y) <= 3)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:ExprGraph])
    xidx = x.index
    yidx = y.index
    @test MOI.objective_expr(d) == :(x[$xidx]^x[$yidx])
    @test MOI.constraint_expr(d, 1) ==
          :(sin(x[$xidx]) * cos(x[$yidx]) - 5 == 0.0)
    @test MOI.constraint_expr(d, 2) == :(1.0 * x[$xidx]^2 - 1.0 == 0.0)
    @test MOI.constraint_expr(d, 3) == :(2.0 * x[$xidx]^2 - 2.0 == 0.0)
    @test MOI.constraint_expr(d, 4) == :(-0.5 <= sin(x[$xidx]) <= 0.5)
    @test MOI.constraint_expr(d, 5) ==
          :(ψ(x[$xidx]) + t(x[$xidx], x[$yidx]) - 3.0 <= 0.0)
    return
end

function test_expression_graph_for_ifelse()
    m = Model()
    @variable(m, x)
    @NLobjective(m, Min, ifelse(x <= 1, x^2, x))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:ExprGraph])
    xidx = x.index
    @test MOI.objective_expr(d) ==
          :(ifelse(x[$xidx] <= 1, x[$xidx]^2, x[$xidx]))
    return
end

function test_expression_graph_for_empty_sum_and_prod()
    m = Model()
    @variable(m, x)
    @NLconstraint(m, x <= sum(0 for i in []) + prod(1 for i in []))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:ExprGraph])
    xidx = x.index
    @test MOI.constraint_expr(d, 1) == :((x[$xidx] - (0.0 + 1.0)) - 0.0 <= 0.0)
    return
end

function test_NLparameter()
    model = Model()
    @NLparameter(model, p == 1.0)
    @test JuMP.value(p) == 1.0
    return
end

function test_NLparameter_set_value()
    model = Model()
    @NLparameter(model, p == 1.0)
    JuMP.set_value(p, 10.0)
    @test JuMP.value(p) == 10.0
    return
end

function test_NLconstraints()
    model = Model()
    @variable(model, 0 <= x <= 1)
    @variable(model, y[1:3])
    @objective(model, Max, x)
    @NLconstraints(model, begin
        ref[i = 1:3], y[i] == 0
        x + y[1] * y[2] * y[3] <= 0.5
    end)
    @test JuMP.num_nonlinear_constraints(model) == 4
    evaluator = JuMP.NLPEvaluator(model)
    MOI.initialize(evaluator, [:ExprGraph])
    for i in 1:3
        @test MOI.constraint_expr(evaluator, i) ==
              :(x[$(y[i].index)] - 0.0 == 0.0)
    end
    @test MOI.constraint_expr(evaluator, 4) == :(
        (
            x[$(x.index)] +
            x[$(y[1].index)] * x[$(y[2].index)] * x[$(y[3].index)]
        ) - 0.5 <= 0.0
    )
    return
end

# This covers the code that computes Hessians in odd chunks of Hess-vec
# products.
function test_dense_Hessian()
    m = Model()
    @variable(m, x[1:18])
    @NLobjective(m, Min, prod(x[i] for i in 1:18))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)
    V = zeros(length(hessian_sparsity))
    values = ones(18)
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 18) ≈
          ones(18, 18) - LinearAlgebra.diagm(0 => ones(18))
    values[1] = 0.5
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    @test _dense_hessian(hessian_sparsity, V, 18) ≈ [
        0 ones(17)'
        ones(17) (ones(17, 17)-LinearAlgebra.diagm(0 => ones(17)))/2
    ]
    return
end

function test_eval_objective_and_eval_objective_gradient()
    m = Model()
    @variable(m, x[1:4])
    @NLparameter(m, p == 2)
    @NLexpression(m, ex, p * x[1])

    ψ(x) = sin(x)
    t(x, y) = x + 3y
    JuMP.register(m, :ψ, 1, ψ, autodiff = true)
    JuMP.register(m, :t, 2, t, autodiff = true)

    @NLobjective(m, Min, ex / 2 + sin(x[2]) / ψ(x[2]) + t(x[3], x[4]))
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Grad])
    variable_values = fill(2.0, (4,))
    @test MOI.eval_objective(d, variable_values) ≈
          variable_values[1] + 1 + variable_values[3] + 3variable_values[4]
    grad = zeros(4)
    MOI.eval_objective_gradient(d, grad, variable_values)
    @test grad ≈ [1.0, 0.0, 1.0, 3.0]
    return
end

function test_eval_constraint_and_jacobians()
    m = Model()
    @variable(m, x[1:4])
    @NLparameter(m, p == 2)
    @NLexpression(m, ex, p * x[1])

    ψ(x) = sin(x)
    t(x, y) = x + 3y
    JuMP.register(m, :ψ, 1, ψ, autodiff = true)
    JuMP.register(m, :t, 2, t, autodiff = true)

    @NLconstraint(m, Min, ex / 2 + sin(x[2]) / ψ(x[2]) + t(x[3], x[4]) <= 0.0)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Jac])
    variable_values = fill(2.0, (4,))
    constraint_value = zeros(1)
    MOI.eval_constraint(d, constraint_value, variable_values)
    @test constraint_value[1] ≈
          variable_values[1] + 1 + variable_values[3] + 3variable_values[4]
    jacobian_sparsity = MOI.jacobian_structure(d)
    I = [i for (i, j) in jacobian_sparsity]
    J = [j for (i, j) in jacobian_sparsity]
    @test all(I .== 1)
    jac_nonzeros = zeros(length(J))
    MOI.eval_constraint_jacobian(d, jac_nonzeros, variable_values)
    jac_values = zeros(4)
    jac_values[J] = jac_nonzeros
    @test jac_values ≈ [1.0, 0.0, 1.0, 3.0]
    return
end

function test_non_macro_nonlinear_functions()
    model = Model()
    @variable(model, x)
    @variable(model, y)
    @expression(model, aff, x + 2y - 3)
    @expression(model, quad, x^2 + 2y^2 - x)
    nlexpr = JuMP.add_nonlinear_expression(model, :($x^2 + $y^2))
    JuMP.set_nonlinear_objective(model, MIN_SENSE, :(2 * $nlexpr))
    JuMP.add_nonlinear_constraint(model, :($x + $y <= 1))
    JuMP.add_nonlinear_constraint(model, :($x + $y >= 1))
    JuMP.add_nonlinear_constraint(model, :($x + $y == 1))
    JuMP.add_nonlinear_constraint(model, :(0 <= $x + $y <= 1))
    JuMP.add_nonlinear_constraint(model, :($aff == 1))
    JuMP.add_nonlinear_constraint(model, :($quad == 1))

    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, [:ExprGraph])
    xidx = x.index
    yidx = y.index
    @test MOI.objective_expr(d) == :(2.0 * (x[$xidx]^2.0 + x[$yidx]^2.0))
    @test MOI.constraint_expr(d, 1) == :((x[$xidx] + x[$yidx]) - 1.0 <= 0.0)
    @test MOI.constraint_expr(d, 2) == :((x[$xidx] + x[$yidx]) - 1.0 >= 0.0)
    @test MOI.constraint_expr(d, 3) == :((x[$xidx] + x[$yidx]) - 1.0 == 0.0)
    @test MOI.constraint_expr(d, 4) == :(0.0 <= x[$xidx] + x[$yidx] <= 1.0)
    @test MOI.constraint_expr(d, 5) ==
          :((-3.0 + x[$xidx] + 2.0 * x[$yidx]) - 1.0 == 0.0)
    @test MOI.constraint_expr(d, 6) == :(
        (+(-1.0 * x[$xidx]) + x[$xidx] * x[$xidx] + x[$yidx] * x[$yidx] * 2.0) -
        1.0 == 0.0
    )
    return
end

function test_views_on_Hessian_functions()
    m = Model()
    @variable(m, a)
    @variable(m, b)
    @variable(m, c)

    @NLexpression(m, foo, a * b + c^2)

    @NLobjective(m, Min, foo)
    d = JuMP.NLPEvaluator(m)
    MOI.initialize(d, [:Hess, :HessVec])
    hessian_sparsity = MOI.hessian_lagrangian_structure(d)

    V = zeros(4)
    values = [1.0, 2.0, 3.0] # Values for a, b, and c, respectively.
    MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
    h = ones(3)
    v = [2.4; 3.5; 4.6]
    MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, Float64[])

    values2 = zeros(10)
    values2[5:7] = values
    values_view = @view values2[5:7]
    V2 = zeros(10)
    V_view = @view V2[4:7]
    MOI.eval_hessian_lagrangian(d, V_view, values_view, 1.0, Float64[])
    @test V_view == V

    h2 = zeros(10)
    h_view = @view h2[3:5]
    v2 = zeros(10)
    v2[4:6] = v
    v_view = @view v2[4:6]
    MOI.eval_hessian_lagrangian_product(
        d,
        h_view,
        values_view,
        v_view,
        1.0,
        Float64[],
    )
    @test h_view == h
    return
end

function test_constant_expressions()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, 10)
    @NLobjective(model, Min, expr + x)
    d = JuMP.NLPEvaluator(model)
    MOI.initialize(d, [:Grad])
    grad = zeros(1)
    MOI.eval_objective_gradient(d, grad, [2.0])
    @test grad == [1.0]
    return
end

function test_user_defined_function_with_variable_closure()
    model = Model()
    @variable(model, x[1:2])
    f(x1) = x1 + x[2]
    @test_throws(
        ErrorException(
            "Expected return type of `Float64` from the user-defined " *
            "function :f, but got `$(AffExpr)`.",
        ),
        register(model, :f, 1, f; autodiff = true),
    )
    return
end

function test_user_defined_function_with_variable_closure_after_register()
    model = Model()
    @variable(model, x)
    f(y) = y < 1 ? y : y + x
    register(model, :f, 1, f; autodiff = true)
    @NLobjective(model, Min, f(x))
    d = NLPEvaluator(model)
    MOI.initialize(d, [:Grad])
    err = ErrorException(
        "Expected return type of Float64 from a user-defined function, " *
        "but got $(typeof(1.0 + x)). Make sure your user-defined " *
        "function only depends on variables passed as arguments.",
    )
    @test_throws(err, MOI.eval_objective(d, [2.0]))
    return
end

function test_user_defined_function_returning_bad_type_after_register()
    model = Model()
    @variable(model, x)
    f(x) = x < 1 ? x : string(x)
    register(model, :f, 1, f; autodiff = true)
    @NLobjective(model, Min, f(x))
    d = NLPEvaluator(model)
    MOI.initialize(d, [:Grad])
    err = ErrorException(
        "Expected return type of Float64 from a user-defined function, " *
        "but got String.",
    )
    @test_throws(err, MOI.eval_objective(d, [2.0]))
    return
end

function test_user_defined_function_returning_bad_type()
    model = Model()
    @variable(model, x)
    f(x) = string(x)
    @test_throws(
        ErrorException(
            "Expected return type of `Float64` from the user-defined " *
            "function :f, but got `String`.",
        ),
        register(model, :f, 1, f; autodiff = true),
    )
    return
end

function test_value_on_NonlinearExpressions()
    model = Model()
    @variable(model, x)
    @NLexpression(model, ex1, sin(x))
    @NLexpression(model, ex2, ex1 + x^2)
    @NLexpression(model, ex3, 2ex1 + ex2 / 2)
    JuMP.set_start_value(x, 2.0)
    @test JuMP.value(JuMP.start_value, ex1) ≈ sin(2.0)
    @test JuMP.value(JuMP.start_value, ex2) ≈ sin(2.0) + 4.0
    @test JuMP.value(JuMP.start_value, ex3) ≈ 2.5 * sin(2.0) + 2.0
    return
end

function test_hessians_disabled_with_user_defined_multivariate_functions()
    model = Model()
    my_f(x, y) = (x - 1)^2 + (y - 2)^2
    JuMP.register(model, :my_f, 2, my_f, autodiff = true)
    @variable(model, x[1:2])
    @NLobjective(model, Min, my_f(x[1], x[2]))
    evaluator = JuMP.NLPEvaluator(model)
    @test !(:Hess in MOI.features_available(evaluator))
    return
end

function test_AffExpr_in_nonlinear()
    model = Model()
    @variable(model, x, start = 1.1)
    @variable(model, y, start = 1.2)
    @expression(model, ex, 2 * x + y + 1)
    nl_ex = @NLexpression(model, ex^2)
    @test isapprox(
        value(start_value, nl_ex),
        (2 * 1.1 + 1.2 + 1)^2,
        atol = 1e-4,
    )
    return
end

function test_QuadExpr_in_nonlinear()
    model = Model()
    @variable(model, x, start = 1.1)
    @variable(model, y, start = 1.2)
    @expression(model, ex, 0.5 * x^2 + y^2 + 2 * x + 1)
    nl_ex = @NLexpression(model, sqrt(ex))
    @test isapprox(
        value(start_value, nl_ex),
        sqrt(0.5 * 1.1^2 + 1.2^2 + 2 * 1.1 + 1),
        atol = 1e-4,
    )
    return
end

function test_error_on_complex_values()
    model = Model()
    @variable(model, x)
    c = sqrt(Complex(-1))
    expected_exception = ErrorException(
        "Unexpected object $c of type $(typeof(c)) in nonlinear expression.",
    )
    @test_throws expected_exception @NLobjective(model, Min, c * x)
    return
end

function test_SpecialFunctions()
    model = Model()
    @variable(model, x)
    @NLconstraint(model, c1, erf(x) <= 0.0)
    d = NLPEvaluator(model)
    MOI.initialize(d, Symbol[:Grad])
    out = zeros(1)
    MOI.eval_constraint(d, out, [2.0])
    @test out[1] ≈ 0.9953222 atol = 1e-6
    return
end

function test_JuMP_extensions()
    model = JuMPExtension.MyModel()
    @variable(model, x)
    err = ErrorException(
        "Encountered an error parsing nonlinear expression: we don't support " *
        "models of type $(typeof(model)). In general, JuMP's nonlinear features " *
        "don't work with JuMP-extensions.",
    )
    @test_throws(err, @NLexpression(model, sqrt(x)))
    return
end

function test_rad2deg_and_deg2rad()
    data = JuMP.Nonlinear.NonlinearData()
    x = 1.0
    operators = data.operators
    @test JuMP.Nonlinear.eval_univariate_hessian(operators, :rad2deg, x) == 0.0
    @test JuMP.Nonlinear.eval_univariate_hessian(operators, :deg2rad, x) == 0.0
    return
end

function test_reregister_univariate()
    model = Model()
    @variable(model, x >= 0)
    f(x) = x^2
    register(model, :f, 1, f; autodiff = true)
    @test_throws ErrorException register(model, :f, 1, f; autodiff = true)
    return
end

function test_reregister_multivariate()
    model = Model()
    @variable(model, x >= 0)
    f(x, y) = x + y
    register(model, :f, 2, f; autodiff = true)
    @test_throws ErrorException register(model, :f, 2, f; autodiff = true)
    return
end

function test_univariate_NLconstraint_is_valid()
    model = Model()
    model2 = Model()
    @variable(model, x >= 0)
    c = @NLconstraint(model, exp(x) <= 1)
    @test is_valid(model, c)
    @test !is_valid(model2, c)
    return
end

function test_multivariate_NLconstraint_is_valid()
    model = Model()
    model2 = Model()
    @variable(model, x >= 0)
    @variable(model, y >= 0)
    c = @NLconstraint(model, exp(x) + log(y) <= 1)
    @test is_valid(model, c)
    @test !is_valid(model2, c)
    return
end

function test_broadcast_error_NLexpression()
    model = Model()
    @variable(model, x)
    @NLexpression(model, expr, sin(x))
    @test_throws(
        ErrorException(
            "`JuMP.value` is not defined for collections of JuMP types. " *
            "Use Julia's broadcast syntax instead: `JuMP.value.(x)`.",
        ),
        value([x, x]),
    )
    return
end

function test_interval_errors()
    model = Model()
    @variable(model, x)
    err = ErrorException(
        "Interval constraint contains non-constant left- or right-hand " *
        "sides. Reformulate as two separate constraints, or move all " *
        "variables into the central term.",
    )
    @test_throws err add_nonlinear_constraint(model, :($x <= $x <= 2 * $x))
    return
end

function test_dual_start_value()
    model = Model()
    @variable(model, x)
    @variable(model, t >= 0)
    @NLconstraint(model, 2x >= 1)
    @NLconstraint(model, t <= sqrt(x))
    @test nonlinear_dual_start_value(model) === nothing
    set_nonlinear_dual_start_value(model, [1.0, -1.0])
    @test nonlinear_dual_start_value(model) == [1.0, -1.0]
    set_nonlinear_dual_start_value(model, nothing)
    @test nonlinear_dual_start_value(model) === nothing
    @test_throws(ArgumentError, set_nonlinear_dual_start_value(model, [1.0]),)
    return
end

function test_user_defined_function_checked_error_univariate()
    function f(x)
        if x >= 1
            y = zeros(1)
            y[1] = (x - 1)^2
            return y[1]
        else
            return (x - 1)^2
        end
    end
    model = Model()
    @variable(model, x >= 0)
    register(model, :f, 1, f; autodiff = true)
    @NLobjective(model, Min, f(x))
    nlp = NLPEvaluator(model)
    MOI.initialize(nlp, Symbol[:Grad])
    g = [NaN]
    MOI.eval_objective_gradient(nlp, g, [0.0])
    @test g == [-2.0]
    err = ErrorException(
        "JuMP's autodiff of the user-defined function f failed with a " *
        "MethodError.\n\n$(JuMP.Nonlinear._FORWARD_DIFF_METHOD_ERROR_HELPER)",
    )
    @test_throws(err, MOI.eval_objective_gradient(nlp, g, [2.0]))
    return
end

function test_user_defined_function_checked_error_univariate_other_error()
    function f(x)
        if typeof(x) != Float64
            error("not differentiable")
        end
        return (x - 1)^2
    end
    model = Model()
    @variable(model, x)
    register(model, :f, 1, f; autodiff = true)
    @NLobjective(model, Min, f(x))
    nlp = NLPEvaluator(model)
    MOI.initialize(nlp, Symbol[:Grad])
    err = ErrorException("not differentiable")
    @test_throws(err, MOI.eval_objective_gradient(nlp, [NaN], [2.0]))
    return
end

function test_user_defined_function_checked_error_univariate()
    f(x) = (x - 1)^2
    function ∇f(x)
        if x >= 1
            y = zeros(1)
            y[1] = 2 * (x - 1)
            return y[1]
        else
            return 2 * (x - 1)
        end
    end
    model = Model()
    @variable(model, x >= 0)
    register(model, :f, 1, f, ∇f; autodiff = true)
    @NLobjective(model, Min, f(x))
    nlp = NLPEvaluator(model)
    MOI.initialize(nlp, Symbol[:Hess])
    H = [NaN]
    MOI.eval_hessian_lagrangian(nlp, H, [0.0], 1.0, Float64[])
    @test H == [2.0]
    err = ErrorException(
        "JuMP's autodiff of the user-defined function f failed with a " *
        "MethodError.\n\n$(JuMP.Nonlinear._FORWARD_DIFF_METHOD_ERROR_HELPER)",
    )
    @test_throws(
        err,
        MOI.eval_hessian_lagrangian(nlp, H, [2.0], 1.0, Float64[]),
    )
    return
end

function test_user_defined_function_checked_error_multivariate()
    function f(x...)
        if x[1] >= 1
            y = zeros(1)
            y[1] = (x[1] - 1)^2
            return y[1]
        else
            return (x[1] - 1)^2 + x[2]
        end
    end
    model = Model()
    @variable(model, x[1:2] >= 0)
    register(model, :f, 2, f; autodiff = true)
    @NLobjective(model, Min, f(x...))
    nlp = NLPEvaluator(model)
    MOI.initialize(nlp, Symbol[])
    g = [NaN, NaN]
    MOI.eval_objective_gradient(nlp, g, [0.0, 1.0])
    @test g == [-2.0, 1.0]
    err = ErrorException(
        "JuMP's autodiff of the user-defined function f failed with a " *
        "MethodError.\n\n$(JuMP.Nonlinear._FORWARD_DIFF_METHOD_ERROR_HELPER)",
    )
    @test_throws(err, MOI.eval_objective_gradient(nlp, g, [2.0, 1.0]))
    return
end

end

TestNLP.runtests()
