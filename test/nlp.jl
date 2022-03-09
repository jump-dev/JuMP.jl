using JuMP
using LinearAlgebra
using SparseArrays
using Test

include(joinpath(@__DIR__, "utilities.jl"))

@static if !(:JuMPExtension in names(Main))
    include(joinpath(@__DIR__, "JuMPExtension.jl"))
end

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

@testset "register-ForwardDiff-incompatibility" begin
    test_register_check_forwarddiff_univariate_f()
    test_register_check_forwarddiff_univariate_gradf()
    test_register_check_forwarddiff_multivariate()
    test_register_check_forwarddiff_multivariate_gradf()
end

@testset "Auto-register-univariate" begin
    test_univariate_error()
    test_univariate_error_existing()
    test_univariate()
    test_univariate_existing_nlpdata()
    test_univariate_redefine()
    test_univariate_register_twice()
    test_univariate_register_twice_error()
end

@testset "Auto-register-multivariate" begin
    test_multivariate_error()
    test_multivariate_error_existing()
    test_multivariate()
    test_multivariate_existing_nlpdata()
    test_multivariate_redefine()
    test_multivariate_register_warn()
    test_multivariate_register_splat()
    test_multivariate_register_splat_existing()
    test_multivariate_max()
    test_multivariate_min()
end

@testset "Nonlinear" begin
    import JuMP: _NonlinearExprData

    function expressions_equal(ex1::_NonlinearExprData, ex2::_NonlinearExprData)
        return ex1.nd == ex2.nd && ex1.const_values == ex2.const_values
    end

    # TODO: These are poorly designed tests because they would not catch errors
    # that affect both _process_NL_expr and _NonlinearExprData's parsing (although
    # they use different pathways). It would be better to check the NodeData
    # representation directly.

    @testset "all_nonlinear_constraints" begin
        model = Model()
        @variable(model, x)
        @NLconstraint(model, c1, x^2 <= 1)
        c2 = @NLconstraint(model, [i = 1:2], x^i >= -1)
        @test all_nonlinear_constraints(model) == [c1; c2]
    end

    @testset "Parse + (binary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, x + y),
            _NonlinearExprData(m, :($x + $y)),
        )
    end

    @testset "Parse + (ternary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, x + y + z),
            _NonlinearExprData(m, :($x + $y + $z)),
        )
    end

    @testset "Parse * (binary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, x * y),
            _NonlinearExprData(m, :($x * $y)),
        )
    end

    @testset "Parse * (ternary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, x * y * z),
            _NonlinearExprData(m, :($x * $y * $z)),
        )
    end

    @testset "Parse ^ (binary)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, x^3),
            _NonlinearExprData(m, :($x^3)),
        )
    end

    @testset "Parse sin (univariate)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, sin(x)),
            _NonlinearExprData(m, :(sin($x))),
        )
    end

    @testset "Parse ifelse" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(
            JuMP.@_process_NL_expr(
                m,
                ifelse(1 == 2 || 3 == 4 && 5 == 6, x, 0.0)
            ),
            _NonlinearExprData(
                m,
                :(ifelse(1 == 2 || 3 == 4 && 5 == 6, $x, 0.0)),
            ),
        )
    end

    @testset "Parse ifelse (3-way comparison)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, ifelse(1 <= 2 <= 3, x, 0.0)),
            _NonlinearExprData(m, :(ifelse(1 <= 2 <= 3, $x, 0.0))),
        )
    end

    @testset "Parse sum" begin
        m = Model()
        @variable(m, x[1:2])
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, sum(x[i] for i in 1:2)),
            _NonlinearExprData(m, :($(x[1]) + $(x[2]))),
        )
    end

    @testset "Parse prod" begin
        m = Model()
        @variable(m, x[1:2])
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, prod(x[i] for i in 1:2)),
            _NonlinearExprData(m, :($(x[1]) * $(x[2]))),
        )
    end

    @testset "Parse subexpressions" begin
        m = Model()
        @variable(m, x)
        @NLexpression(m, ex, x^2)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, ex + 1),
            _NonlinearExprData(m, :($ex + 1)),
        )
    end

    @testset "Parse parameters" begin
        m = Model()
        @NLparameter(m, param == 10)
        @test expressions_equal(
            JuMP.@_process_NL_expr(m, param + 1),
            _NonlinearExprData(m, :($param + 1)),
        )
    end

    @testset "Parse user-defined function (univariate)" begin
        model = Model()
        @variable(model, x)
        user_function = x -> x
        JuMP.register(model, :f, 1, user_function, autodiff = true)
        @test expressions_equal(
            JuMP.@_process_NL_expr(model, f(x)),
            _NonlinearExprData(model, :(f($x))),
        )
    end

    @testset "Parse user-defined function (multivariate)" begin
        model = Model()
        @variable(model, x)
        @variable(model, y)
        user_function = (x, y) -> x
        JuMP.register(model, :f, 2, user_function, autodiff = true)
        @test expressions_equal(
            JuMP.@_process_NL_expr(model, f(x, y)),
            _NonlinearExprData(model, :(f($x, $y))),
        )
    end

    @testset "Parse splatting" begin
        model = Model()
        @variable(model, x[1:2])
        user_function = (x, y) -> x
        JuMP.register(model, :f, 2, user_function, autodiff = true)
        @test expressions_equal(
            JuMP.@_process_NL_expr(model, f(x...)),
            _NonlinearExprData(model, :(f($(x[1]), $(x[2])))),
        )
    end

    @testset "Parse mixed splatting" begin
        model = Model()
        @variable(model, x[1:2])
        @variable(model, y)
        @variable(model, z[1:1])
        @test expressions_equal(
            JuMP.@_process_NL_expr(model, (*)(x..., y, z...)),
            _NonlinearExprData(model, :((*)($(x[1]), $(x[2]), $y, $(z[1])))),
        )
    end

    @testset "Error on splatting non-symbols" begin
        model = Model()
        @variable(model, x[1:2])
        @test_macro_throws ErrorException @NLexpression(model, (*)((x / 2)...))
    end

    @testset "Error on begin...end" begin
        model = Model()
        @variable(model, x)
        err = ErrorException(
            "`begin...end` blocks are not supported in nonlinear macros. The " *
            "nonlinear expression must be a single statement.",
        )
        @test_macro_throws(err, @NLobjective(model, Max, begin
            sin(x) + 1
        end))
    end

    @testset "Error on unexpected splatting" begin
        model = Model()
        @variable(model, x[1:2])
        @test_macro_throws ErrorException @NLexpression(model, x...)
    end

    @testset "Error on x.f(y) in NL expression" begin
        model = Model()
        @variable(model, x[1:2])
        @test_throws(
            ErrorException,
            @NLexpression(model, sum(foo.bar(i) * x[i] for i in 1:2))
        )
    end

    @testset "Error on sum(x)" begin
        m = Model()
        x = [1, 2, 3]
        @test_throws ErrorException @NLexpression(m, sum(x))
    end

    @testset "Error on non-scalar expressions" begin
        m = Model()
        x = [1, 2, 3]
        @test_throws ErrorException @NLexpression(m, x + 1)
    end

    # Converts the lower-triangular sparse Hessian in MOI format into a dense
    # matrix.
    function dense_hessian(hessian_sparsity, V, n)
        I = [i for (i, j) in hessian_sparsity]
        J = [j for (i, j) in hessian_sparsity]
        raw = sparse(I, J, V, n, n)
        return Matrix(raw + raw' - sparse(diagm(0 => diag(raw))))
    end

    # Converts the sparse Jacobian in MOI format into a dense matrix.
    function dense_jacobian(jacobian_sparsity, V, m, n)
        I = [i for (i, j) in jacobian_sparsity]
        J = [j for (i, j) in jacobian_sparsity]
        raw = sparse(I, J, V, m, n)
        return Matrix(raw)
    end

    @testset "Hessian evaluation (Issue #435)" begin
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
        @test dense_hessian(hessian_sparsity, V, 3) ≈
              [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0]

        # make sure we don't get NaNs in this case
        @NLobjective(m, Min, a * b + 3 * c^2)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:Hess])
        values = [1.0, 2.0, -1.0]
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 3) ≈
              [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0]

        # Initialize again
        MOI.initialize(d, [:Hess])
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 3) ≈
              [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0]
    end

    @testset "NaN corner case (Issue #695)" begin
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
    end

    @testset "NaN corner case (Issue #1205)" begin
        m = Model()
        @variable(m, x)

        @NLobjective(m, Min, x^1.0)

        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = zeros(1)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 1) ≈ [0.0]
    end

    @testset "NaN corner case - ifelse (Issue #1205)" begin
        m = Model()
        @variable(m, x)

        @NLobjective(m, Min, ifelse(true, x, x^1.0))

        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = zeros(1)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 1) ≈ [0.0]
    end

    @testset "Product corner case (issue #1181)" begin
        model = Model()
        @variable(model, x[1:2])

        @NLobjective(model, Min, x[1] * x[2])
        d = JuMP.NLPEvaluator(model)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, [0.659, 0.702], 1.0, Float64[])
        @test V == [0.0, 0.0, 1.0]
    end

    @testset "Constant ifelse (Issue #2115)" begin
        model = Model()
        @variable(model, x)
        @NLobjective(model, Min, ifelse(x >= 1, 1, 0))

        d = JuMP.NLPEvaluator(model)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        @test length(hessian_sparsity) == 0
    end

    @testset "Hessians and Hess-vec" begin
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
        @test dense_hessian(hessian_sparsity, V, 3) ≈ correct_hessian

        h = ones(3) # The input values should be overwritten.
        v = [2.4, 3.5, 1.2]
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0, 3.0])
        @test h ≈ correct_hessian * v
    end

    @testset "Hessians and Hess-vec with subexpressions" begin
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
        @test dense_hessian(hessian_sparsity, V, 3) ≈ correct_hessian

        h = ones(3) # The input values should be overwritten.
        v = [2.4, 3.5, 1.2]
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0, 3.0])
        @test h ≈ correct_hessian * v
    end

    @testset "Jacobians and Jac-vec" begin
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
        @test dense_jacobian(jacobian_sparsity, V, 2, 3) ≈ correct_jacobian

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
    end

    @testset "Jacobians and Jac-vec with subexpressions" begin
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
        @test dense_jacobian(jacobian_sparsity, V, 2, 3) ≈ correct_jacobian

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
    end

    @testset "Expression graphs" begin
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
    end

    @testset "More expression graphs" begin
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
    end

    @testset "Expression graph for ifelse" begin
        m = Model()
        @variable(m, x)
        @NLobjective(m, Min, ifelse(x <= 1, x^2, x))
        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:ExprGraph])
        xidx = x.index
        @test MOI.objective_expr(d) ==
              :(ifelse(x[$xidx] <= 1, x[$xidx]^2, x[$xidx]))
    end

    @testset "Expression graph for empty sum and prod" begin
        m = Model()
        @variable(m, x)
        @NLconstraint(m, x <= sum(0 for i in []) + prod(1 for i in []))
        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:ExprGraph])
        xidx = x.index
        @test MOI.constraint_expr(d, 1) == :(x[$xidx] - (0 + 1) <= 0.0)
    end

    @testset "NLparameter" begin
        model = Model()
        @NLparameter(model, p == 1.0)
        @test JuMP.value(p) == 1.0
    end

    @testset "NLparameter set_value" begin
        model = Model()
        @NLparameter(model, p == 1.0)
        JuMP.set_value(p, 10.0)
        @test JuMP.value(p) == 10.0
    end

    @testset "@NLconstraints" begin
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
    end

    # This covers the code that computes Hessians in odd chunks of Hess-vec
    # products.
    @testset "Dense Hessian" begin
        m = Model()
        @variable(m, x[1:18])
        @NLobjective(m, Min, prod(x[i] for i in 1:18))

        d = JuMP.NLPEvaluator(m)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = ones(18)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 18) ≈
              ones(18, 18) - diagm(0 => ones(18))

        values[1] = 0.5
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test dense_hessian(hessian_sparsity, V, 18) ≈ [
            0 ones(17)'
            ones(17) (ones(17, 17)-diagm(0 => ones(17)))/2
        ]
    end

    @testset "eval_objective and eval_objective_gradient" begin
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
    end

    @testset "eval_constraint and Jacobians" begin
        m = Model()
        @variable(m, x[1:4])
        @NLparameter(m, p == 2)
        @NLexpression(m, ex, p * x[1])

        ψ(x) = sin(x)
        t(x, y) = x + 3y
        JuMP.register(m, :ψ, 1, ψ, autodiff = true)
        JuMP.register(m, :t, 2, t, autodiff = true)

        @NLconstraint(
            m,
            Min,
            ex / 2 + sin(x[2]) / ψ(x[2]) + t(x[3], x[4]) <= 0.0
        )
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
    end

    @testset "add_nonlinear_expression, set_nonlinear_objective, and add_nonlinear_constraint" begin
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
            (
                +(-1.0 * x[$xidx]) +
                x[$xidx] * x[$xidx] +
                x[$yidx] * x[$yidx] * 2.0
            ) - 1.0 == 0.0
        )
    end

    @testset "Test views on Hessian functions" begin
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
    end

    @testset "Constant expressions" begin
        model = Model()
        @variable(model, x)
        @NLexpression(model, expr, 10)
        @NLobjective(model, Min, expr + x)
        d = JuMP.NLPEvaluator(model)
        MOI.initialize(d, [:Grad])
        grad = zeros(1)
        MOI.eval_objective_gradient(d, grad, [2.0])
        @test grad == [1.0]
    end

    @testset "User-defined function with variable closure" begin
        model = Model()
        @variable(model, x[1:2])
        f(x1) = x1 + x[2]
        @test_throws(
            ErrorException(
                "Expected return type of `Float64` from the user-defined " *
                "function :f, but got `AffExpr`.",
            ),
            register(model, :f, 1, f; autodiff = true),
        )
    end

    @testset "User-defined function with variable closure after register" begin
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
    end

    @testset "User-defined function returning bad type after register" begin
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
    end

    @testset "User-defined function returning bad type" begin
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
    end

    @testset "JuMP.value on NonlinearExpressions" begin
        model = Model()
        @variable(model, x)
        @NLexpression(model, ex1, sin(x))
        @NLexpression(model, ex2, ex1 + x^2)
        @NLexpression(model, ex3, 2ex1 + ex2 / 2)
        JuMP.set_start_value(x, 2.0)
        @test JuMP.value(JuMP.start_value, ex1) ≈ sin(2.0)
        @test JuMP.value(JuMP.start_value, ex2) ≈ sin(2.0) + 4.0
        @test JuMP.value(JuMP.start_value, ex3) ≈ 2.5 * sin(2.0) + 2.0
    end

    @testset "Hessians disabled with user-defined multivariate functions" begin
        model = Model()
        my_f(x, y) = (x - 1)^2 + (y - 2)^2
        JuMP.register(model, :my_f, 2, my_f, autodiff = true)
        @variable(model, x[1:2])
        @NLobjective(model, Min, my_f(x[1], x[2]))
        evaluator = JuMP.NLPEvaluator(model)
        @test !(:Hess in MOI.features_available(evaluator))
    end

    @testset "AffExpr in nonlinear" begin
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
    end

    @testset "QuadExpr in nonlinear" begin
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
    end

    @testset "Error on complex values" begin
        model = Model()
        @variable(model, x)
        c = sqrt(Complex(-1))
        expected_exception = ErrorException(
            "Unexpected object $c (of type $(typeof(c)) in nonlinear expression.",
        )
        @test_throws expected_exception @NLobjective(model, Min, c * x)
    end
    @testset "Special functions" begin
        model = Model()
        @variable(model, x)
        @NLconstraint(model, c1, erf(x) <= 0.0)
        d = NLPEvaluator(model)
        MOI.initialize(d, Symbol[:Grad])
        out = zeros(1)
        MOI.eval_constraint(d, out, [2.0])
        @test out[1] ≈ 0.9953222 atol = 1e-6
    end

    @testset "JuMP_extensions" begin
        model = JuMPExtension.MyModel()
        @variable(model, x)
        err = ErrorException(
            "Encountered an error parsing nonlinear expression: we don't support " *
            "models of type $(typeof(model)). In general, JuMP's nonlinear features " *
            "don't work with JuMP-extensions.",
        )
        @test_throws(err, @NLexpression(model, sqrt(x)))
    end

    @testset "rad2deg and deg2rad" begin
        x = 1.0
        id = JuMP._Derivatives.univariate_operator_to_id[:rad2deg]
        y = rad2deg(x)
        @test JuMP._Derivatives.eval_univariate_2nd_deriv(id, x, y) == 0.0
        id = JuMP._Derivatives.univariate_operator_to_id[:deg2rad]
        y = deg2rad(x)
        @test JuMP._Derivatives.eval_univariate_2nd_deriv(id, x, y) == 0.0
    end

    @testset "Re-register univariate" begin
        model = Model()
        @variable(model, x >= 0)
        f(x) = x^2
        register(model, :f, 1, f; autodiff = true)
        @test_throws ErrorException register(model, :f, 1, f; autodiff = true)
    end

    @testset "Re-register multivariate" begin
        model = Model()
        @variable(model, x >= 0)
        f(x, y) = x + y
        register(model, :f, 2, f; autodiff = true)
        @test_throws ErrorException register(model, :f, 2, f; autodiff = true)
    end

    @testset "Check univariate NLconstraint is valid" begin
        model = Model()
        model2 = Model()
        @variable(model, x >= 0)
        c = @NLconstraint(model, exp(x) <= 1)
        @test is_valid(model, c)
        @test !is_valid(model2, c)
    end

    @testset "Check multivariate NLconstraint is valid" begin
        model = Model()
        model2 = Model()
        @variable(model, x >= 0)
        @variable(model, y >= 0)
        c = @NLconstraint(model, exp(x) + log(y) <= 1)
        @test is_valid(model, c)
        @test !is_valid(model2, c)
    end

    @testset "Check multivariate NLconstraint is valid" begin
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
    end

    @testset "Interval errors" begin
        model = Model()
        @variable(model, x)
        err = ErrorException(
            "Interval constraint contains non-constant left- or right-hand " *
            "sides. Reformulate as two separate constraints, or move all " *
            "variables into the central term.",
        )
        @test_throws err add_nonlinear_constraint(model, :($x <= $x <= 2 * $x))
    end

    @testset "dual_start_value" begin
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
        @test_throws(
            ArgumentError,
            set_nonlinear_dual_start_value(model, [1.0]),
        )
    end

    @testset "nonlinear_deprecations" begin
        model = Model()
        @variable(model, x)
        @test_throws(ErrorException, num_nl_constraints(model))
        @test_throws(ErrorException, all_nl_constraints(model))
        @test_throws(ErrorException, add_NL_constraint(model, :(x[$x] <= 1)))
        @test_throws(ErrorException, add_NL_expression(model, :(x[$x] + 1)))
        @test_throws(
            ErrorException,
            set_NL_objective(model, MOI.MIN_SENSE, :(x[$x] + 1)),
        )
    end
end
