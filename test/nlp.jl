# TODO: Replace isapprox with ≈ everywhere.
@testset "Nonlinear" begin

    import JuMP: NonlinearExprData

    function expressions_equal(ex1::NonlinearExprData, ex2::NonlinearExprData)
        return ex1.nd == ex2.nd && ex1.const_values == ex2.const_values
    end

    # TODO: These are poorly designed tests because they would not catch errors
    # that affect both processNLExpr and NonlinearExprData's parsing (although
    # they use different pathways). It would be better to check the NodeData
    # representation directly.

    @testset "Parse + (binary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @test expressions_equal(@JuMP.processNLExpr(m, x + y),
                                NonlinearExprData(m, :($x + $y)))
    end

    @testset "Parse + (ternary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        @test expressions_equal(@JuMP.processNLExpr(m, x + y + z),
                                NonlinearExprData(m, :($x + $y + $z)))
    end

    @testset "Parse * (binary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @test expressions_equal(@JuMP.processNLExpr(m, x * y),
                                NonlinearExprData(m, :($x * $y)))
    end

    @testset "Parse * (ternary)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        @test expressions_equal(@JuMP.processNLExpr(m, x * y * z),
                                NonlinearExprData(m, :($x * $y * $z)))
    end

    @testset "Parse ^ (binary)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(@JuMP.processNLExpr(m, x^3),
                                NonlinearExprData(m, :($x^3)))
    end

    @testset "Parse sin (univariate)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(@JuMP.processNLExpr(m, sin(x)),
                                NonlinearExprData(m, :(sin($x))))
    end

    @testset "Parse ifelse" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(@JuMP.processNLExpr(m, ifelse(1 == 2 || 3 == 4 && 5 == 6, x, 0.0)),
                                NonlinearExprData(m, :(ifelse(1 == 2 || 3 == 4 && 5 == 6, $x, 0.0))))
    end

    @testset "Parse ifelse (3-way comparison)" begin
        m = Model()
        @variable(m, x)
        @test expressions_equal(@JuMP.processNLExpr(m, ifelse(1 <= 2 <= 3, x, 0.0)),
                                NonlinearExprData(m, :(ifelse(1 <= 2 <= 3, $x, 0.0))))
    end

    @testset "Parse sum" begin
        m = Model()
        @variable(m, x[1:2])
        @test expressions_equal(@JuMP.processNLExpr(m, sum(x[i] for i in 1:2)),
                                NonlinearExprData(m, :($(x[1]) + $(x[2]))))
    end

    @testset "Parse prod" begin
        m = Model()
        @variable(m, x[1:2])
        @test expressions_equal(@JuMP.processNLExpr(m, prod(x[i] for i in 1:2)),
                                NonlinearExprData(m, :($(x[1]) * $(x[2]))))
    end

    @testset "Parse subexpressions" begin
        m = Model()
        @variable(m, x)
        @NLexpression(m, ex, x^2)
        @test expressions_equal(@JuMP.processNLExpr(m, ex + 1),
                                NonlinearExprData(m, :($ex + 1)))
    end

    @testset "Parse subexpressions" begin
        m = Model()
        @NLparameter(m, param == 10)
        @test expressions_equal(@JuMP.processNLExpr(m, param + 1),
                                NonlinearExprData(m, :($param + 1)))
    end

    @testset "Parse user-defined function (univariate)" begin
        m = Model()
        @variable(m, x)
        f(x) = x
        JuMP.register(m, :f, 1, f, autodiff=true)
        @test expressions_equal(@JuMP.processNLExpr(m, f(x)),
                                NonlinearExprData(m, :(f($x))))
    end

    @testset "Parse user-defined function (multivariate)" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)
        f(x,y) = x
        JuMP.register(m, :f, 2, f, autodiff=true)
        @test expressions_equal(@JuMP.processNLExpr(m, f(x,y)),
                                NonlinearExprData(m, :(f($x,$y))))
    end

    @testset "Error on sum(x)" begin
        m = Model()
        x = [1,2,3]
        @test_throws ErrorException @NLexpression(m, sum(x))
    end

    @testset "Error on non-scalar expressions" begin
        m = Model()
        x = [1,2,3]
        @test_throws ErrorException @NLexpression(m, x + 1)
    end

    # Converts the lower-triangular sparse Hessian in MOI format into a dense
    # matrix.
    function dense_hessian(hessian_sparsity, V, n)
        I = [i for (i,j) in hessian_sparsity]
        J = [j for (i,j) in hessian_sparsity]
        raw = sparse(I, J, V, n, n)
        return Matrix(raw + raw' - sparse(diagm(0=>diag(raw))))
    end

    @testset "Hessian evaluation (Issue #435)" begin
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)

        @NLexpression(m, foo, a * b + c^2)

        @NLobjective(m, Min, foo)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = [1.0, 2.0, 3.0] # Values for a, b, and c, respectively.
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 3), [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0])

        # make sure we don't get NaNs in this case
        @NLobjective(m, Min, a * b + 3*c^2)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        values = [1.0, 2.0, -1.0]
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 3), [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])

        # Initialize again
        MOI.initialize!(d, [:Hess])
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 3), [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])
    end

    @testset "NaN corner case (Issue #695)" begin
        m = Model()
        x0 = 0.0
        y0 = 0.0
        @variable(m, x)
        @variable(m, y)

        @NLobjective(m, Min, (x - x0) /(sqrt(y0) + sqrt(y)))

        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:HessVec])
        h = ones(2)
        v = [2.4,3.5]
        values = [1.0, 2.0] # For x and y.
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, Float64[])
        correct = [0.0 -1/(2*2^(3/2)); -1/(2*2^(3/2)) 3/(4*2^(5/2))]*v
        @test isapprox(h, correct)
    end

    @testset "NaN corner case (Issue #1205)" begin
        m = Model()
        @variable(m, x)

        @NLobjective(m, Min, x^1.0)

        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = zeros(1)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 1), [0.0])
    end

    @testset "NaN corner case - ifelse (Issue #1205)" begin
        m = Model()
        @variable(m, x)

        @NLobjective(m, Min, ifelse(true, x, x^1.0))

        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = zeros(1)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 1), [0.0])
    end

    @testset "Hessians and Hess-vec" begin
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)

        @NLobjective(m, Min, a*b + c^2)
        @NLconstraint(m, c*b <= 1)
        @NLconstraint(m, a^2/2 <= 1)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:HessVec, :Hess])

        values = [1.0, 2.0, 3.0] # For a, b, c.
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, [2.0, 3.0])
        correct_hessian = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]
        @test isapprox(dense_hessian(hessian_sparsity, V, 3), correct_hessian)

        h = ones(3) # The input values should be overwritten.
        v = [2.4,3.5,1.2]
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0,3.0])
        @test isapprox(h, correct_hessian*v)
    end

    @testset "Hessians and Hess-vec with subexpressions" begin
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)

        @NLexpression(m, ab, a*b)
        @NLobjective(m, Min, ab + c^2)
        @NLconstraint(m, c*b <= 1)
        @NLconstraint(m, a^2/2 <= 1)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:HessVec, :Hess])

        values = [1.0, 2.0, 3.0] # For a, b, c.
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, [2.0, 3.0])
        correct_hessian = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]
        @test dense_hessian(hessian_sparsity, V, 3) ≈ correct_hessian

        h = ones(3) # The input values should be overwritten.
        v = [2.4,3.5,1.2]
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0,3.0])
        @test h ≈ correct_hessian*v
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
        MOI.initialize!(d, [:ExprGraph])
        xidx = x.index
        yidx = y.index
        @test MOI.objective_expr(d) == :(x[$xidx]^2.0 + x[$yidx]^2.0)
        @test MOI.constraint_expr(d,1) == :((exp(x[$xidx]) - x[$yidx]) - 0.0 == 0.0)
        @test MOI.constraint_expr(d,2) == :((exp(x[$xidx]) + 1) - 0.0 == 0.0)
    end

    @testset "More expression graphs" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        ψ(x) = 1
        t(x,y) = 2
        JuMP.register(m, :ψ, 1, ψ, autodiff=true)
        JuMP.register(m, :t, 2, t, autodiff=true)

        @NLobjective(m, Min, x^y)
        @NLconstraint(m, sin(x)*cos(y) == 5)
        @NLconstraint(m, nlconstr[i=1:2], i*x^2 == i)
        @NLconstraint(m, -0.5 <= sin(x) <= 0.5)
        @NLconstraint(m, ψ(x) + t(x,y) <= 3)

        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:ExprGraph])
        xidx = x.index
        yidx = y.index
        @test MOI.objective_expr(d) == :(x[$xidx]^x[$yidx])
        @test MOI.constraint_expr(d,1) == :(sin(x[$xidx]) * cos(x[$yidx]) - 5 == 0.0)
        @test MOI.constraint_expr(d,2) == :(1.0*x[$xidx]^2 - 1.0 == 0.0)
        @test MOI.constraint_expr(d,3) == :(2.0*x[$xidx]^2 - 2.0 == 0.0)
        @test MOI.constraint_expr(d,4) == :(-0.5 <= sin(x[$xidx]) <= 0.5)
        @test MOI.constraint_expr(d,5) == :(ψ(x[$xidx]) + t(x[$xidx],x[$yidx]) - 3.0 <= 0.0)
    end

    @testset "Expression graph for ifelse" begin
        m = Model()
        @variable(m, x)
        @NLobjective(m, Min, ifelse( x <= 1, x^2, x) )
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:ExprGraph])
        xidx = x.index
        @test MOI.objective_expr(d) == :(ifelse( x[$xidx] <= 1, x[$xidx]^2, x[$xidx]))
    end

    @testset "Expression graph for empty sum and prod" begin
        m = Model()
        @variable(m, x)
        @NLconstraint(m, x <= sum(0 for i in []) + prod(1 for i in []))
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:ExprGraph])
        xidx = x.index
        @test MOI.constraint_expr(d,1) == :(x[$xidx] - (0 + 1) <= 0.0)
    end

    # This covers the code that computes Hessians in odd chunks of Hess-vec
    # products.
    @testset "Dense Hessian" begin
        m = Model()
        @variable(m, x[1:18])
        @NLobjective(m, Min, prod(x[i] for i=1:18))

        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(hessian_sparsity))
        values = ones(18)
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 18), ones(18,18)-diagm(0=>ones(18)))

        values[1] = 0.5
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        @test isapprox(dense_hessian(hessian_sparsity, V, 18),
                       [0 ones(17)'
                        ones(17)  (ones(17,17) - diagm(0=>ones(17)))/2 ])
    end

    @testset "eval_objective and eval_objective_gradient" begin
        m = Model()
        @variable(m, x[1:4])
        @NLparameter(m, p == 2)
        @NLexpression(m, ex, p*x[1])

        ψ(x) = sin(x)
        t(x,y) = x+3y
        JuMP.register(m, :ψ, 1, ψ, autodiff=true)
        JuMP.register(m, :t, 2, t, autodiff=true)

        @NLobjective(m, Min, ex/2 + sin(x[2])/ψ(x[2]) + t(x[3],x[4]))
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Grad])
        variable_values = fill(2.0, (4,))
        @test isapprox(MOI.eval_objective(d, variable_values), variable_values[1] + 1 + variable_values[3] + 3variable_values[4])
        grad = zeros(4)
        MOI.eval_objective_gradient(d, grad, variable_values)
        @test isapprox(grad, [1.0, 0.0, 1.0, 3.0])
    end

    @testset "eval_constraint and Jacobians" begin
        m = Model()
        @variable(m, x[1:4])
        @NLparameter(m, p == 2)
        @NLexpression(m, ex, p*x[1])

        ψ(x) = sin(x)
        t(x,y) = x+3y
        JuMP.register(m, :ψ, 1, ψ, autodiff=true)
        JuMP.register(m, :t, 2, t, autodiff=true)

        @NLconstraint(m, Min, ex/2 + sin(x[2])/ψ(x[2]) + t(x[3],x[4]) <= 0.0)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Jac])
        variable_values = fill(2.0, (4,))
        constraint_value = zeros(1)
        MOI.eval_constraint(d, constraint_value, variable_values)
        @test isapprox(constraint_value[1], variable_values[1] + 1 + variable_values[3] + 3variable_values[4])
        jacobian_sparsity = MOI.jacobian_structure(d)
        I = [i for (i,j) in jacobian_sparsity]
        J = [j for (i,j) in jacobian_sparsity]
        @test all(I .== 1)
        jac_nonzeros = zeros(length(J))
        MOI.eval_constraint_jacobian(d, jac_nonzeros, variable_values)
        jac_values = zeros(4)
        jac_values[J] = jac_nonzeros
        @test isapprox(jac_values, [1.0, 0.0, 1.0, 3.0])
    end


end
