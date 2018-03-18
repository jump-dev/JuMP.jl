
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

    @testset "Hessian evaluation (Issue #435)" begin
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)

        @NLexpression(m, foo, a * b + c^2)

        @NLobjective(m, Min, foo)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        I,J = MOI.hessian_lagrangian_structure(d)
        V = zeros(length(I))
        values = [1.0, 2.0, 3.0] # Values for a, b, and c, respectively.
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        hess_raw = sparse(I,J,V)
        # Convert from lower triangular
        hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 2.0])

        # make sure we don't get NaNs in this case
        @NLobjective(m, Min, a * b + 3*c^2)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:Hess])
        values = [1.0, 2.0, -1.0]
        V = zeros(length(I))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        hess_raw = sparse(I,J,V)
        hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])

        # Initialize again
        MOI.initialize!(d, [:Hess])
        V = zeros(length(I))
        MOI.eval_hessian_lagrangian(d, V, values, 1.0, Float64[])
        hess_raw = sparse(I,J,V)
        hess_sparse = hess_raw + hess_raw' - sparse(diagm(diag(hess_raw)))
        @test isapprox(hess_sparse, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 6.0])
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

    @testset "Hess-vec" begin
        m = Model()
        @variable(m, a)
        @variable(m, b)
        @variable(m, c)

        @NLobjective(m, Min, a*b + c^2)
        @NLconstraint(m, c*b <= 1)
        @NLconstraint(m, a^2/2 <= 1)
        d = JuMP.NLPEvaluator(m)
        MOI.initialize!(d, [:HessVec])
        h = ones(3) # test that input values are overwritten
        v = [2.4,3.5,1.2]
        values = [1.0, 2.0, 3.0] # For a, b, c.
        MOI.eval_hessian_lagrangian_product(d, h, values, v, 1.0, [2.0,3.0])
        correct = [3.0 1.0 0.0; 1.0 0.0 2.0; 0.0 2.0 2.0]*v
        @test isapprox(h, correct)
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

end
