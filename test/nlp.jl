
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

end
