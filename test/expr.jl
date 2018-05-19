# For "expression^3 and unary*"
struct PowVariable <: JuMP.AbstractVariableRef
    pow::Int
end
Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow*i)
Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)
Base.copy(x::PowVariable) = x

@testset "Expression" begin
    @testset "value for GenericAffExpr" begin
        expr1 = JuMP.GenericAffExpr(3.0, 3 => -5.0, 2 => 4.0)
        @test @inferred(JuMP.value(expr1, -)) == 10.0
        expr2 = JuMP.GenericAffExpr{Int,Int}(2)
        @test typeof(@inferred(JuMP.value(expr2, i -> 1.0))) == Float64
        @test @inferred(JuMP.value(expr2, i -> 1.0)) == 2.0
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::V)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, :b),
                                     JuMP.GenericAffExpr(1.0, :a => 2.0, :b => 1.0))
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::C)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(JuMP.add_to_expression!(aff, 1.0),
                                     JuMP.GenericAffExpr(2.0, :a => 2.0))
    end

    @testset "linearterms(::AffExpr)" begin
        m = Model()
        @variable(m, x[1:10])

        aff = 1*x[1] + 2*x[2]
        k = 0
        @test length(linearterms(aff)) == 2
        for (coeff, var) in linearterms(aff)
            if k == 0
                @test coeff == 1
                @test var === x[1]
            elseif k == 1
                @test coeff == 2
                @test var === x[2]
            end
            k += 1
        end
        @test k == 2
    end

    @testset "linearterms(::AffExpr) for empty expression" begin
        k = 0
        aff = zero(AffExpr)
        @test length(linearterms(aff)) == 0
        for (coeff, var) in linearterms(aff)
            k += 1
        end
        @test k == 0
    end

    @testset "Copy AffExpr between models" begin
        m = Model()
        @variable(m, x)
        m2 = Model()
        aff = copy(2x + 1, m2)
        aff_expected = 2*copy(x, m2) + 1
        @test JuMP.isequal_canonical(aff, aff_expected)
    end

    @testset "expression^3 and unary*" begin
        m = Model()
        x = PowVariable(1)
        # Calls (*)((x*x)^6)
        y = @expression m (x*x)^3
        @test y.pow == 6
        z = @inferred (x*x)^3
        @test z.pow == 6
    end
end
