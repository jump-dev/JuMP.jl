# For "expression^3 and unary*"
struct PowVariable <: JuMP.AbstractJuMPScalar
    pow::Int
end
Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow*i)
Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)
Base.copy(x::PowVariable) = x

@testset "Expression" begin
    @testset "value for GenericAffExpr" begin
        expr1 = JuMP.GenericAffExpr([3, 2], [-5., 4.], 3.)
        @test @inferred(JuMP.value(expr1, -)) == 10.
        expr2 = JuMP.GenericAffExpr(Int[], Int[], 2)
        @test typeof(@inferred(JuMP.value(expr2, i -> 1.0))) == Float64
        @test @inferred(JuMP.value(expr2, i -> 1.0)) == 2.0
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
