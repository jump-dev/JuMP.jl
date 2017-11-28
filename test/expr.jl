@testset "value for GenericAffExpr" begin
    expr1 = JuMP.GenericAffExpr([3, 2], [-5., 4.], 3.)
    @test @inferred(JuMP.value(expr1, -)) == 10.
    expr2 = JuMP.GenericAffExpr(Int[], Int[], 2)
    @test typeof(@inferred(JuMP.value(expr2, i -> 1.0))) == Float64
    @test @inferred(JuMP.value(expr2, i -> 1.0)) == 2.0
end
