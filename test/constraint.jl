

@testset "Constraints" begin
    @testset "Linear constraints" begin
        m = Model()
        @variable(m, x)

        cref = @constraint(m, 2x <= 10)
        c = JuMP.getconstraint(cref, AffExpr)
        @test JuMP.isequal_canonical(c.func, 2x)
        @test c.set == MOI.LessThan(10.0)
        
        cref = @constraint(m, 3x + 1 ≥ 10)
        c = JuMP.getconstraint(cref, AffExpr)
        @test JuMP.isequal_canonical(c.func, 3x)
        @test c.set == MOI.GreaterThan(9.0)
        
        cref = @constraint(m, 1 == -x)
        c = JuMP.getconstraint(cref, AffExpr)
        @test JuMP.isequal_canonical(c.func, 1.0x)
        @test c.set == MOI.EqualTo(-1.0)

    end
end
