

@testset "Constraints" begin
    @testset "AffExpr constraints" begin
        m = Model()
        @variable(m, x)

        cref = @constraint(m, 2x <= 10)
        c = JuMP.constraintobject(cref, AffExpr, MOI.LessThan)
        @test JuMP.isequal_canonical(c.func, 2x)
        @test c.set == MOI.LessThan(10.0)
        @test_throws TypeError JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        cref = @constraint(m, 3x + 1 ≥ 10)
        c = JuMP.constraintobject(cref, AffExpr, MOI.GreaterThan)
        @test JuMP.isequal_canonical(c.func, 3x)
        @test c.set == MOI.GreaterThan(9.0)

        cref = @constraint(m, 1 == -x)
        c = JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, 1.0x)
        @test c.set == MOI.EqualTo(-1.0)

    end

    @testset "QuadExpr constraints" begin
        m = Model()
        @variable(m, x)
        @variable(m, y)

        cref = @constraint(m, x^2 + x <= 1)
        c = JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test JuMP.isequal_canonical(c.func, x^2 + x)
        @test c.set == MOI.LessThan(1.0)

        cref = @constraint(m, y*x - 1.0 == 0.0)
        c = JuMP.constraintobject(cref, QuadExpr, MOI.EqualTo)
        @test JuMP.isequal_canonical(c.func, x*y)
        @test c.set == MOI.EqualTo(1.0)
        @test_throws TypeError JuMP.constraintobject(cref, QuadExpr, MOI.LessThan)
        @test_throws TypeError JuMP.constraintobject(cref, AffExpr, MOI.EqualTo)

        # cref = @constraint(m, [x^2 - 1] in MOI.SecondOrderCone(1))
        # c = JuMP.constraintobject(cref, QuadExpr, MOI.SecondOrderCone)
        # @test JuMP.isequal_canonical(c.func, -1 + x^2)
        # @test c.set == MOI.SecondOrderCone(1)
    end
end
