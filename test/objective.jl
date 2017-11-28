

@testset "Objectives" begin
    @testset "Linear objectives" begin
        m = Model()
        @variable(m, x)

        @objective(m, Min, 2x)
        @test JuMP.objectivesense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, AffExpr), 2x)

        @objective(m, Max, x + 3x + 1)
        @test JuMP.objectivesense(m) == :Max
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, AffExpr), 4x + 1)
    end

    @testset "Quadratic objectives" begin
        m = Model()
        @variable(m, x)

        @objective(m, Min, x^2 + 2x)
        @test JuMP.objectivesense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, QuadExpr), x^2 + 2x)
        @test_throws TypeError JuMP.objectivefunction(m, AffExpr)
    end
end
