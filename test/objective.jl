

@testset "Objectives" begin
    @testset "Linear objectives" begin
        m = Model()
        @variable(m, x)

        @objective(m, Min, 2x)
        @test getobjectivesense(m) == :Min
        @test JuMP.isequal_canonical(getobjective(m, AffExpr), 2x)

        @objective(m, Max, x + 3x + 1)
        @test getobjectivesense(m) == :Max
        @test JuMP.isequal_canonical(getobjective(m, AffExpr), 4x + 1)
    end
end
