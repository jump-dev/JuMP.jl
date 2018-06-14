function objectives_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "SingleVariable objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x)
        @test JuMP.objectivesense(m) == :Min
        @test JuMP.objectivefunction(m, VariableRefType) == x

        @objective(m, Max, x)
        @test JuMP.objectivesense(m) == :Max
        @test JuMP.objectivefunction(m, VariableRefType) == x
    end

    @testset "Linear objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, 2x)
        @test JuMP.objectivesense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, AffExprType), 2x)

        @objective(m, Max, x + 3x + 1)
        @test JuMP.objectivesense(m) == :Max
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, AffExprType), 4x + 1)
    end

    @testset "Quadratic objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x^2 + 2x)
        @test JuMP.objectivesense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objectivefunction(m, QuadExprType), x^2 + 2x)
        @test_throws ErrorException JuMP.objectivefunction(m, AffExprType)
    end
end


@testset "Objectives for JuMP.Model" begin
    objectives_test(Model, VariableRef{Model{JuMP.NonDirectBackendType}})
end

@testset "Objectives for JuMPExtension.MyModel" begin
    objectives_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
