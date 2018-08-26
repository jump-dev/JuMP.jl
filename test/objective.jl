function objectives_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "SingleVariable objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x)
        @test JuMP.objective_sense(m) == :Min
        @test JuMP.objective_function(m, VariableRefType) == x

        @objective(m, Max, x)
        @test JuMP.objective_sense(m) == :Max
        @test JuMP.objective_function(m, VariableRefType) == x
    end

    @testset "Linear objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, 2x)
        @test JuMP.objective_sense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 2x)

        @objective(m, Max, x + 3x + 1)
        @test JuMP.objective_sense(m) == :Max
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 4x + 1)
    end

    @testset "Quadratic objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x^2 + 2x)
        @test JuMP.objective_sense(m) == :Min
        @test JuMP.isequal_canonical(JuMP.objective_function(m, QuadExprType), x^2 + 2x)
        @test_throws ErrorException JuMP.objective_function(m, AffExprType)
    end
end


@testset "Objectives for JuMP.Model" begin
    objectives_test(Model, VariableRef)
end

@testset "Objectives for JuMPExtension.MyModel" begin
    objectives_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
