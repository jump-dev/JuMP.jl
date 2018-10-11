function objectives_test(ModelType::Type{<:JuMP.AbstractModel}, VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "objective_sense set and get" begin
        model = ModelType()
        JuMP.set_objective_sense(model, MOI.FeasibilitySense)
        @test JuMP.objective_sense(model) == MOI.FeasibilitySense
    end

    @testset "SingleVariable objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x)
        @test JuMP.objective_sense(m) == MOI.MinSense
        @test JuMP.objective_function(m, VariableRefType) == x

        @objective(m, Max, x)
        @test JuMP.objective_sense(m) == MOI.MaxSense
        @test JuMP.objective_function(m, VariableRefType) == x
    end

    @testset "Linear objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, 2x)
        @test JuMP.objective_sense(m) == MOI.MinSense
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 2x)

        @objective(m, Max, x + 3x + 1)
        @test JuMP.objective_sense(m) == MOI.MaxSense
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 4x + 1)
    end

    @testset "Quadratic objectives" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, Min, x^2 + 2x)
        @test JuMP.objective_sense(m) == MOI.MinSense
        @test JuMP.isequal_canonical(JuMP.objective_function(m, QuadExprType), x^2 + 2x)
        @test_throws InexactError JuMP.objective_function(m, AffExprType)
    end

    @testset "Sense as symbol" begin
        m = ModelType()
        @variable(m, x)

        @objective(m, :Min, 2x)
        @test JuMP.objective_sense(m) == MOI.MinSense
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 2x)
    end

    @testset "Sense in variable" begin
        m = ModelType()
        @variable(m, x)

        sense = :Min
        @objective(m, sense, 2x)
        @test JuMP.objective_sense(m) == MOI.MinSense
        @test JuMP.isequal_canonical(JuMP.objective_function(m, AffExprType), 2x)

        sense = :Man
        @test_throws ErrorException @objective(m, sense, 2x)
    end
end


@testset "Objectives for JuMP.Model" begin
    objectives_test(Model, VariableRef)
end

@testset "Objectives for JuMPExtension.MyModel" begin
    objectives_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
