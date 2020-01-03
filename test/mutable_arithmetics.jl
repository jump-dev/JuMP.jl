using LinearAlgebra, Test

import MutableArithmetics
const MA = MutableArithmetics

using JuMP

function promote_operation_test(op::Function, x::Type, y::Type)
    f() = JuMP._MA.promote_operation(op, x, y)
    @test typeof(op(zero(x), zero(y))) == f()
    @test 0 == @allocated f()
end

function mutable_arithmetics_test(ModelType::Type{<:JuMP.AbstractModel},
                                  VariableRefType::Type{<:JuMP.AbstractVariableRef})
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}

    @testset "promote_operation" begin
        for op in [+, -, *]
            for T in [Int, Float64]
                promote_operation_test(op, T, VariableRefType)
                promote_operation_test(op, VariableRefType, T)
                promote_operation_test(op, T, AffExprType)
                promote_operation_test(op, AffExprType, T)
                promote_operation_test(op, T, QuadExprType)
                promote_operation_test(op, QuadExprType, T)
            end
            promote_operation_test(op, VariableRefType, VariableRefType)
            promote_operation_test(op, VariableRefType, AffExprType)
            promote_operation_test(op, AffExprType, VariableRefType)
            if op != *
                promote_operation_test(op, VariableRefType, QuadExprType)
                promote_operation_test(op, QuadExprType, VariableRefType)
                promote_operation_test(op, AffExprType, QuadExprType)
                promote_operation_test(op, QuadExprType, AffExprType)
            end
        end
    end

    @testset "Int" begin
        model = ModelType()
        @variable(model, x)
        @testset "Affine" begin
            MA.Test.int_test(typeof(1x), exclude = ["int_mul", "int_add", "int_add_mul"])
        end
        @testset "Quadratic" begin
            MA.Test.int_test(typeof(1x^2), exclude = ["int_mul", "int_add", "int_add_mul"])
        end
    end

    @testset "Scalar" begin
        model = ModelType()
        @variable(model, x)
        exclude = ["cube"]
        MA.Test.scalar_test(x, exclude = exclude)
        MA.Test.scalar_test(2x + 3, exclude = exclude)
        MA.Test.scalar_test(2x^2 + 4x + 1, exclude = exclude)
    end
    @testset "Quadratic" begin
        model = ModelType()
        @variable(model, w)
        @variable(model, x)
        @variable(model, y)
        @variable(model, z)
        MA.Test.quadratic_test(w, x, y, z)
    end
    @testset "Sparse" begin
        model = ModelType()
        @variable(model, X11)
        @variable(model, X23)
        @variable(model, Xd[1:3, 1:3])
        MA.Test.sparse_test(X11, X23, Xd)
    end
    @testset "Vector" begin
        model = ModelType()
        @variable(model, x[1:3])
        MA.Test.array_test(x)
    end
    @testset "Matrix" begin
        model = ModelType()
        @variable(model, x[1:2, 1:2])
        MA.Test.array_test(x)
        @variable(model, y[1:2, 1:2], Symmetric)
        MA.Test.array_test(y)
        @variable(model, z[1:2, 1:3])
        MA.Test.array_test(z)
    end
    @testset "DenseAxisVector" begin
        model = ModelType()
        @variable(model, y[2:5])
        MA.Test.array_test(y, exclude = ["matrix_vector", "non_array"])
    end

end

@testset "Operators for JuMP.Model" begin
    mutable_arithmetics_test(Model, VariableRef)
end

@testset "Operators for JuMPExtension.MyModel" begin
    mutable_arithmetics_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end
