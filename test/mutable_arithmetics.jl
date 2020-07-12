module TestMutableArithmetics

using LinearAlgebra
using JuMP
using Test

const MA = JuMP._MA

# @static if !(:JuMPExtension in names(Main))
#     include(joinpath(@__DIR__, "JuMPExtension.jl"))
# end

struct DummyVariableRef <: JuMP.AbstractVariableRef end
JuMP.name(::DummyVariableRef) = "dummy"
struct DummyExpr end
Base.:(+)(::GenericAffExpr{Float64,DummyVariableRef}, ::JuMP.AbstractJuMPScalar) = DummyExpr()
Base.:(+)(::JuMP.AbstractJuMPScalar, ::JuMP.GenericAffExpr{Float64,DummyVariableRef}) = DummyExpr()
Base.:(-)(::GenericAffExpr{Float64,DummyVariableRef}, ::JuMP.AbstractJuMPScalar) = DummyExpr()
Base.:(-)(::JuMP.AbstractJuMPScalar, ::JuMP.GenericAffExpr{Float64,DummyVariableRef}) = DummyExpr()

function promote_operation_test(op::Function, x::Type, y::Type)
    f() = JuMP._MA.promote_operation(op, x, y)
    @test typeof(op(zero(x), zero(y))) == f()
    @test 0 == @allocated f()
end

function test_promote_operation(ModelType, VariableRefType)
    AffExprType = JuMP.GenericAffExpr{Float64, VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64, VariableRefType}
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

function test_int(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    MA.Test.int_test(typeof(1x), exclude = ["int_mul", "int_add", "int_add_mul"])
    MA.Test.int_test(typeof(1x^2), exclude = ["int_mul", "int_add", "int_add_mul"])
end

function test_scalar(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    exclude = ["cube"]
    MA.Test.scalar_test(x, exclude = exclude)
    MA.Test.scalar_test(2x + 3, exclude = exclude)
    MA.Test.scalar_test(2x^2 + 4x + 1, exclude = exclude)
end

function test_quadratic(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    MA.Test.quadratic_test(w, x, y, z)
end

function test_sparse(ModelType, ::Any)
    model = ModelType()
    @variable(model, X11)
    @variable(model, X23)
    @variable(model, Xd[1:3, 1:3])
    MA.Test.sparse_test(X11, X23, Xd)
end

function test_vector(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:3])
    MA.Test.array_test(x)
end

function test_symmetric_matrix(ModelType, ::Any)
    model = ModelType()
    @variable(model, y[1:2, 1:2], Symmetric)
    MA.Test.array_test(y)
end

function test_nonsquare_matrix(ModelType, ::Any)
    model = ModelType()
    @variable(model, z[1:2, 1:3])
    MA.Test.array_test(z)
end

function test_DenseAxisVector(ModelType, ::Any)
    model = ModelType()
    @variable(model, y[2:5])
    MA.Test.array_test(y, exclude = ["matrix_vector", "non_array"])
end

function test_different_variables(ModelType, ::Any)
    model = ModelType()
    x = @variable(model)
    y = DummyVariableRef()
    aff = x + 1
    function _promote_test(a, b)
        A = typeof(a)
        B = typeof(b)
        @test MA.promote_operation(+, A, B) == DummyExpr
        @test MA.promote_operation(-, A, B) == DummyExpr
        @test MA.promote_operation(+, B, A) == DummyExpr
        @test MA.promote_operation(-, B, A) == DummyExpr
    end
    _promote_test(x, y)
    _promote_test(aff, y)
end

function runtests()
    @testset "mutable_arithmetics.jl" begin
        for name in names(@__MODULE__; all = true)
            if !startswith("$(name)", "test_")
                continue
            end
            f = getfield(@__MODULE__, name)
            @testset "$(name)" begin
                f(Model, VariableRef)
            end
            # Note(odow): We disable JuMPExtension tests for MutableArithmetics
            # because they take far to looooooong. MutableArithmetics is a
            # tested package. We're also testing that it works with JuMP. We
            # don't need to double up on our tests.
            # @testset "$(name)-JuMPExtension" begin
            #     f(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
            # end
        end
    end
end

end

TestMutableArithmetics.runtests()
