#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestMutableArithmetics

using LinearAlgebra
using JuMP
using Test

const MA = JuMP._MA

include(joinpath(@__DIR__, "JuMPExtension.jl"))

struct DummyVariableRef <: JuMP.AbstractVariableRef end
JuMP.name(::DummyVariableRef) = "dummy"
struct DummyExpr end
function Base.:(+)(
    ::GenericAffExpr{Float64,DummyVariableRef},
    ::JuMP.AbstractJuMPScalar,
)
    return DummyExpr()
end
function Base.:(+)(
    ::JuMP.AbstractJuMPScalar,
    ::JuMP.GenericAffExpr{Float64,DummyVariableRef},
)
    return DummyExpr()
end
function Base.:(-)(
    ::GenericAffExpr{Float64,DummyVariableRef},
    ::JuMP.AbstractJuMPScalar,
)
    return DummyExpr()
end
function Base.:(-)(
    ::JuMP.AbstractJuMPScalar,
    ::JuMP.GenericAffExpr{Float64,DummyVariableRef},
)
    return DummyExpr()
end

function promote_operation_test(op::Function, x::Type, y::Type)
    f() = JuMP._MA.promote_operation(op, x, y)
    @test typeof(op(zero(x), zero(y))) == f()
    @test 0 == @allocated f()
end

function test_promote_operation(ModelType, VariableRefType)
    for S in [Float64, ComplexF64]
        AffExprType = JuMP.GenericAffExpr{S,VariableRefType}
        QuadExprType = JuMP.GenericQuadExpr{S,VariableRefType}
        for op in [+, -, *]
            for T in [Int, S]
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
end

function test_int(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    for a in [1, 1im]
        MA.Test.int_test(
            typeof(a * x);
            exclude = ["int_mul", "int_add", "int_add_mul"],
        )
        MA.Test.int_test(
            typeof(a * x^2);
            exclude = ["int_mul", "int_add", "int_add_mul"],
        )
    end
    return
end

function test_scalar(ModelType, ::Any)
    model = ModelType()
    @variable(model, x)
    exclude = ["cube"]
    MA.Test.scalar_test(x; exclude = exclude)
    for a in [2, 2im]
        MA.Test.scalar_test(a * x + 3; exclude = exclude)
        MA.Test.scalar_test(a * x^2 + 4x + 1; exclude = exclude)
    end
    return
end

function test_quadratic(ModelType, ::Any)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    return MA.Test.quadratic_test(w, x, y, z)
end

function test_sparse(ModelType, ::Any)
    model = ModelType()
    @variable(model, X11)
    @variable(model, X23)
    @variable(model, Xd[1:3, 1:3])
    return MA.Test.sparse_test(X11, X23, Xd)
end

function test_vector(ModelType, ::Any)
    model = ModelType()
    @variable(model, x[1:3])
    return MA.Test.array_test(x)
end

function test_symmetric_matrix(ModelType, ::Any)
    model = ModelType()
    @variable(model, y[1:2, 1:2], Symmetric)
    return MA.Test.array_test(y)
end

function test_nonsquare_matrix(ModelType, ::Any)
    model = ModelType()
    @variable(model, z[1:2, 1:3])
    return MA.Test.array_test(z)
end

function test_DenseAxisVector(ModelType, ::Any)
    model = ModelType()
    @variable(model, y[2:5])
    return MA.Test.array_test(y; exclude = ["matrix_vector", "non_array"])
end

function test_different_variables(ModelType, ::Any)
    model = ModelType()
    x = @variable(model)
    y = DummyVariableRef()
    function _promote_test(a, b)
        A = typeof(a)
        B = typeof(b)
        @test MA.promote_operation(+, A, B) == DummyExpr
        @test MA.promote_operation(-, A, B) == DummyExpr
        @test MA.promote_operation(+, B, A) == DummyExpr
        @test MA.promote_operation(-, B, A) == DummyExpr
    end
    _promote_test(x, y)
    for aff in [x + 1, x + im]
        _promote_test(aff, y)
    end
    return
end

function runtests()
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
    @testset "test_promote_operation-JuMPExtension" begin
        test_promote_operation(
            JuMPExtension.MyModel,
            JuMPExtension.MyVariableRef,
        )
    end
end

end

TestMutableArithmetics.runtests()
