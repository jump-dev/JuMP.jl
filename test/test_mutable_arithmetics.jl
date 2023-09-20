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

using JuMP
using Test

struct DummyVariableRef <: AbstractVariableRef end

JuMP.name(::DummyVariableRef) = "dummy"

struct DummyExpr end

function Base.:(+)(
    ::GenericAffExpr{Float64,DummyVariableRef},
    ::AbstractJuMPScalar,
)
    return DummyExpr()
end

function Base.:(+)(
    ::AbstractJuMPScalar,
    ::GenericAffExpr{Float64,DummyVariableRef},
)
    return DummyExpr()
end

function Base.:(-)(
    ::GenericAffExpr{Float64,DummyVariableRef},
    ::AbstractJuMPScalar,
)
    return DummyExpr()
end

function Base.:(-)(
    ::AbstractJuMPScalar,
    ::GenericAffExpr{Float64,DummyVariableRef},
)
    return DummyExpr()
end

function promote_operation_test(op::Function, x::Type, y::Type)
    f() = JuMP._MA.promote_operation(op, x, y)
    @test typeof(op(zero(x), zero(y))) == f()
    @test 0 == @allocated f()
    return
end

function test_extension_promote_operation(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    for S in [Float64, ComplexF64]
        AffExprType = GenericAffExpr{S,VariableRefType}
        QuadExprType = GenericQuadExpr{S,VariableRefType}
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
    return
end

function test_extension_int(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x)
    for a in [1, 1im]
        JuMP._MA.Test.int_test(
            typeof(a * x);
            exclude = ["int_mul", "int_add", "int_add_mul"],
        )
        JuMP._MA.Test.int_test(
            typeof(a * x^2);
            exclude = ["int_mul", "int_add", "int_add_mul"],
        )
    end
    return
end

function test_extension_scalar(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x)
    exclude = ["cube"]
    JuMP._MA.Test.scalar_test(x; exclude = exclude)
    for a in [2, 2im]
        JuMP._MA.Test.scalar_test(a * x + 3; exclude = exclude)
        JuMP._MA.Test.scalar_test(a * x^2 + 4x + 1; exclude = exclude)
    end
    return
end

function test_extension_quadratic(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:4])
    # Test is excluded because of https://github.com/jump-dev/MutableArithmetics.jl/issues/227
    JuMP._MA.Test.quadratic_test(x...; exclude = ["quadratic_add_canonical"])
    return
end

function test_extension_sparse(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, X11)
    @variable(model, X23)
    @variable(model, Xd[1:3, 1:3])
    JuMP._MA.Test.sparse_test(X11, X23, Xd)
    return
end

function test_extension_vector(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, x[1:3])
    JuMP._MA.Test.array_test(x)
    return
end

function test_extension_symmetric_matrix(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, y[1:2, 1:2], Symmetric)
    JuMP._MA.Test.array_test(y)
    return
end

function test_extension_nonsquare_matrix(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, z[1:2, 1:3])
    JuMP._MA.Test.array_test(z)
    return
end

function test_extension_DenseAxisVector(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, y[2:5])
    JuMP._MA.Test.array_test(y; exclude = ["matrix_vector", "non_array"])
    return
end

function test_extension_different_variables(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    x = @variable(model)
    y = DummyVariableRef()
    for a in (x, x + 1, x + im)
        @test JuMP._MA.promote_operation(+, typeof(a), typeof(y)) == DummyExpr
        @test JuMP._MA.promote_operation(-, typeof(a), typeof(y)) == DummyExpr
        @test JuMP._MA.promote_operation(+, typeof(y), typeof(a)) == DummyExpr
        @test JuMP._MA.promote_operation(-, typeof(y), typeof(a)) == DummyExpr
    end
    return
end

function test_scaling_errors()
    model = Model()
    @variable(model, x)
    @test_throws InexactError JuMP._MA.scaling(x + 1.0)
    @test_throws InexactError JuMP._MA.scaling(x^2 + 1.0)
    return
end

end
