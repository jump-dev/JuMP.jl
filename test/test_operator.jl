#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestOperators

using JuMP
using Test

import LinearAlgebra
import SparseArrays

include(joinpath(@__DIR__, "utilities.jl"))

struct MyType{T}
    a::T
end

struct MySumType{T}
    a::T
end

Base.copy(t::MyType) = t

Base.zero(::Type{MyType{T}}) where {T} = MyType(zero(T))

Base.one(::Type{MyType{T}}) where {T} = MyType(one(T))

Base.zero(::Type{MySumType{T}}) where {T} = MySumType(zero(T))

Base.zero(::MySumType{T}) where {T} = MySumType(zero(T))

Base.transpose(t::MyType) = MyType(t.a)

Base.transpose(t::MySumType) = MySumType(t.a)

LinearAlgebra.adjoint(t::Union{MyType,MySumType}) = t

function Base.:(+)(
    t1::MyT,
    t2::MyS,
) where {MyT<:Union{MyType,MySumType},MyS<:Union{MyType,MySumType}}
    return MySumType(t1.a + t2.a)
end

Base.:(*)(t1::MyType{S}, t2::T) where {S,T} = MyType(t1.a * t2)

Base.:(*)(t1::S, t2::MyType{T}) where {S,T} = MyType(t1 * t2.a)

Base.:(*)(t1::MyType{S}, t2::MyType{T}) where {S,T} = MyType(t1.a * t2.a)

function JuMP.isequal_canonical(t::MySumType, s::MySumType)
    return isequal_canonical(t.a, s.a)
end

function JuMP.isequal_canonical(
    x::AbstractArray{<:MySumType},
    y::AbstractArray{<:MySumType},
)
    return size(x) == size(y) && all(isequal_canonical.(x, y))
end

function test_extension_promotion(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    I = Int
    V = VariableRefType
    A = GenericAffExpr{T,VariableRefType}
    Q = GenericQuadExpr{T,VariableRefType}
    @test promote_type(V, I) == A
    @test promote_type(I, V) == A
    @test promote_type(A, I) == A
    @test promote_type(I, A) == A
    @test promote_type(A, V) == A
    @test promote_type(V, A) == A
    @test promote_type(Q, I) == Q
    @test promote_type(I, Q) == Q
    @test promote_type(Q, A) == Q
    @test promote_type(A, Q) == Q
    @test promote_type(Q, V) == Q
    @test promote_type(V, Q) == Q
    return
end

function test_extension_broadcast_division_error(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2, 1:2])
    A = [1 2; 3 4]
    B = SparseArrays.sparse(A)
    y = SparseArrays.SparseMatrixCSC(
        2,
        2,
        copy(B.colptr),
        copy(B.rowval),
        vec(x),
    )
    @test_throws ErrorException A ./ x
    @test_throws ErrorException B ./ x
    @test_throws ErrorException A ./ y
    @test_throws ErrorException B ./ y
    # TODO: Refactor to avoid calling the internal JuMP function
    # `_densify_with_jump_eltype`.
    #z = _densify_with_jump_eltype((2 .* y) ./ 3)
    #@test isequal_canonical((2 .* x) ./ 3, z)
    #z = _densify_with_jump_eltype(2 * (y ./ 3))
    #@test isequal_canonical(2 .* (x ./ 3), z)
    #z = _densify_with_jump_eltype((x[1,1],) .* B)
    #@test isequal_canonical((x[1,1],) .* A, z)
    return
end

function test_extension_vectorized_comparisons(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x[1:3])
    A = [
        1 2 3
        0 4 5
        6 0 7
    ]
    B = SparseArrays.sparse(A)
    # force vector output
    cref1 = @constraint(m, reshape(x, (1, 3)) * A * x .>= 1)
    c1 = constraint_object.(cref1)
    f1 = map(c -> c.func, c1)
    @test isequal_canonical(
        f1,
        [
            x[1] * x[1] +
            2x[1] * x[2] +
            4x[2] * x[2] +
            9x[1] * x[3] +
            5x[2] * x[3] +
            7x[3] * x[3],
        ],
    )
    @test all(c -> c.set.lower == 1, c1)
    cref2 = @constraint(m, x' * A * x >= 1)
    c2 = constraint_object.(cref2)
    @test isequal_canonical(f1[1], c2.func)
    mat = [
        3x[1] + 12x[3] + 4x[2]
        2x[1] + 12x[2] + 10x[3]
        15x[1] + 5x[2] + 21x[3]
    ]
    cref3 = @constraint(m, (x'A)' + 2A * x .<= 1)
    c3 = constraint_object.(cref3)
    f3 = map(c -> c.func, c3)
    @test isequal_canonical(f3, mat)
    @test all(c -> c.set.upper == 1, c3)
    @test isequal_canonical((x'A)' + 2A * x, (x'A)' + 2B * x)
    @test isequal_canonical((x'A)' + 2A * x, (x'B)' + 2A * x)
    @test isequal_canonical((x'A)' + 2A * x, (x'B)' + 2B * x)
    @test isequal_canonical((x'A)' + 2A * x, JuMP._MA.@rewrite((x'A)' + 2A * x))
    @test isequal_canonical((x'A)' + 2A * x, JuMP._MA.@rewrite((x'B)' + 2A * x))
    @test isequal_canonical((x'A)' + 2A * x, JuMP._MA.@rewrite((x'A)' + 2B * x))
    @test isequal_canonical((x'A)' + 2A * x, JuMP._MA.@rewrite((x'B)' + 2B * x))
    cref4 = @constraint(m, -1 .<= (x'A)' + 2A * x .<= 1)
    c4 = constraint_object.(cref4)
    f4 = map(c -> c.func, c4)
    @test isequal_canonical(f4, mat)
    @test all(c -> c.set.lower == -1, c4)
    @test all(c -> c.set.upper == 1, c4)
    cref5 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= 1)
    c5 = constraint_object.(cref5)
    f5 = map(c -> c.func, c5)
    @test isequal_canonical(f5, mat)
    @test map(c -> c.set.lower, c5) == -[1:3;]
    @test all(c -> c.set.upper == 1, c4)
    cref6 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= [3:-1:1;])
    c6 = constraint_object.(cref6)
    f6 = map(c -> c.func, c6)
    @test isequal_canonical(f6, mat)
    @test map(c -> c.set.lower, c6) == -[1:3;]
    @test map(c -> c.set.upper, c6) == [3:-1:1;]
    cref7 = @constraint(m, -[1:3;] .<= (x'A)' + 2A * x .<= 3)
    c7 = constraint_object.(cref7)
    f7 = map(c -> c.func, c7)
    @test isequal_canonical(f7, mat)
    @test map(c -> c.set.lower, c7) == -[1:3;]
    @test all(c -> c.set.upper == 3, c7)
    return
end

function test_extension_custom_dimension_mismatch(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    ElemT = MySumType{GenericAffExpr{T,VariableRefType}}
    model = ModelType()
    @variable(model, Q[1:3, 1:3], PSD)
    x = [MyType(1), MyType(2), MyType(3)]
    y = Q * x
    z = x' * Q
    @test y isa Vector{ElemT}
    @test size(y) == (3,)
    @test z isa LinearAlgebra.Adjoint{ElemT,<:Any}
    @test size(z) == (1, 3)
    for i in 1:3
        # Q is symmetric
        a = zero(GenericAffExpr{T,VariableRefType})
        a += Q[1, i]
        a += 2Q[2, i]
        a += 3Q[3, i]
        # Q[1,i] + 2Q[2,i] + 3Q[3,i] is rearranged as 2 Q[2,3] + Q[1,3] + 3 Q[3,3]
        @test z[i].a == y[i].a == a
    end
    return
end

function test_extension_matrix_multiplication(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    ElemT = MySumType{GenericAffExpr{T,VariableRefType}}
    model = ModelType()
    @variable(model, Q[1:3, 1:3], PSD)
    X = MyType.((1:3)' .+ (1:3))
    @test_expression Q * X
    Y = Q * X
    @test_expression X' * Q
    Z = X' * Q
    @test Y isa Matrix{ElemT}
    @test size(Y) == (3, 3)
    @test Z isa Matrix{ElemT}
    @test size(Z) == (3, 3)
    @test isequal_canonical(Z', Y)
    @test_expression Q * X'
    Y = Q * X'
    @test_expression X * Q
    Z = X * Q
    @test Y isa Matrix{ElemT}
    @test size(Y) == (3, 3)
    @test Z isa Matrix{ElemT}
    @test size(Z) == (3, 3)
    @test isequal_canonical(Z', Y)
    return
end

function test_extension_operator_warn(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable model x[1:51]
    # JuMPExtension does not have the `operator_counter` field
    if ModelType <: Model
        @test model.operator_counter == 0
    end
    # Triggers the increment of operator_counter since sum(x) has more than 50 terms
    @test_expression(sum(x) + 2x[1])
    if ModelType <: Model
        # The following check verifies that this test covers the code incrementing `operator_counter`
        @test model.operator_counter == 1
    end
    return
end

function test_extension_uniform_scaling(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test_expression_with_string x + 2 * LinearAlgebra.I "x + 2"
    @test_expression_with_string (x + 1) + LinearAlgebra.I "x + 2"
    @test_expression_with_string x - 2 * LinearAlgebra.I "x - 2"
    @test_expression_with_string (x - 1) - LinearAlgebra.I "x - 2"
    @test_expression_with_string 2 * LinearAlgebra.I + x "x + 2"
    @test_expression_with_string LinearAlgebra.I + (x + 1) "x + 2"
    @test_expression_with_string 2 * LinearAlgebra.I - x "-x + 2"
    @test_expression_with_string (2im * LinearAlgebra.I) - x "-x + 2im"
    @test_expression_with_string LinearAlgebra.I - (x - 1) "-x + 2"
    @test_expression_with_string (LinearAlgebra.I * im) - (x - 1) "-x + (1 + im)"
    @test_expression_with_string LinearAlgebra.I * x "x"
    @test_expression_with_string (LinearAlgebra.I * im) * x "x im"
    @test_expression_with_string LinearAlgebra.I * (x + 1) "x + 1"
    @test_expression_with_string (LinearAlgebra.I * im) * (x + 1) "x im + im"
    @test_expression_with_string (x + 1) * LinearAlgebra.I "x + 1"
    @test_expression_with_string (x + 1) * (LinearAlgebra.I * im) "x im + im"
    return
end

function test_extension_basic_operators(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    aff = @inferred 7.1 * x + 2.5
    @test_expression_with_string 7.1 * x + 2.5 "7.1 x + 2.5"
    aff2 = @inferred 1.2 * y + 1.2
    @test_expression_with_string 1.2 * y + 1.2 "1.2 y + 1.2"
    q = @inferred 2.5 * y * z + aff
    @test_expression_with_string 2.5 * y * z + aff "2.5 y*z + 7.1 x + 2.5"
    q2 = @inferred 8 * x * z + aff2
    @test_expression_with_string 8 * x * z + aff2 "8 x*z + 1.2 y + 1.2"
    @test_expression_with_string 2 * x * x + 1 * y * y + z + 3 "2 x² + y² + z + 3"
    return
end

function test_extension_basic_operators_number(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    aff = @inferred 7.1 * x + 2.5
    q = @inferred 2.5 * y * z + aff
    # 1-1 Number--Number
    # ???
    # 1-2 Number--Variable
    @test_expression_with_string 4.13 + w "w + 4.13"
    @test_expression_with_string 3.16 - w "-w + 3.16"
    @test_expression_with_string 5.23 * w "5.23 w"
    @test_throws ErrorException 2.94 / w
    # 1-3 Number--AffExpr
    @test_expression_with_string 1.5 + aff "7.1 x + 4"
    @test_expression_with_string 1.5 - aff "-7.1 x - 1"
    @test_expression_with_string 2 * aff "14.2 x + 5"
    @test_throws ErrorException 2 / aff
    # 1-4 Number--QuadExpr
    @test_expression_with_string 1.5 + q "2.5 y*z + 7.1 x + 4"
    @test_expression_with_string 1.5 - q "-2.5 y*z - 7.1 x - 1"
    @test_expression_with_string 2 * q "5 y*z + 14.2 x + 5"
    @test_throws ErrorException 2 / q
    return
end

function test_extension_basic_operators_variable(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    β = T(71) / T(10)
    aff = @inferred β * x + T(25) / T(10)
    q = @inferred T(25) / T(10) * y * z + aff
    # 2-0 Variable unary
    @test (+x) === x
    @test_expression_with_string -x "-x"
    # 2-1 Variable--Number
    α = T(413) / T(100)
    @test_expression_with_string w + α "w + 4.13"
    @test_expression_with_string w - α "w - 4.13"
    @test_expression_with_string w * α "4.13 w"
    @test_expression_with_string w / T(2) "0.5 w"
    @test w == w
    @test_expression_with_string x * y - 1 "x*y - 1"
    @test_expression_with_string x^2 "x²"
    @test_expression_with_string x^1 "x"
    @test_expression_with_string x^0 "1"
    @test_throws ErrorException x^3
    @test_throws ErrorException x^(T(15) / T(10))
    # 2-2 Variable--Variable
    @test_expression_with_string w + x "w + x"
    @test_expression_with_string w - x "w - x"
    @test_expression_with_string w * x "w*x"
    @test_expression_with_string x - x "0"
    @test_throws ErrorException w / x
    @test_expression_with_string y * z - x "y*z - x"
    # 2-3 Variable--AffExpr
    @test_expression_with_string z + aff "z + 7.1 x + 2.5"
    @test_expression_with_string z - aff "z - 7.1 x - 2.5"
    @test_expression_with_string z * aff "7.1 z*x + 2.5 z"
    @test_throws ErrorException z / aff
    @test_throws MethodError z ≤ aff
    @test_expression_with_string β * x - aff "0 x - 2.5"
    # 2-4 Variable--QuadExpr
    @test_expression_with_string w + q "2.5 y*z + w + 7.1 x + 2.5"
    @test_expression_with_string w - q "-2.5 y*z + w - 7.1 x - 2.5"
    @test_throws ErrorException w * q
    @test_throws ErrorException w / q
    @test transpose(x) === x
    @test conj(x) === x
    return
end

function test_extension_basic_operators_affexpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    aff = @inferred T(71) / T(10) * x + T(25) / T(10)
    aff2 = @inferred T(12) / T(10) * y + T(12) / T(10)
    q = @inferred T(25) / T(10) * y * z + aff
    # 3-0 AffExpr unary
    @test_expression_with_string +aff "7.1 x + 2.5"
    @test_expression_with_string -aff "-7.1 x - 2.5"
    # 3-1 AffExpr--Number
    @test_expression_with_string aff + 1.5 "7.1 x + 4"
    @test_expression_with_string aff - 1.5 "7.1 x + 1"
    @test_expression_with_string aff * 2 "14.2 x + 5"
    @test_expression_with_string aff / 2 "3.55 x + 1.25"
    @test_throws MethodError aff ≤ 1
    @test aff == aff
    @test_throws MethodError aff ≥ 1
    @test_expression_with_string aff - 1 "7.1 x + 1.5"
    @test_expression_with_string aff^2 "50.41 x² + 35.5 x + 6.25"
    @test_expression_with_string (7.1 * x + 2.5)^2 "50.41 x² + 35.5 x + 6.25"
    @test_expression_with_string aff^1 "7.1 x + 2.5"
    @test_expression_with_string (7.1 * x + 2.5)^1 "7.1 x + 2.5"
    @test_expression_with_string aff^0 "1"
    @test_expression_with_string (7.1 * x + 2.5)^0 "1"
    @test_throws ErrorException aff^3
    @test_throws ErrorException (7.1 * x + 2.5)^3
    @test_throws ErrorException aff^1.5
    @test_throws ErrorException (7.1 * x + 2.5)^1.5
    # 3-2 AffExpr--Variable
    @test_expression_with_string aff + z "7.1 x + z + 2.5"
    @test_expression_with_string aff - z "7.1 x - z + 2.5"
    @test_expression_with_string aff * z "7.1 x*z + 2.5 z"
    @test_throws ErrorException aff / z
    @test_expression_with_string aff - 7.1 * x "0 x + 2.5"
    # 3-3 AffExpr--AffExpr
    @test_expression_with_string aff + aff2 "7.1 x + 1.2 y + 3.7"
    @test_expression_with_string aff - aff2 "7.1 x - 1.2 y + 1.3"
    @test_expression_with_string aff * aff2 "8.52 x*y + 3 y + 8.52 x + 3"
    @test string((x + x) * (x + 3)) == string((x + 3) * (x + x))  # Issue #288
    @test_throws ErrorException aff / aff2
    @test_expression_with_string aff - aff "0 x"
    # 4-4 AffExpr--QuadExpr
    @test_expression_with_string aff2 + q "2.5 y*z + 1.2 y + 7.1 x + 3.7"
    @test_expression_with_string aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    @test_throws ErrorException aff2 * q
    @test_throws ErrorException aff2 / q
    @test transpose(aff) === aff
    @test conj(aff) === aff
    return
end

function test_extension_basic_operators_quadexpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, w)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z)
    aff = @inferred 7.1 * x + 2.5
    aff2 = @inferred 1.2 * y + 1.2
    q = @inferred 2.5 * y * z + aff
    q2 = @inferred 8 * x * z + aff2
    # 4-0 QuadExpr unary
    @test_expression_with_string +q "2.5 y*z + 7.1 x + 2.5"
    @test_expression_with_string -q "-2.5 y*z - 7.1 x - 2.5"
    # 4-1 QuadExpr--Number
    @test_expression_with_string q + 1.5 "2.5 y*z + 7.1 x + 4"
    @test_expression_with_string q - 1.5 "2.5 y*z + 7.1 x + 1"
    @test_expression_with_string q * 2 "5 y*z + 14.2 x + 5"
    @test_expression_with_string q / 2 "1.25 y*z + 3.55 x + 1.25"
    @test q == q
    @test_expression_with_string aff2 - q "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    # 4-2 QuadExpr--Variable
    @test_expression_with_string q + w "2.5 y*z + 7.1 x + w + 2.5"
    @test_expression_with_string q - w "2.5 y*z + 7.1 x - w + 2.5"
    @test_throws ErrorException q * w
    @test_throws ErrorException q / w
    # 4-3 QuadExpr--AffExpr
    @test_expression_with_string q + aff2 "2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @test_expression_with_string q - aff2 "2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @test_throws ErrorException q * aff2
    @test_throws ErrorException q / aff2
    # 4-4 QuadExpr--QuadExpr
    @test_expression_with_string q + q2 "2.5 y*z + 8 x*z + 7.1 x + 1.2 y + 3.7"
    @test_expression_with_string q - q2 "2.5 y*z - 8 x*z + 7.1 x - 1.2 y + 1.3"
    @test_throws ErrorException q * q2
    @test_throws ErrorException q / q2
    @test transpose(q) === q
    @test conj(q) === q
    return
end

function test_extension_dot(ModelType = Model, VariableRefType = VariableRef)
    model = ModelType()
    @variable(model, 0 ≤ x[1:3] ≤ 1)
    @test_expression_with_string LinearAlgebra.dot(x[1], x[1]) "x[1]²"
    @test_expression_with_string LinearAlgebra.dot(2, x[1]) "2 x[1]"
    @test_expression_with_string LinearAlgebra.dot(x[1], 2) "2 x[1]"
    c = vcat(1:3)
    @test_expression_with_string LinearAlgebra.dot(c, x) "x[1] + 2 x[2] + 3 x[3]"
    @test_expression_with_string LinearAlgebra.dot(x, c) "x[1] + 2 x[2] + 3 x[3]"
    A = [1 3; 2 4]
    @variable(model, 1 ≤ y[1:2, 1:2] ≤ 1)
    @test_expression_with_string LinearAlgebra.dot(A, y) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    @test_expression_with_string LinearAlgebra.dot(y, A) "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    B = ones(2, 2, 2)
    @variable(model, 0 ≤ z[1:2, 1:2, 1:2] ≤ 1)
    @test_expression_with_string LinearAlgebra.dot(B, z) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @test_expression_with_string LinearAlgebra.dot(z, B) "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @objective(
        model,
        Max,
        LinearAlgebra.dot(x, ones(3)) - LinearAlgebra.dot(y, ones(2, 2))
    )
    for i in 1:3
        set_start_value(x[i], 1)
    end
    for i in 1:2, j in 1:2
        set_start_value(y[i, j], 1)
    end
    for i in 1:2, j in 1:2, k in 1:2
        set_start_value(z[i, j, k], 1)
    end
    @test LinearAlgebra.dot(c, start_value.(x)) ≈ 6
    @test LinearAlgebra.dot(A, start_value.(y)) ≈ 10
    @test LinearAlgebra.dot(B, start_value.(z)) ≈ 8
    return
end

function test_model_dot_fallback()
    # Check that dot is not falling back to default, inefficient
    # addition (JuMP PR #943).
    model = Model()
    @test model.operator_counter == 0
    @variable(model, x[1:100])
    set_start_value.(x, 1:100)
    @expression(model, test_sum, sum(x[i] * i for i in 1:100))
    @expression(model, test_dot1, LinearAlgebra.dot(x, 1:100))
    @expression(model, test_dot2, LinearAlgebra.dot(1:100, x))
    @test model.operator_counter == 0
    test_add = test_dot1 + test_dot2
    @test model.operator_counter == 1  # Check triggerable.
    test_sum_value = value(start_value, test_sum)
    @test test_sum_value ≈ value(start_value, test_dot1)
    @test test_sum_value ≈ value(start_value, test_dot2)
    return
end

function test_extension_higher_level(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, 0 ≤ matrix[1:3, 1:3] ≤ 1, start = 1)
    # "sum(j::DenseAxisArray{Variable})"
    @test_expression_with_string sum(matrix) "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    # "sum(j::DenseAxisArray{T}) where T<:Real"
    @test sum(start_value.(matrix)) ≈ 9
    # "sum(j::Array{VariableRef})"
    @test string(sum(matrix[1:3, 1:3])) == string(sum(matrix))
    # "sum(affs::Array{AffExpr})"
    @test_expression_with_string sum([2 * matrix[i, j] for i in 1:3, j in 1:3]) "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"
    # "sum(quads::Array{QuadExpr})"
    @test_expression_with_string sum([
        2 * matrix[i, j]^2 for i in 1:3, j in 1:3
    ]) "2 matrix[1,1]² + 2 matrix[2,1]² + 2 matrix[3,1]² + 2 matrix[1,2]² + 2 matrix[2,2]² + 2 matrix[3,2]² + 2 matrix[1,3]² + 2 matrix[2,3]² + 2 matrix[3,3]²"
    S = [1, 3]
    @variable(model, x[S], start = 1)
    # "sum(j::JuMPDict{VariableRef})"
    @test_expression sum(x)
    @test length(string(sum(x))) == 11 # order depends on hashing
    @test occursin("x[1]", string(sum(x)))
    @test occursin("x[3]", string(sum(x)))
    # "sum(j::JuMPDict{T}) where T<:Real"
    @test sum(start_value.(x)) == 2
    return
end

function test_complex_mult_variable()
    model = Model()
    @variable(model, x[1:3])
    A = rand(ComplexF64, 3, 3)
    @test (@inferred A * x) isa Vector{GenericAffExpr{ComplexF64,VariableRef}}
    return
end

function test_complex_pow()
    model = Model()
    @variable(model, x)
    y = (1.0 + 2.0im) * x
    @test y^0 == (1.0 + 0im)
    @test y^1 == 0 * y * y + y
    @test y^2 == y * y
    @test_throws ErrorException y^3
    return
end

end
