#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestExpr

using JuMP
using Test

import LinearAlgebra
import SparseArrays

const MA = JuMP._MA

include(joinpath(@__DIR__, "utilities.jl"))

# For "expression^3 and unary*"
struct PowVariable <: AbstractVariableRef
    pow::Int
end

Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow * i)

Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)

Base.copy(x::PowVariable) = x

function test_extension_isequal_GenericAffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test isequal(x + 1, x + 1)
    return
end

function test_extension_hash_GenericAffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test hash(x + 1) == hash(x + 1)
    return
end

function test_extension_drop_zeros!_GenericAffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x[1:2])
    expr = x[1] + x[2] - x[2] + 1
    @test !isequal(expr, x[1] + 1)
    drop_zeros!(expr)
    @test isequal(expr, x[1] + 1)
    return
end

function test_extension_iszero_GenericAffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test !iszero(x + 1)
    @test !iszero(x + 0)
    @test iszero(0 * x + 0)
    @test iszero(x - x)
    return
end

function test_extension_isequal_GenericQuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test isequal(x^2 + 1, x^2 + 1)
    return
end

function test_extension_hash_GenericQuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test hash(x^2 + 1) == hash(x^2 + 1)
    return
end

function test_extension_drop_zeros!_GenericQuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x[1:2])
    expr = x[1]^2 + x[2]^2 - x[2]^2 + x[1] + x[2] - x[2] + 1
    @test !isequal(expr, x[1]^2 + x[1] + 1)
    drop_zeros!(expr)
    @test isequal(expr, x[1]^2 + x[1] + 1)
    return
end

function test_extension_iszero_GenericQuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x)
    @test !iszero(x^2 + 1)
    @test !iszero(x^2 + 0)
    @test !iszero(x^2 + 0 * x + 0)
    @test iszero(0 * x^2 + 0 * x + 0)
    @test iszero(x^2 - x^2)
    return
end

function test_value_GenericAffExpr()
    expr1 = GenericAffExpr(3.0, 3 => -5.0, 2 => 4.0)
    @test @inferred(value(-, expr1)) == 10.0
    expr2 = GenericAffExpr{Int,Int}(2)
    @test typeof(@inferred(value(i -> 1.0, expr2))) == Float64
    @test @inferred(value(i -> 1.0, expr2)) == 2.0
    return
end

function test_value_GenericQuadExpr()
    # 1 + 2x(1) + 3x(2)
    affine_term = GenericAffExpr(1.0, 1 => 2.0, 2 => 3.0)
    # 1 + 2x(1) + 3x(2) + 4x(1)^2 + 5x(1)*x(2) + 6x(2)^2
    expr = GenericQuadExpr(
        affine_term,
        UnorderedPair(1, 1) => 4.0,
        UnorderedPair(1, 2) => 5.0,
        UnorderedPair(2, 2) => 6.0,
    )
    @test typeof(@inferred(value(i -> 1.0, expr))) == Float64
    @test @inferred(value(i -> 1.0, expr)) == 21
    @test @inferred(value(i -> 2.0, expr)) == 71
    return
end

function test_add_to_expression_GenericAffExpr_V()
    aff = GenericAffExpr(1.0, :a => 2.0)
    @test isequal_canonical(
        add_to_expression!(aff, :b),
        GenericAffExpr(1.0, :a => 2.0, :b => 1.0),
    )
    return
end

function test_add_to_expression_GenericAffExpr_C()
    aff = GenericAffExpr(1.0, :a => 2.0)
    @test isequal_canonical(
        add_to_expression!(aff, 1.0),
        GenericAffExpr(2.0, :a => 2.0),
    )
    return
end

function test_extension_linear_terms_AffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    @variable(m, x[1:10])

    aff = 1 * x[1] + 2 * x[2]
    k = 0
    @test length(linear_terms(aff)) == 2
    for (coeff, var) in linear_terms(aff)
        if k == 0
            @test coeff == 1
            @test var === x[1]
        elseif k == 1
            @test coeff == 2
            @test var === x[2]
        end
        k += 1
    end
    @test k == 2
    return
end

function test_extension_linear_terms_empty_AffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    k = 0
    aff = zero(GenericAffExpr{T,VariableRefType})
    @test length(linear_terms(aff)) == 0
    for (coeff, var) in linear_terms(aff)
        k += 1
    end
    @test k == 0
    return
end

function test_extension_coefficient_AffExpr_VariableRefType(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    x = @variable(m, x)
    y = @variable(m, y)
    aff = @expression(m, 1.0 * x)
    @test coefficient(aff, x) == 1.0
    @test coefficient(aff, y) == 0.0
    return
end

function test_extension_coefficient_AffExpr_VariableRefType_VariableRefType(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    m = ModelType()
    x = @variable(m, x)
    aff = @expression(m, 1.0 * x)
    @test coefficient(aff, x, x) == 0.0
    return
end

function test_extension_coefficient_QuadExpr_VariableRefType(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    x = @variable(m, x)
    y = @variable(m, y)
    z = @variable(m, z)
    quad = @expression(m, T(6) * x^2 + T(5) * x * y + T(2) * y + T(3) * x)
    @test coefficient(quad, x) == T(3)
    @test coefficient(quad, y) == T(2)
    @test coefficient(quad, z) == T(0)
    return
end

function test_extension_coefficient_QuadExpr_VariableRefType_VariableRefType(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    m = ModelType()
    x = @variable(m, x)
    y = @variable(m, y)
    z = @variable(m, z)
    quad = @expression(m, T(6) * x^2 + T(5) * x * y + T(2) * y + T(3) * x)
    @test coefficient(quad, x, y) == T(5)
    @test coefficient(quad, x, x) == T(6)
    @test coefficient(quad, x, y) == coefficient(quad, y, x)
    @test coefficient(quad, z, z) == T(0)
    return
end

function test_extension_MA_add_mul(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    T = value_type(ModelType)
    model = ModelType()
    @variable(model, x)
    @variable(model, y)
    # MA.add_mul!!(ex::Number, c::Number, x::GenericAffExpr)
    aff = MA.add_mul!!(1, 2, GenericAffExpr(T(1), x => T(1)))
    @test isequal_canonical(aff, GenericAffExpr(T(3), x => T(2)))
    # MA.add_mul!!(ex::Number, c::Number, x::GenericQuadExpr) with c == 0
    QuadExprType = GenericQuadExpr{T,VariableRefType}
    quad = MA.add_mul!!(2, 0, QuadExprType())
    @test isequal_canonical(quad, convert(QuadExprType, 2))
    # MA.add_mul!!(ex::Number, c::VariableRef, x::VariableRef)"
    @test_expression_with_string MA.add_mul(5, x, y) "x*y + 5"
    @test_expression_with_string MA.add_mul!!(5, x, y) "x*y + 5"
    # MA.add_mul!!(ex::Number, c::T, x::T) where T<:GenericAffExpr" begin
    @test_expression_with_string MA.add_mul(1, 2x, x + 1) "2 x² + 2 x + 1"
    @test_expression_with_string MA.add_mul!!(1, 2x, x + 1) "2 x² + 2 x + 1"
    # MA.add_mul!!(ex::Number, c::GenericAffExpr{C,V}, x::V) where {C,V}" begin
    @test_expression_with_string MA.add_mul(1, 2x, x) "2 x² + 1"
    @test_expression_with_string MA.add_mul!!(1, 2x, x) "2 x² + 1"
    # MA.add_mul!!(ex::Number, c::GenericQuadExpr, x::Number)" begin
    @test_expression_with_string MA.add_mul(0, x^2, 1) "x²"
    @test_expression_with_string MA.add_mul!!(0, x^2, 1) "x²"
    # MA.add_mul!!(ex::Number, c::GenericQuadExpr, x::Number) with c == 0" begin
    @test_expression_with_string MA.add_mul(0, x^2, 0) "0"
    @test_expression_with_string MA.add_mul!!(0, x^2, 0) "0"
    # MA.add_mul!!(aff::AffExpr,c::VariableRef,x::AffExpr)" begin
    @test_expression_with_string MA.add_mul(2x, x, x + 1) "x² + 3 x"
    @test_expression_with_string MA.add_mul!!(2x, x, x + 1) "x² + 3 x"
    # MA.add_mul!!(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
    @test_expression_with_string MA.add_mul(2x, x, 1) "3 x"
    @test_expression_with_string MA.add_mul!!(2x, x, 1) "3 x"
    # MA.add_mul!!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V}" begin
    @test_expression_with_string MA.add_mul(2x, x^2, 1) "x² + 2 x"
    @test_expression_with_string MA.add_mul!!(2x, x^2, 1) "x² + 2 x"
    # MA.add_mul!!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V} with x == 0" begin
    @test_expression_with_string MA.add_mul(2x, x^2, 0) "2 x"
    @test_expression_with_string MA.add_mul!!(2x, x^2, 0) "2 x"
    # MA.add_mul!!(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) where {C,V} with c == 0" begin
    @test_expression_with_string MA.add_mul(2x, 0, x^2) "2 x"
    @test_expression_with_string MA.add_mul!!(2x, 0, x^2) "2 x"
    # MA.add_mul!!(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
    @test_expression_with_string MA.add_mul(2x, x + 1, x + 0) "x² + 3 x"
    @test_expression_with_string MA.add_mul!!(2x, x + 1, x + 0) "x² + 3 x"
    # MA.add_mul!!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
    @test_expression_with_string MA.add_mul(x^2, x + 1, 1) "x² + x + 1"
    @test_expression_with_string MA.add_mul!!(x^2, x + 1, 1) "x² + x + 1"
    # MA.add_mul!!(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V}) where {C,V}" begin
    @test_expression_with_string MA.add_mul(x^2, x, x + 1) "2 x² + x"
    @test_expression_with_string MA.add_mul!!(x^2, x, x + 1) "2 x² + x"
    # MA.add_mul!!(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number) where {C,V}" begin
    @test_expression_with_string MA.add_mul(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
    @test_expression_with_string MA.add_mul!!(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
    # MA.add_mul!!(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
    @test_expression_with_string MA.add_mul(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
    @test_expression_with_string MA.add_mul!!(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
    return
end

function test_extension_unary_plus_AffExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test_expression_with_string (+)(x + 1) "x + 1"
    return
end

function test_extension_unary_plus_QuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test_expression_with_string (+)(x^2 + 1) "x² + 1"
    return
end

function test_extension_sum_VectorVariableRef(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x[1:2])
    @test_expression_with_string sum(x) "x[1] + x[2]"
    return
end

function test_extension_expression_cubed(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    x = PowVariable(1)
    # Calls (*)((x*x)^6)
    y = @expression model (x * x)^3
    @test y.pow == 6
    z = @inferred (x * x)^3
    @test z.pow == 6
    return
end

function test_extension_ndims_QuadExpr(
    ModelType = Model,
    VariableRefType = VariableRef,
)
    model = ModelType()
    @variable(model, x)
    @test ndims(x^2 + 1) == 0
    return
end

function test_equal_0()
    model = Model()
    @variable(model, x)
    @test x + 0.0 != 0.0
    @test AffExpr(0.0) == 0.0
    @test AffExpr(1.0) == 1.0
    @test QuadExpr(AffExpr(0.0)) == 0.0
    @test QuadExpr(AffExpr(1.0)) == 1.0
    @test x^2 + 0.0 != 0.0
    return
end

function test_issue_2309()
    model = Model()
    @variable(model, x[1:10])
    I = SparseArrays.sparse(LinearAlgebra.Diagonal(ones(10)))
    A = I + LinearAlgebra.Diagonal(x)
    @test A isa SparseArrays.SparseMatrixCSC
    @test SparseArrays.nnz(A) == 10
    return
end

function test_eltype_QuadTermIterator()
    model = Model()
    @variable(model, x[1:2])
    y = x[1]^2 + x[2]^2
    iterator = quad_terms(y)
    @test eltype(iterator) == Tuple{Float64,VariableRef,VariableRef}
    return
end

function test_GenericQuadExpr_constructor()
    model = Model()
    @variable(model, x[1:2])
    y = 1 * x[1] + 2 * x[2] + 3 * x[1]^2 + 4 * x[2]^2
    map_coefficients_inplace!(c -> 2c, y)
    @test y == 2 * x[1] + 4 * x[2] + 6 * x[1]^2 + 8 * x[2]^2
    return
end

function test_GenericQuadExpr_map_coefficients_inplace!()
    model = Model()
    @variable(model, x[1:2])
    y = 1 * x[1] + 2 * x[2] + 3 * x[1]^2 + 4 * x[2]^2
    map_coefficients_inplace!(c -> 2c, y)
    @test y == 2 * x[1] + 4 * x[2] + 6 * x[1]^2 + 8 * x[2]^2
    return
end

function test_expression_ambiguities()
    # These tests use expressions with unusual key types so that we can test
    # the fallback methods needed to avoid method ambiguities.
    model = Model()
    quad = GenericQuadExpr{Int,Int}()
    aff = GenericAffExpr{Int,Int}(0, 1 => 1)
    @test add_to_expression!(quad, aff, 1) isa GenericQuadExpr
    @test quad == GenericQuadExpr{Int,Int}(aff)
    quad = GenericQuadExpr{Int,Int}()
    @test add_to_expression!(quad, 0, 0) isa GenericQuadExpr
    @variable(model, x)
    @test add_to_expression!(x + 1.0, 1.0, 2) isa GenericAffExpr
    @test add_to_expression!(x^2, 1.0, 2) isa GenericQuadExpr
    return
end

end  # TestExpr
