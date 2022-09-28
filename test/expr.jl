#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

using JuMP
using LinearAlgebra
using SparseArrays
using Test

const MA = JuMP._MA

include(joinpath(@__DIR__, "utilities.jl"))

@static if !(:JuMPExtension in names(Main))
    include(joinpath(@__DIR__, "JuMPExtension.jl"))
end

# For "expression^3 and unary*"
struct PowVariable <: JuMP.AbstractVariableRef
    pow::Int
end
Base.:^(x::PowVariable, i::Int) = PowVariable(x.pow * i)
Base.:*(x::PowVariable, y::PowVariable) = PowVariable(x.pow + y.pow)
Base.copy(x::PowVariable) = x

function expressions_test(
    ModelType::Type{<:JuMP.AbstractModel},
    VariableRefType::Type{<:JuMP.AbstractVariableRef},
)
    AffExprType = JuMP.GenericAffExpr{Float64,VariableRefType}
    QuadExprType = JuMP.GenericQuadExpr{Float64,VariableRefType}

    @testset "isequal(::GenericAffExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test isequal(x + 1, x + 1)
    end

    @testset "hash(::GenericAffExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test hash(x + 1) == hash(x + 1)
    end

    @testset "drop_zeros!(::GenericAffExpr)" begin
        m = ModelType()
        @variable(m, x[1:2])
        expr = x[1] + x[2] - x[2] + 1
        @test !isequal(expr, x[1] + 1)
        JuMP.drop_zeros!(expr)
        @test isequal(expr, x[1] + 1)
    end

    @testset "iszero(::GenericAffExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test !iszero(x + 1)
        @test !iszero(x + 0)
        @test iszero(0 * x + 0)
        @test iszero(x - x)
    end

    @testset "isequal(::GenericQuadExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test isequal(x^2 + 1, x^2 + 1)
    end

    @testset "hash(::GenericQuadExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test hash(x^2 + 1) == hash(x^2 + 1)
    end

    @testset "drop_zeros!(::GenericQuadExpr)" begin
        m = ModelType()
        @variable(m, x[1:2])
        expr = x[1]^2 + x[2]^2 - x[2]^2 + x[1] + x[2] - x[2] + 1
        @test !isequal(expr, x[1]^2 + x[1] + 1)
        JuMP.drop_zeros!(expr)
        @test isequal(expr, x[1]^2 + x[1] + 1)
    end

    @testset "iszero(::GenericQuadExpr)" begin
        m = ModelType()
        @variable(m, x)
        @test !iszero(x^2 + 1)
        @test !iszero(x^2 + 0)
        @test !iszero(x^2 + 0 * x + 0)
        @test iszero(0 * x^2 + 0 * x + 0)
        @test iszero(x^2 - x^2)
    end

    @testset "value for GenericAffExpr" begin
        expr1 = JuMP.GenericAffExpr(3.0, 3 => -5.0, 2 => 4.0)
        @test @inferred(JuMP.value(-, expr1)) == 10.0
        expr2 = JuMP.GenericAffExpr{Int,Int}(2)
        @test typeof(@inferred(JuMP.value(i -> 1.0, expr2))) == Float64
        @test @inferred(JuMP.value(i -> 1.0, expr2)) == 2.0
    end

    @testset "value for GenericQuadExpr" begin
        # 1 + 2x(1) + 3x(2)
        affine_term = JuMP.GenericAffExpr(1.0, 1 => 2.0, 2 => 3.0)
        # 1 + 2x(1) + 3x(2) + 4x(1)^2 + 5x(1)*x(2) + 6x(2)^2
        expr = JuMP.GenericQuadExpr(
            affine_term,
            JuMP.UnorderedPair(1, 1) => 4.0,
            JuMP.UnorderedPair(1, 2) => 5.0,
            JuMP.UnorderedPair(2, 2) => 6.0,
        )
        @test typeof(@inferred(JuMP.value(i -> 1.0, expr))) == Float64
        @test @inferred(JuMP.value(i -> 1.0, expr)) == 21
        @test @inferred(JuMP.value(i -> 2.0, expr)) == 71
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::V)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(
            JuMP.add_to_expression!(aff, :b),
            JuMP.GenericAffExpr(1.0, :a => 2.0, :b => 1.0),
        )
    end

    @testset "add_to_expression!(::GenericAffExpr{C,V}, ::C)" begin
        aff = JuMP.GenericAffExpr(1.0, :a => 2.0)
        @test JuMP.isequal_canonical(
            JuMP.add_to_expression!(aff, 1.0),
            JuMP.GenericAffExpr(2.0, :a => 2.0),
        )
    end

    @testset "linear_terms(::AffExpr)" begin
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
    end

    @testset "linear_terms(::AffExpr) for empty expression" begin
        k = 0
        aff = zero(AffExprType)
        @test length(linear_terms(aff)) == 0
        for (coeff, var) in linear_terms(aff)
            k += 1
        end
        @test k == 0
    end

    @testset "coefficient(aff::AffExpr, v::VariableRef)" begin
        m = ModelType()
        x = @variable(m, x)
        y = @variable(m, y)
        aff = @expression(m, 1.0 * x)
        @test coefficient(aff, x) == 1.0
        @test coefficient(aff, y) == 0.0
    end

    @testset "coefficient(aff::AffExpr, v1::VariableRef, v2::VariableRef)" begin
        m = ModelType()
        x = @variable(m, x)
        aff = @expression(m, 1.0 * x)
        @test coefficient(aff, x, x) == 0.0
    end

    @testset "coefficient(quad::QuadExpr, v::VariableRef)" begin
        m = ModelType()
        x = @variable(m, x)
        y = @variable(m, y)
        z = @variable(m, z)
        quad = @expression(m, 6.0 * x^2 + 5.0 * x * y + 2.0 * y + 3.0 * x)
        @test coefficient(quad, x) == 3.0
        @test coefficient(quad, y) == 2.0
        @test coefficient(quad, z) == 0.0
    end

    @testset "coefficient(quad::Quad, v1::VariableRef, v2::VariableRef)" begin
        m = ModelType()
        x = @variable(m, x)
        y = @variable(m, y)
        z = @variable(m, z)
        quad = @expression(m, 6.0 * x^2 + 5.0 * x * y + 2.0 * y + 3.0 * x)
        @test coefficient(quad, x, y) == 5.0
        @test coefficient(quad, x, x) == 6.0
        @test coefficient(quad, x, y) == coefficient(quad, y, x)
        @test coefficient(quad, z, z) == 0.0
    end

    @testset "MA.add_mul!!(ex::Number, c::Number, x::GenericAffExpr)" begin
        aff = MA.add_mul!!(1.0, 2.0, JuMP.GenericAffExpr(1.0, :a => 1.0))
        @test JuMP.isequal_canonical(aff, JuMP.GenericAffExpr(3.0, :a => 2.0))
    end

    @testset "MA.add_mul!!(ex::Number, c::Number, x::GenericQuadExpr) with c == 0" begin
        quad = MA.add_mul!!(2.0, 0.0, QuadExprType())
        @test JuMP.isequal_canonical(quad, convert(QuadExprType, 2.0))
    end

    @testset "MA.add_mul!!(ex::Number, c::VariableRef, x::VariableRef)" begin
        model = ModelType()
        @variable(model, x)
        @variable(model, y)
        @test_expression_with_string MA.add_mul(5.0, x, y) "x*y + 5"
        @test_expression_with_string MA.add_mul!!(5.0, x, y) "x*y + 5"
    end

    @testset "MA.add_mul!!(ex::Number, c::T, x::T) where T<:GenericAffExpr" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(1.0, 2x, x + 1) "2 x² + 2 x + 1"
        @test_expression_with_string MA.add_mul!!(1.0, 2x, x + 1) "2 x² + 2 x + 1"
    end

    @testset "MA.add_mul!!(ex::Number, c::GenericAffExpr{C,V}, x::V) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(1.0, 2x, x) "2 x² + 1"
        @test_expression_with_string MA.add_mul!!(1.0, 2x, x) "2 x² + 1"
    end

    @testset "MA.add_mul!!(ex::Number, c::GenericQuadExpr, x::Number)" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(0.0, x^2, 1.0) "x²"
        @test_expression_with_string MA.add_mul!!(0.0, x^2, 1.0) "x²"
    end

    @testset "MA.add_mul!!(ex::Number, c::GenericQuadExpr, x::Number) with c == 0" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(0.0, x^2, 0.0) "0"
        @test_expression_with_string MA.add_mul!!(0.0, x^2, 0.0) "0"
    end

    @testset "MA.add_mul!!(aff::AffExpr,c::VariableRef,x::AffExpr)" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, x, x + 1) "x² + 3 x"
        @test_expression_with_string MA.add_mul!!(2x, x, x + 1) "x² + 3 x"
    end

    @testset "MA.add_mul!!(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, x, 1) "3 x"
        @test_expression_with_string MA.add_mul!!(2x, x, 1) "3 x"
    end

    @testset "MA.add_mul!!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, x^2, 1) "x² + 2 x"
        @test_expression_with_string MA.add_mul!!(2x, x^2, 1) "x² + 2 x"
    end

    @testset "MA.add_mul!!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) where {C,V} with x == 0" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, x^2, 0) "2 x"
        @test_expression_with_string MA.add_mul!!(2x, x^2, 0) "2 x"
    end

    @testset "MA.add_mul!!(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) where {C,V} with c == 0" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, 0, x^2) "2 x"
        @test_expression_with_string MA.add_mul!!(2x, 0, x^2) "2 x"
    end

    @testset "MA.add_mul!!(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(2x, x + 1, x + 0) "x² + 3 x"
        @test_expression_with_string MA.add_mul!!(2x, x + 1, x + 0) "x² + 3 x"
    end

    @testset "MA.add_mul!!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(x^2, x + 1, 1) "x² + x + 1"
        @test_expression_with_string MA.add_mul!!(x^2, x + 1, 1) "x² + x + 1"
    end

    @testset "MA.add_mul!!(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V}) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(x^2, x, x + 1) "2 x² + x"
        @test_expression_with_string MA.add_mul!!(x^2, x, x + 1) "2 x² + x"
    end

    @testset "MA.add_mul!!(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
        @test_expression_with_string MA.add_mul!!(x^2 + x, x^2 + x, 2.0) "3 x² + 3 x"
    end

    @testset "MA.add_mul!!(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string MA.add_mul(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
        @test_expression_with_string MA.add_mul!!(x^2 + x, x + 0, x + 1) "2 x² + 2 x"
    end

    @testset "(+)(::AffExpr)" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string (+)(x + 1) "x + 1"
    end

    @testset "(+)(::QuadExpr)" begin
        model = ModelType()
        @variable(model, x)
        @test_expression_with_string (+)(x^2 + 1) "x² + 1"
    end

    @testset "sum(::Vector{VariableRef})" begin
        model = ModelType()
        @variable(model, x[1:2])
        @test_expression_with_string sum(x) "x[1] + x[2]"
    end

    @testset "expression^3 and unary*" begin
        model = ModelType()
        x = PowVariable(1)
        # Calls (*)((x*x)^6)
        y = @expression model (x * x)^3
        @test y.pow == 6
        z = @inferred (x * x)^3
        @test z.pow == 6
    end

    @testset "ndims(::QuadExpr)" begin
        model = ModelType()
        @variable(model, x)
        @test ndims(x^2 + 1) == 0
    end
end

@testset "Expressions for JuMP.Model" begin
    expressions_test(Model, VariableRef)
end

@testset "Expressions for JuMPExtension.MyModel" begin
    expressions_test(JuMPExtension.MyModel, JuMPExtension.MyVariableRef)
end

@testset "==0" begin
    model = Model()
    @variable(model, x)
    @test x + 0.0 != 0.0
    @test AffExpr(0.0) == 0.0
    @test AffExpr(1.0) == 1.0
    @test QuadExpr(AffExpr(0.0)) == 0.0
    @test QuadExpr(AffExpr(1.0)) == 1.0
    @test x^2 + 0.0 != 0.0
end

if VERSION >= v"1.6"
    # Don't test this on Julia 1.0. There are some issues adding the sparse and
    # diagonal matrices.
    @testset "issue_2309" begin
        model = Model()
        @variable(model, x[1:10])
        I = SparseArrays.sparse(LinearAlgebra.Diagonal(ones(10)))
        A = I + LinearAlgebra.Diagonal(x)
        @test A isa SparseArrays.SparseMatrixCSC
        @test SparseArrays.nnz(A) == 10
    end
end

@testset "eltype-QuadTermIterator" begin
    model = Model()
    @variable(model, x[1:2])
    y = x[1]^2 + x[2]^2
    iterator = quad_terms(y)
    @test eltype(iterator) == Tuple{Float64,VariableRef,VariableRef}
end

@testset "GenericQuadExpr_constructor" begin
    model = Model()
    @variable(model, x[1:2])
    y = 1 * x[1] + 2 * x[2] + 3 * x[1]^2 + 4 * x[2]^2
    map_coefficients_inplace!(c -> 2c, y)
    @test y == 2 * x[1] + 4 * x[2] + 6 * x[1]^2 + 8 * x[2]^2
end

@testset "GenericQuadExpr_map_coefficients_inplace!" begin
    model = Model()
    @variable(model, x[1:2])
    y = 1 * x[1] + 2 * x[2] + 3 * x[1]^2 + 4 * x[2]^2
    map_coefficients_inplace!(c -> 2c, y)
    @test y == 2 * x[1] + 4 * x[2] + 6 * x[1]^2 + 8 * x[2]^2
end

@testset "expression_ambiguities" begin
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
end
