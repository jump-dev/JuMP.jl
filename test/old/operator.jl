#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/operator.jl
# Testing operator overloading is correct
#############################################################################
using JuMP
using Base.Test
using OffsetArrays

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const  leq = JuMP.repl[:leq]
const  geq = JuMP.repl[:geq]
const   eq = JuMP.repl[:eq]
const Vert = JuMP.repl[:Vert]
const sub2 = JuMP.repl[:sub2]

# For "DimensionMismatch when performing vector-matrix multiplication with custom types #988"
import Base: +, *
struct MyType{T}
    a::T
end
struct MySumType{T}
    a::T
end
Base.one(::Type{MyType{T}}) where {T} = MyType(one(T))
Base.zero(::Type{MySumType{T}}) where {T} = MySumType(zero(T))
Base.zero(::MySumType{T}) where {T} = MySumType(zero(T))
Base.transpose(t::MyType) = MyType(t.a)
Base.transpose(t::MySumType) = MySumType(t.a)
+(t1::MyT, t2::MyS) where {MyT<:Union{MyType, MySumType}, MyS<:Union{MyType, MySumType}} = MySumType(t1.a+t2.a)
*(t1::MyType{S}, t2::T) where {S, T} = MyType(t1.a*t2)
*(t1::S, t2::MyType{T}) where {S, T} = MyType(t1*t2.a)
*(t1::MyType{S}, t2::MyType{T}) where {S, T} = MyType(t1.a*t2.a)


@testset "Operator overloads" begin

    _lt(x,y) = (x.col < y.col)
    function sort_expr!(x::AffExpr)
        idx = sortperm(x.vars, lt=_lt)
        x.vars = x.vars[idx]
        x.coeffs = x.coeffs[idx]
        return x
    end

    vec_eq(x,y) = vec_eq([x;], [y;])

    function vec_eq(x::AbstractArray, y::AbstractArray)
        size(x) == size(y) || return false
        for i in 1:length(x)
            v, w = convert(AffExpr,x[i]), convert(AffExpr,y[i])
            sort_expr!(v)
            sort_expr!(w)
            string(v) == string(w) || return false
        end
        return true
    end

    function vec_eq(x::Array{QuadExpr}, y::Array{QuadExpr})
        size(x) == size(y) || return false
        for i in 1:length(x)
            string(x[i]) == string(y[i]) || return false
        end
        return true
    end

    @testset "Operators for non-Array AbstractArrays" begin
        m = Model()
        @variable(m, x[1:3])

        # This is needed to compare arrays that have nonstandard indexing
        elements_equal(A::AbstractArray{T, N}, B::AbstractArray{T, N}) where {T, N} = all(a == b for (a, b) in zip(A, B))

        for x2 in (OffsetArray(x, -length(x)), view(x, :), sparse(x))
            @test elements_equal(+x, +x2)
            @test elements_equal(-x, -x2)
            @test elements_equal(x + first(x), x2 + first(x2))
            @test elements_equal(x - first(x), x2 - first(x2))
            @test elements_equal(first(x) - x, first(x2) - x2)
            @test elements_equal(first(x) + x, first(x2) + x2)
            @test elements_equal(2 * x, 2 * x2)
            @test elements_equal(first(x) + x2, first(x2) + x)
            @test sum(x) == sum(x2)
            if !JuMP.one_indexed(x2)
                @test_throws DimensionMismatch x + x2
            end
        end
    end

    @testset "Norm and diagm for non-Array AbstractArrays" begin
        m = Model()
        @variable(m, x[1:3])

        for x2 in (OffsetArray(x, -length(x)), view(x, :), sparse(x))
            if !JuMP.one_indexed(x2)
                @test_throws AssertionError diagm(x2)
            else
                @test diagm(x) == diagm(x2)
            end
            @test norm(x).terms == norm(x2).terms
        end
    end

    @testset "DimensionMismatch when performing vector-matrix multiplication with custom types #988" begin
        m = Model()
        @variable m Q[1:3, 1:3] SDP

        x = [MyType(1), MyType(2), MyType(3)]
        y = Q * x
        z = x' * Q
        ElemT = MySumType{JuMP.GenericAffExpr{Float64,JuMP.Variable}}
        @test typeof(y) == Vector{ElemT}
        @test size(y) == (3,)
        @test typeof(z) == (isdefined(Base, :RowVector) ? RowVector{ElemT, ConjArray{ElemT, 1, Vector{ElemT}}} : Matrix{ElemT})
        @test size(z) == (1, 3)
        for i in 1:3
            # Q is symmetric
            a = zero(JuMP.GenericAffExpr{Float64,JuMP.Variable})
            a += Q[1,i]
            a += 2Q[2,i]
            a += 3Q[3,i]
            # Q[1,i] + 2Q[2,i] + 3Q[3,i] is rearranged as 2 Q[2,3] + Q[1,3] + 3 Q[3,3]
            @test z[i].a == y[i].a == a
        end
    end
end
