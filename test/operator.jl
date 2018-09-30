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
using JuMP, Compat.Test, Compat, Compat.LinearAlgebra, Compat.SparseArrays
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
Base.convert(::Type{MySumType{T}}, t::MyType{T}) where {T} = MySumType(t.a)
if VERSION ≥ v"0.7-"
    Base.adjoint(t::MyType) = t
    Base.adjoint(t::MySumType) = t
end


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

    function vec_eq(x::AbstractArray{QuadExpr}, y::AbstractArray{QuadExpr})
        size(x) == size(y) || return false
        for i in 1:length(x)
            string(x[i]) == string(y[i]) || return false
        end
        return true
    end

    @testset "Testing basic operator overloads" begin
        m = Model()
        @variable(m, w)
        @variable(m, x)
        @variable(m, y)
        @variable(m, z)
        aff = 7.1 * x + 2.5
        @test string(aff) == "7.1 x + 2.5"
        aff2 = 1.2 * y + 1.2
        @test string(aff2) == "1.2 y + 1.2"
        q = 2.5 * y * z + aff
        @test string(q) == "2.5 y*z + 7.1 x + 2.5"
        q2 = 8 * x * z + aff2
        @test string(q2) == "8 x*z + 1.2 y + 1.2"
        q3 = 2 * x * x + 1 * y * y + z + 3
        @test string(q3) == "2 x² + y² + z + 3"
        nrm = norm([w,1-w])
        @test string(nrm) == "$Vert[w,-w + 1]$Vert$sub2"
        socexpr = 1.5*nrm - 2 - w
        @test string(socexpr) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 2"

        @test isequal(3w + 2y, 3w +2y) == true
        @test isequal(3w + 2y + 1, 3w + 2y) == false

        # Different objects that must all interact:
        # 1. Number
        # 2. Variable
        # 3. AffExpr
        # 4. QuadExpr

        # 1. Number tests
        @testset "Number--???" begin
        # 1-1 Number--Number - nope!
        # 1-2 Number--Variable
        @test string(4.13 + w) == "w + 4.13"
        @test string(3.16 - w) == "-w + 3.16"
        @test string(5.23 * w) == "5.23 w"
        @test_throws ErrorException 2.94 / w
        @test string(@LinearConstraint(2.1 ≤ w)) == "-w $leq -2.1"
        @test string(@LinearConstraint(2.1 == w)) == "-w $eq -2.1"
        @test string(@LinearConstraint(2.1 ≥ w)) == "-w $geq -2.1"
        # 1-3 Number--Norm
        @test string(4.13 + nrm) == "$Vert[w,-w + 1]$Vert$sub2 + 4.13"
        @test string(3.16 - nrm) == "-1.0 $Vert[w,-w + 1]$Vert$sub2 + 3.16"
        @test string(5.23 * nrm) == "5.23 $Vert[w,-w + 1]$Vert$sub2"
        @test_throws MethodError 2.94 / nrm
        @test_throws ErrorException @SOCConstraint(2.1 ≤ nrm)
        @test_throws ErrorException @SOCConstraint(2.1 == nrm)
        @test string(@SOCConstraint(2.1 ≥ nrm)) == "$Vert[w,-w + 1]$Vert$sub2 $leq 2.1"
        # 1-4 Number--AffExpr
        @test string(1.5 + aff) == "7.1 x + 4"
        @test string(1.5 - aff) == "-7.1 x - 1"
        @test string(2 * aff) == "14.2 x + 5"
        @test_throws ErrorException 2 / aff
        @test string(@LinearConstraint(1 ≤ aff)) == "-7.1 x $leq 1.5"
        @test string(@LinearConstraint(1 == aff)) == "-7.1 x $eq 1.5"
        @test string(@LinearConstraint(1 ≥ aff)) == "-7.1 x $geq 1.5"
        # 1-5 Number--QuadExpr
        @test string(1.5 + q) == "2.5 y*z + 7.1 x + 4"
        @test string(1.5 - q) == "-2.5 y*z - 7.1 x - 1"
        @test string(2 * q) == "5 y*z + 14.2 x + 5"
        @test_throws ErrorException 2 / q
        @test string(@QuadConstraint(1 ≤ q)) == "-2.5 y*z - 7.1 x - 1.5 $leq 0"
        @test string(@QuadConstraint(1 == q)) == "-2.5 y*z - 7.1 x - 1.5 $eq 0"
        @test string(@QuadConstraint(1 ≥ q)) == "-2.5 y*z - 7.1 x - 1.5 $geq 0"
        # 1-6 Number--SOCExpr
        @test string(1.5 + socexpr) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 0.5"
        @test string(1.5 - socexpr) == "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + 3.5"
        @test string(1.5 * socexpr) == "2.25 $Vert[w,-w + 1]$Vert$sub2 - 1.5 w - 3"
        @test_throws MethodError 1.5 / socexpr
        @test_throws ErrorException @SOCConstraint(2.1 ≤ socexpr)
        @test_throws ErrorException @SOCConstraint(2.1 == socexpr)
        @test string(@SOCConstraint(2.1 ≥ socexpr)) == "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 4.1"
        end

        # 2. Variable tests
        @testset "Variable--???" begin
        # 2-0 Variable unary
        @test (+x) === x
        @test string(-x) == "-x"
        # 2-1 Variable--Number
        @test string(w + 4.13) == "w + 4.13"
        @test string(w - 4.13) == "w - 4.13"
        @test string(w * 4.13) == "4.13 w"
        @test string(w / 2.00) == "0.5 w"
        @test w == w
        @test string(@LinearConstraint(w ≤ 1)) == "w $leq 1"
        @test string(@LinearConstraint(w == 1)) == "w $eq 1"
        @test string(@LinearConstraint(w ≥ 1)) == "w $geq 1"
        @test string(@QuadConstraint(x*y ≤ 1)) == "x*y - 1 $leq 0"
        @test string(@QuadConstraint(x*y == 1)) == "x*y - 1 $eq 0"
        @test string(@QuadConstraint(x*y ≥ 1)) == "x*y - 1 $geq 0"
        # 2-2 Variable--Variable
        @test string(w + x) == "w + x"
        @test string(w - x) == "w - x"
        @test string(w * x) == "w*x"
        @test string(x - x) == "0"
        @test_throws ErrorException w / x
        @test string(@LinearConstraint(w ≤ x)) == "w - x $leq 0"
        @test string(@LinearConstraint(w == x)) == "w - x $eq 0"
        @test string(@LinearConstraint(w ≥ x)) == "w - x $geq 0"
        @test string(@QuadConstraint(y*z ≤ x)) == "y*z - x $leq 0"
        @test string(@QuadConstraint(y*z == x)) == "y*z - x $eq 0"
        @test string(@QuadConstraint(y*z ≥ x)) == "y*z - x $geq 0"
        @test string(@LinearConstraint(x ≤ x)) == "0 $leq 0"
        @test string(@LinearConstraint(x == x)) == "0 $eq 0"
        @test string(@LinearConstraint(x ≥ x)) == "0 $geq 0"
        # 2-3 Variable--Norm
        @test string(w + nrm) == "$Vert[w,-w + 1]$Vert$sub2 + w"
        @test string(w - nrm) == "-1.0 $Vert[w,-w + 1]$Vert$sub2 + w"
        @test_throws MethodError w * nrm
        @test_throws MethodError w / nrm
        @test_throws ErrorException @SOCConstraint(w ≤ nrm)
        @test_throws ErrorException @SOCConstraint(w == nrm)
        @test string(@SOCConstraint(w ≥ nrm)) == "$Vert[w,-w + 1]$Vert$sub2 $leq w"
        # 2-4 Variable--AffExpr
        @test string(z + aff) == "7.1 x + z + 2.5"
        @test string(z - aff) == "-7.1 x + z - 2.5"
        @test string(z * aff) == "7.1 x*z + 2.5 z"
        @test_throws ErrorException z / aff
        @test_throws MethodError z ≤ aff
        @test string(@LinearConstraint(z ≤ aff)) ∈ ["z - 7.1 x $leq 2.5", "-7.1 x + z $leq 2.5"]
        @test string(@LinearConstraint(z == aff)) ∈ ["z - 7.1 x $eq 2.5", "-7.1 x + z $eq 2.5"]
        @test string(@LinearConstraint(z ≥ aff)) ∈ ["z - 7.1 x $geq 2.5", "-7.1 x + z $geq 2.5"]
        @test string(@LinearConstraint(7.1 * x - aff ≤ 0)) == "0 $leq 2.5"
        @test string(@LinearConstraint(7.1 * x - aff == 0)) == "0 $eq 2.5"
        @test string(@LinearConstraint(7.1 * x - aff ≥ 0)) == "0 $geq 2.5"
        # 2-5 Variable--QuadExpr
        @test string(w + q) == "2.5 y*z + 7.1 x + w + 2.5"
        @test string(w - q) == "-2.5 y*z - 7.1 x + w - 2.5"
        @test_throws ErrorException w*q
        @test_throws ErrorException w/q
        @test string(@QuadConstraint(w ≤ q)) ∈ ["-2.5 y*z + w - 7.1 x - 2.5 $leq 0", "-2.5 y*z - 7.1 x + w - 2.5 $leq 0"]
        @test string(@QuadConstraint(w == q)) ∈  ["-2.5 y*z + w - 7.1 x - 2.5 $eq 0", "-2.5 y*z - 7.1 x + w - 2.5 $eq 0"]
        @test string(@QuadConstraint(w ≥ q)) ∈  ["-2.5 y*z + w - 7.1 x - 2.5 $geq 0", "-2.5 y*z - 7.1 x + w - 2.5 $geq 0"]
        # 2-6 Variable--SOCExpr
        @test string(y + socexpr) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + y - 2"
        @test string(y - socexpr) == "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + y + 2"
        @test_throws MethodError y * socexpr
        @test_throws MethodError y / socexpr
        @test_throws ErrorException @SOCConstraint(y ≤ socexpr)
        @test_throws ErrorException @SOCConstraint(y == socexpr)
        @test string(@SOCConstraint(y ≥ socexpr)) ∈ ["1.5 $Vert[w,-w + 1]$Vert$sub2 $leq y + w + 2", "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + y + 2"]

        end

        # 3. Norm tests
        @testset "Norm--???" begin
        # 3-0 Norm unary
        @test string(+nrm) == "$Vert[w,-w + 1]$Vert$sub2"
        @test string(-nrm) == "-1.0 $Vert[w,-w + 1]$Vert$sub2"
        # 3-1 Norm--Number
        @test string(nrm + 1.5) == "$Vert[w,-w + 1]$Vert$sub2 + 1.5"
        @test string(nrm - 1.5) == "$Vert[w,-w + 1]$Vert$sub2 - 1.5"
        @test string(nrm * 1.5) == "1.5 $Vert[w,-w + 1]$Vert$sub2"
        @test string(nrm / 1.5) == "0.6666666666666666 $Vert[w,-w + 1]$Vert$sub2"
        @test string(@SOCConstraint(nrm ≤ 1.5)) == "$Vert[w,-w + 1]$Vert$sub2 $leq 1.5"
        @test_throws ErrorException @SOCConstraint(nrm == 1.5)
        @test_throws ErrorException @SOCConstraint(nrm ≥ 1.5)
        # 3-2 Norm--Variable
        @test string(nrm + w) == "$Vert[w,-w + 1]$Vert$sub2 + w"
        @test string(nrm - w) == "$Vert[w,-w + 1]$Vert$sub2 - w"
        @test_throws MethodError nrm * w
        @test_throws MethodError nrm / w
        @test string(@SOCConstraint(nrm ≤ w)) == "$Vert[w,-w + 1]$Vert$sub2 $leq w"
        @test_throws ErrorException @SOCConstraint(nrm == w)
        @test_throws ErrorException @SOCConstraint(nrm ≥ w)
        # 3-3 Norm--Norm
        @test_throws MethodError nrm + nrm
        @test_throws MethodError nrm - nrm
        @test_throws MethodError nrm * nrm
        @test_throws MethodError nrm / nrm
        @test_throws MethodError @SOCConstraint(nrm ≤ nrm)
        @test_throws MethodError @SOCConstraint(nrm == nrm)
        @test_throws MethodError @SOCConstraint(nrm ≥ nrm)
        # 3-4 Norm--AffExpr
        @test string(nrm + aff) == "$Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
        @test string(nrm - aff) == "$Vert[w,-w + 1]$Vert$sub2 - 7.1 x - 2.5"
        @test_throws MethodError nrm * aff
        @test_throws MethodError nrm / aff
        @test string(@SOCConstraint(nrm ≤ aff)) == "$Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + 2.5"
        @test_throws ErrorException @SOCConstraint(nrm == aff)
        @test_throws ErrorException @SOCConstraint(nrm ≥ aff)
        # 3-5 Norm--QuadExpr
        @test_throws MethodError nrm + q
        @test_throws MethodError nrm - q
        @test_throws MethodError nrm * q
        @test_throws MethodError nrm / q
        @test_throws MethodError @SOCConstraint(nrm ≤ q)
        @test_throws MethodError @SOCConstraint(nrm == q)
        @test_throws MethodError @SOCConstraint(nrm ≥ q)
        # 3-6 Norm--SOCExpr
        @test_throws MethodError nrm + socexpr
        @test_throws MethodError nrm - socexpr
        @test_throws MethodError nrm * socexpr
        @test_throws MethodError nrm / socexpr
        @test_throws MethodError @SOCConstraint(nrm ≤ socexpr)
        @test_throws MethodError @SOCConstraint(nrm == socexpr)
        @test_throws MethodError @SOCConstraint(nrm ≥ socexpr)
        end

        # 4. AffExpr tests
        @testset "AffExpr--???" begin
        # 4-0 AffExpr unary
        @test string(+aff) == "7.1 x + 2.5"
        @test string(-aff) == "-7.1 x - 2.5"
        # 4-1 AffExpr--Number
        @test string(aff + 1.5) == "7.1 x + 4"
        @test string(aff - 1.5) == "7.1 x + 1"
        @test string(aff * 2) == "14.2 x + 5"
        @test string(aff / 2) == "3.55 x + 1.25"
        @test_throws MethodError aff ≤ 1
        @test aff == aff
        @test_throws MethodError aff ≥ 1
        @test string(@LinearConstraint(aff ≤ 1)) == "7.1 x $leq -1.5"
        @test string(@LinearConstraint(aff == 1)) == "7.1 x $eq -1.5"
        @test string(@LinearConstraint(aff ≥ 1)) == "7.1 x $geq -1.5"
        # 4-2 AffExpr--Variable
        @test string(aff + z) == "7.1 x + z + 2.5"
        @test string(aff - z) == "7.1 x - z + 2.5"
        @test string(aff * z) == "7.1 x*z + 2.5 z"
        @test_throws ErrorException aff/z
        @test string(@LinearConstraint(aff ≤ z)) == "7.1 x - z $leq -2.5"
        @test string(@LinearConstraint(aff == z)) == "7.1 x - z $eq -2.5"
        @test string(@LinearConstraint(aff ≥ z)) == "7.1 x - z $geq -2.5"
        @test string(@LinearConstraint(aff - 7.1 * x ≤ 0)) == "0 $leq -2.5"
        @test string(@LinearConstraint(aff - 7.1 * x == 0)) == "0 $eq -2.5"
        @test string(@LinearConstraint(aff - 7.1 * x ≥ 0)) == "0 $geq -2.5"
        # 4-3 AffExpr--Norm
        @test string(aff + nrm) == "$Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
        @test string(aff - nrm) == "-1.0 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
        @test_throws MethodError aff * nrm
        @test_throws MethodError aff / nrm
        @test_throws ErrorException @SOCConstraint(aff ≤ nrm)
        @test_throws ErrorException @SOCConstraint(aff == nrm)
        @test string(@SOCConstraint(aff ≥ nrm)) == "$Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + 2.5"
        # 4-4 AffExpr--AffExpr
        @test string(aff + aff2) == "7.1 x + 1.2 y + 3.7"
        @test string(aff - aff2) == "7.1 x - 1.2 y + 1.3"
        @test string(aff * aff2) == "8.52 x*y + 3 y + 8.52 x + 3"
        @test string((x+x)*(x+3)) == string((x+3)*(x+x))  # Issue #288
        @test_throws ErrorException aff/aff2
        @test string(@LinearConstraint(aff ≤ aff2)) == "7.1 x - 1.2 y $leq -1.3"
        @test string(@LinearConstraint(aff == aff2)) == "7.1 x - 1.2 y $eq -1.3"
        @test string(@LinearConstraint(aff ≥ aff2)) == "7.1 x - 1.2 y $geq -1.3"
        @test string(@LinearConstraint(aff-aff ≤ 0)) == "0 $leq 0"
        @test string(@LinearConstraint(aff-aff == 0)) == "0 $eq 0"
        @test string(@LinearConstraint(aff-aff ≥ 0)) == "0 $geq 0"
        # 4-5 AffExpr--QuadExpr
        @test string(aff2 + q) == "2.5 y*z + 1.2 y + 7.1 x + 3.7"
        @test string(aff2 - q) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
        @test_throws ErrorException aff2 * q
        @test_throws ErrorException aff2 / q
        @test string(@QuadConstraint(aff2 ≤ q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $leq 0"
        @test string(@QuadConstraint(aff2 == q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $eq 0"
        @test string(@QuadConstraint(aff2 ≥ q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $geq 0"
        # 4-6 AffExpr--SOCExpr
        @test string(aff + socexpr) == "1.5 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x - w + 0.5"
        @test string(aff - socexpr) == "-1.5 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x + w + 4.5"
        @test_throws MethodError aff * socexpr
        @test_throws MethodError aff / socexpr
        @test_throws ErrorException @SOCConstraint(aff ≤ socexpr)
        @test_throws ErrorException @SOCConstraint(aff == socexpr)
        @test string(@SOCConstraint(aff ≥ socexpr)) == "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + w + 4.5"
        end

        # 5. QuadExpr
        @testset "QuadExpr--???" begin
        # 5-0 QuadExpr unary
        @test string(+q) == "2.5 y*z + 7.1 x + 2.5"
        @test string(-q) == "-2.5 y*z - 7.1 x - 2.5"
        # 5-1 QuadExpr--Number
        @test string(q + 1.5) == "2.5 y*z + 7.1 x + 4"
        @test string(q - 1.5) == "2.5 y*z + 7.1 x + 1"
        @test string(q * 2) == "5 y*z + 14.2 x + 5"
        @test string(q / 2) == "1.25 y*z + 3.55 x + 1.25"
        @test q == q
        @test string(@QuadConstraint(aff2 ≤ q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $leq 0"
        @test string(@QuadConstraint(aff2 == q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $eq 0"
        @test string(@QuadConstraint(aff2 ≥ q)) == "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $geq 0"
        # 5-2 QuadExpr--Variable
        @test string(q + w) == "2.5 y*z + 7.1 x + w + 2.5"
        @test string(q - w) == "2.5 y*z + 7.1 x - w + 2.5"
        @test_throws ErrorException q*w
        @test_throws ErrorException q/w
        @test string(@QuadConstraint(q ≤ w)) == "2.5 y*z + 7.1 x - w + 2.5 $leq 0"
        @test string(@QuadConstraint(q == w)) == "2.5 y*z + 7.1 x - w + 2.5 $eq 0"
        @test string(@QuadConstraint(q ≥ w)) == "2.5 y*z + 7.1 x - w + 2.5 $geq 0"
        # 5-3 QuadExpr--Norm
        @test_throws MethodError q + nrm
        @test_throws MethodError q - nrm
        @test_throws MethodError q * nrm
        @test_throws MethodError q / nrm
        @test_throws MethodError @SOCConstraint(q ≤ nrm)
        @test_throws MethodError @SOCConstraint(q == nrm)
        @test_throws MethodError @SOCConstraint(q ≥ nrm)
        # 5-4 QuadExpr--AffExpr
        @test string(q + aff2) == "2.5 y*z + 7.1 x + 1.2 y + 3.7"
        @test string(q - aff2) == "2.5 y*z + 7.1 x - 1.2 y + 1.3"
        @test_throws ErrorException q * aff2
        @test_throws ErrorException q / aff2
        @test string(@QuadConstraint(q ≤ aff2)) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
        @test string(@QuadConstraint(q == aff2)) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
        @test string(@QuadConstraint(q ≥ aff2)) == "2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
        # 5-5 QuadExpr--QuadExpr
        @test string(q + q2) == "8 x*z + 2.5 y*z + 7.1 x + 1.2 y + 3.7"
        @test string(q - q2) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3"
        @test_throws ErrorException q * q2
        @test_throws ErrorException q / q2
        @test string(@QuadConstraint(q ≤ q2)) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
        @test string(@QuadConstraint(q == q2)) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
        @test string(@QuadConstraint(q ≥ q2)) == "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
        # 4-6 QuadExpr--SOCExpr
        @test_throws MethodError q + socexpr
        @test_throws MethodError q - socexpr
        @test_throws MethodError q * socexpr
        @test_throws MethodError q / socexpr
        @test_throws MethodError @SOCConstraint(q ≤ socexpr)
        @test_throws MethodError @SOCConstraint(q == socexpr)
        @test_throws MethodError @SOCConstraint(q ≥ socexpr)
        end

        # 6. SOCExpr tests
        @testset "SOCExpr--???" begin
        # 6-0 SOCExpr unary
        @test string(+socexpr) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 2"
        @test string(-socexpr) == "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + 2"
        # 6-1 SOCExpr--Number
        @test string(socexpr + 1.5) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 0.5"
        @test string(socexpr - 1.5) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 3.5"
        @test string(socexpr * 1.5) == "2.25 $Vert[w,-w + 1]$Vert$sub2 - 1.5 w - 3"
        @test string(socexpr / 1.5) == "$Vert[w,-w + 1]$Vert$sub2 - 0.6666666666666666 w - 1.3333333333333333"
        @test string(@SOCConstraint(socexpr ≤ 1.5)) == "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 3.5"
        @test_throws ErrorException @SOCConstraint(socexpr == 1.5)
        @test_throws ErrorException @SOCConstraint(socexpr ≥ 1.5)
        # 6-2 SOCExpr--Variable
        @test string(socexpr + y) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + y - 2"
        @test string(socexpr - y) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - y - 2"
        @test_throws MethodError socexpr * y
        @test_throws MethodError socexpr / y
        @test string(@SOCConstraint(socexpr ≤ y)) == "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + y + 2"
        @test_throws ErrorException @SOCConstraint(socexpr == y)
        @test_throws ErrorException @SOCConstraint(socexpr ≥ y)
        # 6-3 SOCExpr--Norm
        @test_throws MethodError socexpr + nrm
        @test_throws MethodError socexpr - nrm
        @test_throws MethodError socexpr * nrm
        @test_throws MethodError socexpr / nrm
        @test_throws MethodError @SOCConstraint(socexpr ≤ nrm)
        @test_throws MethodError @SOCConstraint(socexpr == nrm)
        @test_throws MethodError @SOCConstraint(socexpr ≥ nrm)
        # 6-4 SOCExpr--AffExpr
        @test string(socexpr + aff) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + 7.1 x + 0.5"
        @test string(socexpr - aff) == "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 7.1 x - 4.5"
        @test_throws MethodError socexpr * aff
        @test_throws MethodError socexpr / aff
        @test string(@SOCConstraint(socexpr ≤ aff)) == "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 7.1 x + 4.5"
        @test_throws ErrorException @SOCConstraint(socexpr == aff)
        @test_throws ErrorException @SOCConstraint(socexpr ≥ aff)
        # 6-5 SOCExpr--QuadExpr
        @test_throws MethodError socexpr + q
        @test_throws MethodError socexpr - q
        @test_throws MethodError socexpr * q
        @test_throws MethodError socexpr / q
        @test_throws MethodError @SOCConstraint(socexpr ≤ q)
        @test_throws MethodError @SOCConstraint(socexpr == q)
        @test_throws MethodError @SOCConstraint(socexpr ≥ q)
        # 6-6 SOCExpr--SOCExpr
        @test_throws MethodError socexpr + socexpr
        @test_throws MethodError socexpr - socexpr
        @test_throws MethodError socexpr * socexpr
        @test_throws MethodError socexpr / socexpr
        @test_throws MethodError @SOCConstraint(socexpr ≤ socexpr)
        @test_throws MethodError @SOCConstraint(socexpr == socexpr)
        @test_throws MethodError @SOCConstraint(socexpr ≥ socexpr)
        end
    end

    @testset "Higher-level operators" begin
    @testset "sum" begin
        sum_m = Model()
        @variable(sum_m, 0 ≤ matrix[1:3,1:3] ≤ 1, start = 1)
        # sum(j::JuMPArray{Variable})
        @test string(sum(matrix)) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
        # sum(j::JuMPArray{Variable}) in a macro
        @objective(sum_m, Max, sum(matrix))
        @test string(sum_m.obj) == "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"

        # sum{T<:Real}(j::JuMPArray{T})
        @test isapprox(sum(getvalue(matrix)), 9)
        # sum(j::Array{Variable})
        @test string(sum(matrix[1:3,1:3])) == string(sum(matrix))
        # sum(affs::Array{AffExpr})
        @test string(sum([2*matrix[i,j] for i in 1:3, j in 1:3])) == "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"

        S = [1,3]
        @variable(sum_m, x[S], start=1)
        # sum(j::JuMPDict{Variable})
        @test length(string(sum(x))) == 11 # order depends on hashing
        @test occursin("x[1]",string(sum(x))) == true
        @test occursin("x[3]",string(sum(x))) == true
        # sum{T<:Real}(j::JuMPDict{T})
        @test sum(getvalue(x)) == 2
    end

    @testset "dot" begin
        dot_m = Model()
        @variable(dot_m, 0 ≤ x[1:3] ≤ 1)
        c = vcat(1:3)
        @test string(dot(c,x)) == "x[1] + 2 x[2] + 3 x[3]"
        @test string(dot(x,c)) == "x[1] + 2 x[2] + 3 x[3]"

        A = [1 3 ; 2 4]
        @variable(dot_m, 1 ≤ y[1:2,1:2] ≤ 1)
        @test string(dot(A,y)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
        @test string(dot(y,A)) == "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

        B = ones(2,2,2)
        @variable(dot_m, 0 ≤ z[1:2,1:2,1:2] ≤ 1)
        @test string(dot(B,z)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
        @test string(dot(z,B)) == "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

        @objective(dot_m, Max, dot(x, ones(3)) - dot(y, ones(2,2)) )
        #solve(dot_m)
        for i in 1:3
            setvalue(x[i], 1)
        end
        for i in 1:2, j in 1:2
            setvalue(y[i,j], 1)
        end
        @test isapprox(dot(convert(Array{Float64}, c), getvalue(x)), 6.0)
        @test isapprox(dot(convert(Array{Float64}, A), getvalue(y)), 10.0)

        # https://github.com/JuliaOpt/JuMP.jl/issues/656
        issue656 = Model()
        @variable(issue656, x)
        floats = Float64[i for i in 1:2]
        anys   = Array{Any}(undef, 2)
        anys[1] = 10
        anys[2] = 20 + x
        @test dot(floats, anys) == 10 + 40 + 2x

        # https://github.com/JuliaOpt/JuMP.jl/pull/943
        pull943 = Model()
        @variable(pull943, x[1 : 10^6]);
        setvalue.(x, 1 : 10^6)
        @expression(pull943, testsum, sum(x[i] * i for i = 1 : 10^6))
        @expression(pull943, testdot1, dot(x, 1 : 10^6))
        @expression(pull943, testdot2, dot(1 : 10^6, x))
        @test isapprox(getvalue(testsum), getvalue(testdot1))
        @test isapprox(getvalue(testsum), getvalue(testdot2))
    end
    end

    @testset "Vectorized operations" begin

    @testset "Transpose" begin
        m = Model()
        @variable(m, x[1:3])
        @variable(m, y[1:2,1:3])
        @variable(m, z[2:5])
        @test vec_eq(x', [x[1] x[2] x[3]])
        @test vec_eq(transpose(x), [x[1] x[2] x[3]])
        @test vec_eq(y', [y[1,1] y[2,1]
                          y[1,2] y[2,2]
                          y[1,3] y[2,3]])
        @test vec_eq(transpose(y),
                    [y[1,1] y[2,1]
                     y[1,2] y[2,2]
                     y[1,3] y[2,3]])
        @test_throws ErrorException z'
        @test_throws ErrorException transpose(z)
    end

    @testset "Vectorized arithmetic" begin
        m = Model()
        @variable(m, x[1:3])
        A = [2 1 0
             1 2 1
             0 1 2]
        B = sparse(A)
        @variable(m, X11)
        @variable(m, X23)
        X = sparse([1, 2], [1, 3], [X11, X23], 3, 3) # for testing Variable
        @variable(m, Xd[1:3, 1:3])
        Y = sparse([1, 2], [1, 3], [2X11, 4X23], 3, 3) # for testing GenericAffExpr
        Yd = [2X11 0    0
              0    0 4X23
              0    0    0]
        Z = sparse([1, 2], [1, 3], [X11^2, 2X23^2], 3, 3) # for testing GenericQuadExpr
        Zd = [X11^2 0      0
              0     0 2X23^2
              0     0      0]
        v = [4, 5, 6]

        @testset "Sum of matrices" begin
            @test vec_eq(JuMP.@Expression(Xd + Yd),     Xd + Yd)
            @test vec_eq(JuMP.@Expression(Xd + 2Yd),    Xd + 2Yd)
            @test vec_eq(JuMP.@Expression(Xd + Yd * 2), Xd + Yd * 2)
            @test vec_eq(JuMP.@Expression(Yd + Xd),     Yd + Xd)
            @test vec_eq(JuMP.@Expression(Yd + 2Xd),    Yd + 2Xd)
            @test vec_eq(JuMP.@Expression(Yd + Xd * 2), Yd + Xd * 2)
            @test vec_eq(JuMP.@Expression(Yd + Zd),     Yd + Zd)
            @test vec_eq(JuMP.@Expression(Yd + 2Zd),    Yd + 2Zd)
            @test vec_eq(JuMP.@Expression(Yd + Zd * 2), Yd + Zd * 2)
            @test vec_eq(JuMP.@Expression(Zd + Yd),     Zd + Yd)
            @test vec_eq(JuMP.@Expression(Zd + 2Yd),    Zd + 2Yd)
            @test vec_eq(JuMP.@Expression(Zd + Yd * 2), Zd + Yd * 2)
            @test vec_eq(JuMP.@Expression(Zd + Xd),     Zd + Xd)
            @test vec_eq(JuMP.@Expression(Zd + 2Xd),    Zd + 2Xd)
            @test vec_eq(JuMP.@Expression(Zd + Xd * 2), Zd + Xd * 2)
            @test vec_eq(JuMP.@Expression(Xd + Zd),     Xd + Zd)
            @test vec_eq(JuMP.@Expression(Xd + 2Zd),    Xd + 2Zd)
            @test vec_eq(JuMP.@Expression(Xd + Zd * 2), Xd + Zd * 2)
        end

        @test vec_eq(A*x, [2x[1] +  x[2]
                            2x[2] + x[1] + x[3]
                            x[2] + 2x[3]])
        @test vec_eq(A*x, B*x)
        @test vec_eq(A*x, @JuMP.Expression(B*x))
        @test vec_eq(@JuMP.Expression(A*x), @JuMP.Expression(B*x))
        @test vec_eq(x'*A, [2x[1]+x[2]; 2x[2]+x[1]+x[3]; x[2]+2x[3]]')
        @test vec_eq(x'*A, x'*B)
        @test vec_eq(x'*A, @JuMP.Expression(x'*B))
        @test vec_eq(@JuMP.Expression(x'*A), @JuMP.Expression(x'*B))
        @test vec_eq(x'*A*x, [2x[1]*x[1] + 2x[1]*x[2] + 2x[2]*x[2] + 2x[2]*x[3] + 2x[3]*x[3]])
        @test vec_eq(x'A*x, x'*B*x)
        @test vec_eq(x'*A*x, @JuMP.Expression(x'*B*x))
        @test vec_eq(@JuMP.Expression(x'*A*x), @JuMP.Expression(x'*B*x))

        y = A*x
        @test vec_eq(-x, [-x[1], -x[2], -x[3]])
        @test vec_eq(-y, [-2x[1] -  x[2]
                           -x[1] - 2x[2] -  x[3]
                                   -x[2] - 2x[3]])
        @test vec_eq(y + 1, [2x[1] +  x[2]         + 1
                              x[1] + 2x[2] +  x[3] + 1
                                      x[2] + 2x[3] + 1])
        @test vec_eq(y - 1, [2x[1] +  x[2]         - 1
                              x[1] + 2x[2] +  x[3] - 1
                                      x[2] + 2x[3] - 1])
        @test vec_eq(y + 2ones(3), [2x[1] +  x[2]         + 2
                                     x[1] + 2x[2] +  x[3] + 2
                                             x[2] + 2x[3] + 2])
        @test vec_eq(y - 2ones(3), [2x[1] +  x[2]         - 2
                                     x[1] + 2x[2] +  x[3] - 2
                                             x[2] + 2x[3] - 2])
        @test vec_eq(2ones(3) + y, [2x[1] +  x[2]         + 2
                                     x[1] + 2x[2] +  x[3] + 2
                                             x[2] + 2x[3] + 2])
        @test vec_eq(2ones(3) - y, [-2x[1] -  x[2]         + 2
                                     -x[1] - 2x[2] -  x[3] + 2
                                             -x[2] - 2x[3] + 2])
        @test vec_eq(y + x, [3x[1] +  x[2]
                              x[1] + 3x[2] +  x[3]
                                      x[2] + 3x[3]])
        @test vec_eq(x + y, [3x[1] +  x[2]
                              x[1] + 3x[2] +  x[3]
                                      x[2] + 3x[3]])
        @test vec_eq(2y + 2x, [6x[1] + 2x[2]
                               2x[1] + 6x[2] + 2x[3]
                                       2x[2] + 6x[3]])
        @test vec_eq(y - x, [ x[1] + x[2]
                              x[1] + x[2] + x[3]
                                     x[2] + x[3]])
        @test vec_eq(x - y, [-x[1] - x[2]
                             -x[1] - x[2] - x[3]
                                    -x[2] - x[3]])
        @test vec_eq(y + x[:], [3x[1] +  x[2]
                                 x[1] + 3x[2] +  x[3]
                                         x[2] + 3x[3]])
        @test vec_eq(x[:] + y, [3x[1] +  x[2]
                                 x[1] + 3x[2] +  x[3]
                                         x[2] + 3x[3]])

        @test vec_eq(@JuMP.Expression(A*x/2), A*x/2)
        @test vec_eq(X*v,  [4X11; 6X23; 0])
        @test vec_eq(v'*X,  [4X11  0   5X23])
        @test vec_eq(transpose(v)*X, [4X11  0   5X23])
        @test vec_eq(X'*v,  [4X11;  0;  5X23])
        @test vec_eq(transpose(X)*v, [4X11; 0;  5X23])
        @test vec_eq(X*A,  [2X11  X11  0
                            0     X23  2X23
                            0     0    0   ])
        @test vec_eq(A*X,  [2X11  0    X23
                            X11   0    2X23
                            0     0    X23])
        @test vec_eq(A*X', [2X11  0    0
                            X11   X23  0
                            0     2X23 0])
        @test vec_eq(X'*A, [2X11  X11  0
                            0     0    0
                            X23   2X23 X23])
        @test vec_eq(transpose(X)*A, [2X11 X11  0
                             0    0    0
                             X23  2X23 X23])
        @test vec_eq(A'*X, [2X11  0 X23
                            X11   0 2X23
                            0     0 X23])
        @test vec_eq(transpose(X)*A, X'*A)
        @test vec_eq(transpose(A)*X, A'*X)
        @test vec_eq(X*A, X*B)
        @test vec_eq(Y'*A, transpose(Y)*A)
        @test vec_eq(A*Y', A*transpose(Y))
        @test vec_eq(Z'*A, transpose(Z)*A)
        @test vec_eq(Xd'*Y, transpose(Xd)*Y)
        @test vec_eq(Y'*Xd, transpose(Y)*Xd)
        @test vec_eq(Xd'*Xd, transpose(Xd)*Xd)
        # @test_broken vec_eq(A*X, B*X)
        # @test_broken vec_eq(A*X', B*X')
        @test vec_eq(X'*A, X'*B)
        # @test_broken(X'*X, transpose(X)*X) # sparse quadratic known to be broken, see #912
    end

    @testset "Dot-ops" begin
        m = Model()
        @variable(m, x[1:2,1:2])
        A = [1 2;
             3 4]
        B = sparse(A)
        y = SparseMatrixCSC(2, 2, copy(B.colptr), copy(B.rowval), vec(x))
        @test vec_eq(A.+x, [1+x[1,1]  2+x[1,2];
                                       3+x[2,1]  4+x[2,2]])
        @test vec_eq(A.+x, B.+x)
        # @test vec_eq(A.+x, A.+y) == true
        # @test vec_eq(A.+y, B.+y) == true
        @test vec_eq(x.+A, [1+x[1,1]  2+x[1,2];
                                       3+x[2,1]  4+x[2,2]])
        @test vec_eq(x.+A, x.+B) == true
        @test vec_eq(x.+A, y.+A)
        @test vec_eq(x .+ x, [2x[1,1] 2x[1,2]; 2x[2,1] 2x[2,2]])
        # @test vec_eq(y.+A, y.+B) == true
        @test vec_eq(A.-x, [1-x[1,1]  2-x[1,2];
                                       3-x[2,1]  4-x[2,2]])
        @test vec_eq(A.-x, B.-x)
        @test vec_eq(A.-x, A.-y)
        @test vec_eq(x .- x, [zero(AffExpr) for _1 in 1:2, _2 in 1:2])
        # @test vec_eq(A.-y, B.-y) == true
        @test vec_eq(x.-A, [-1+x[1,1]  -2+x[1,2];
                                       -3+x[2,1]  -4+x[2,2]])
        @test vec_eq(x.-A, x.-B)
        @test vec_eq(x.-A, y.-A)
        # @test vec_eq(y.-A, y.-B) == true
        @test vec_eq(A.*x, [1*x[1,1]  2*x[1,2];
                                       3*x[2,1]  4*x[2,2]])
        @test vec_eq(A.*x, B.*x)
        @test vec_eq(A.*x, A.*y)
        # @test vec_eq(A.*y, B.*y) == true

        @test vec_eq(x.*A, [1*x[1,1]  2*x[1,2];
                                       3*x[2,1]  4*x[2,2]])
        @test vec_eq(x.*A, x.*B)
        @test vec_eq(x.*A, y.*A)
        # @test vec_eq(y.*A, y.*B) == true

        @test vec_eq(x .* x, [x[1,1]^2 x[1,2]^2; x[2,1]^2 x[2,2]^2])
        @test_throws ErrorException vec_eq(A./x, [1*x[1,1]  2*x[1,2];
                                              3*x[2,1]  4*x[2,2]])
        @test vec_eq(x./A, [1/1*x[1,1]  1/2*x[1,2];
                                       1/3*x[2,1]  1/4*x[2,2]])
        @test vec_eq(x./A, x./B)
        @test vec_eq(x./A, y./A)
        # @test vec_eq(A./y, B./y) == true

        @test vec_eq((2*x) / 3, Array((2*y) / 3))
        @test vec_eq(2 * (x/3), Array(2 * (y/3)))
        @test vec_eq(x[1,1] * A, Array(x[1,1] * B))
    end

    @testset "Vectorized comparisons" begin
        m = Model()
        @variable(m, x[1:3])
        A = [1 2 3
             0 4 5
             6 0 7]
        B = sparse(A)
        # force vector output
        @constraint(m, reshape(x,(1,3))*A*x .>= 1)
        @test vec_eq(m.quadconstr[1].terms, [x[1]*x[1] + 2x[1]*x[2] + 4x[2]*x[2] + 9x[1]*x[3] + 5x[2]*x[3] + 7x[3]*x[3] - 1])
        @test m.quadconstr[1].sense == :(>=)
        @constraint(m, x'*A*x >= 1)
        @test vec_eq(m.quadconstr[1].terms, m.quadconstr[2].terms)

        mat = [ 3x[1] + 12x[3] +  4x[2]
                2x[1] + 12x[2] + 10x[3]
               15x[1] +  5x[2] + 21x[3]]

        @constraint(m, (x'A)' .+ 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[1:3])
        lbs   = map(v->v.lb,    m.linconstr[1:3])
        ubs   = map(v->v.ub,    m.linconstr[1:3])
        @test vec_eq(terms, mat)
        @test lbs == fill(-Inf, 3)
        @test ubs == fill(   1, 3)
        @test vec_eq((x'A)' .+ 2A*x, (x'A)' .+ 2B*x)
        @test vec_eq((x'A)' .+ 2A*x, (x'B)' .+ 2A*x)
        @test vec_eq((x'A)' .+ 2A*x, (x'B)' .+ 2B*x)
        @test vec_eq((x'A)' .+ 2A*x, @JuMP.Expression((x'A)' .+ 2A*x))
        @test vec_eq((x'A)' .+ 2A*x, @JuMP.Expression((x'B)' .+ 2A*x))
        @test vec_eq((x'A)' .+ 2A*x, @JuMP.Expression((x'A)' .+ 2B*x))
        @test vec_eq((x'A)' .+ 2A*x, @JuMP.Expression((x'B)' .+ 2B*x))

        @constraint(m, -1 .<= (x'A)' .+ 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[4:6])
        lbs   = map(v->v.lb,    m.linconstr[4:6])
        ubs   = map(v->v.ub,    m.linconstr[4:6])
        @test vec_eq(terms, mat) == true
        @test lbs == fill(-1, 3)
        @test ubs == fill( 1, 3)

        @constraint(m, -[1:3;] .<= (x'A)' .+ 2A*x .<= 1)
        terms = map(v->v.terms, m.linconstr[7:9])
        lbs   = map(v->v.lb,    m.linconstr[7:9])
        ubs   = map(v->v.ub,    m.linconstr[7:9])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == fill( 1, 3)

        @constraint(m, -[1:3;] .<= (x'A)' .+ 2A*x .<= [3:-1:1;])
        terms = map(v->v.terms, m.linconstr[10:12])
        lbs   = map(v->v.lb,    m.linconstr[10:12])
        ubs   = map(v->v.ub,    m.linconstr[10:12])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == [3:-1:1;]

        @constraint(m, -[1:3;] .<= (x'A)' .+ 2A*x .<= 3)
        terms = map(v->v.terms, m.linconstr[13:15])
        lbs   = map(v->v.lb,    m.linconstr[13:15])
        ubs   = map(v->v.ub,    m.linconstr[13:15])
        @test vec_eq(terms, mat) == true
        @test lbs == -[1:3;]
        @test ubs == fill(3,3)
    end

    end

    @testset "JuMPArray concatenation" begin
        m = Model()
        @variable(m, x[1:3])
        @variable(m, y[1:3,1:3])
        @variable(m, z[1:1])
        @variable(m, w[1:1,1:3])

        @test vec_eq([x y], [x[1] y[1,1] y[1,2] y[1,3]
                                        x[2] y[2,1] y[2,2] y[2,3]
                                        x[3] y[3,1] y[3,2] y[3,3]])

        @test vec_eq([x 2y+1], [x[1] 2y[1,1]+1 2y[1,2]+1 2y[1,3]+1
                                           x[2] 2y[2,1]+1 2y[2,2]+1 2y[2,3]+1
                                           x[3] 2y[3,1]+1 2y[3,2]+1 2y[3,3]+1])

        @test vec_eq([1 x'], [1 x[1] x[2] x[3]])
        @test vec_eq([2x;x], [2x[1],2x[2],2x[3],x[1],x[2],x[3]])
        # vcat on JuMPArray
        @test vec_eq([x;x], [x[1],x[2],x[3],x[1],x[2],x[3]])
        # hcat on JuMPArray
        @test vec_eq([x x], [x[1] x[1]
                                        x[2] x[2]
                                        x[3] x[3]])
        # hvcat on JuMPArray
        tmp1 = [z w; x y]
        tmp2 = [z[1] w[1,1] w[1,2] w[1,3]
                x[1] y[1,1] y[1,2] y[1,3]
                x[2] y[2,1] y[2,2] y[2,3]
                x[3] y[3,1] y[3,2] y[3,3]]
        @test vec_eq(tmp1, tmp2)
        tmp3 = [1 2x'
                x 2y-x*x']
        tmp4 = [1    2x[1]               2x[2]               2x[3]
                x[1] -x[1]*x[1]+2y[1,1]  -x[1]*x[2]+2y[1,2]  -x[1]*x[3] + 2y[1,3]
                x[2] -x[1]*x[2]+2y[2,1]  -x[2]*x[2]+2y[2,2]  -x[2]*x[3] + 2y[2,3]
                x[3] -x[1]*x[3]+2y[3,1]  -x[2]*x[3]+2y[3,2]  -x[3]*x[3] + 2y[3,3]]
        @test vec_eq(tmp3, tmp4)

        A = sprand(3, 3, 0.2)
        B = Array(A)
        @test vec_eq([A y], [B y])
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
            #@test sum(x) == sum(x2)  # TODO stackoverflow in offset arrays
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
                #@test diagm(x) == diagm(x2)  # TODO fix OffsetArrays
            end
            #@test norm(x).terms == norm(x2).terms  # TODO stackoverflow in offset arrays
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
        if VERSION ≤ v"0.7-"
            @test typeof(z) == (isdefined(Base, :RowVector) ? RowVector{ElemT, ConjArray{ElemT, 1, Vector{ElemT}}} : Matrix{ElemT})
        else
            @test typeof(z) == Adjoint{ElemT,Vector{ElemT}}
        end
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
