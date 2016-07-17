#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/operator.jl
# Testing operator overloading is correct
#############################################################################
using JuMP, FactCheck

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const  leq = JuMP.repl[:leq]
const  geq = JuMP.repl[:geq]
const   eq = JuMP.repl[:eq]
const Vert = JuMP.repl[:Vert]
const sub2 = JuMP.repl[:sub2]

facts("[operator] Testing basic operator overloads") do
    m = Model()
    @variable(m, w)
    @variable(m, x)
    @variable(m, y)
    @variable(m, z)
    aff = 7.1 * x + 2.5
    @fact string(aff) --> "7.1 x + 2.5"
    aff2 = 1.2 * y + 1.2
    @fact string(aff2) --> "1.2 y + 1.2"
    q = 2.5 * y * z + aff
    @fact string(q) --> "2.5 y*z + 7.1 x + 2.5"
    q2 = 8 * x * z + aff2
    @fact string(q2) --> "8 x*z + 1.2 y + 1.2"
    q3 = 2 * x * x + 1 * y * y + z + 3
    @fact string(q3) --> "2 x² + y² + z + 3"

    nrm = norm([w,1-w])
    @fact string(nrm) --> "$Vert[w,-w + 1]$Vert$sub2"
    socexpr = 1.5*nrm - 2 - w
    @fact string(socexpr) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 2"

    @fact isequal(3w + 2y, 3w +2y) --> true
    @fact isequal(3w + 2y + 1, 3w + 2y) --> false

    # Different objects that must all interact:
    # 1. Number
    # 2. Variable
    # 3. AffExpr
    # 4. QuadExpr

    # 1. Number tests
    context("Number--???") do
    # 1-1 Number--Number - nope!
    # 1-2 Number--Variable
    @fact string(4.13 + w) --> "w + 4.13"
    @fact string(3.16 - w) --> "-w + 3.16"
    @fact string(5.23 * w) --> "5.23 w"
    @fact_throws  2.94 / w
    @fact string(@LinearConstraint(2.1 ≤ w)) --> "-w $leq -2.1"
    @fact string(@LinearConstraint(2.1 == w)) --> "-w $eq -2.1"
    @fact string(@LinearConstraint(2.1 ≥ w)) --> "-w $geq -2.1"
    # 1-3 Number--Norm
    @fact string(4.13 + nrm) --> "$Vert[w,-w + 1]$Vert$sub2 + 4.13"
    @fact string(3.16 - nrm) --> "-1.0 $Vert[w,-w + 1]$Vert$sub2 + 3.16"
    @fact string(5.23 * nrm) --> "5.23 $Vert[w,-w + 1]$Vert$sub2"
    @fact_throws 2.94 / nrm
    @fact_throws @SOCConstraint(2.1 ≤ nrm)
    @fact_throws @SOCConstraint(2.1 == nrm)
    @fact string(@SOCConstraint(2.1 ≥ nrm)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq 2.1"
    # 1-4 Number--AffExpr
    @fact string(1.5 + aff) --> "7.1 x + 4"
    @fact string(1.5 - aff) --> "-7.1 x - 1"
    @fact string(2 * aff) --> "14.2 x + 5"
    @fact_throws  2 / aff
    @fact string(@LinearConstraint(1 ≤ aff)) --> "-7.1 x $leq 1.5"
    @fact string(@LinearConstraint(1 == aff)) --> "-7.1 x $eq 1.5"
    @fact string(@LinearConstraint(1 ≥ aff)) --> "-7.1 x $geq 1.5"
    # 1-5 Number--QuadExpr
    @fact string(1.5 + q) --> "2.5 y*z + 7.1 x + 4"
    @fact string(1.5 - q) --> "-2.5 y*z - 7.1 x - 1"
    @fact string(2 * q) --> "5 y*z + 14.2 x + 5"
    @fact_throws  2 / q
    @fact string(@QuadConstraint(1 ≤ q)) --> "-2.5 y*z - 7.1 x - 1.5 $leq 0"
    @fact string(@QuadConstraint(1 == q)) --> "-2.5 y*z - 7.1 x - 1.5 $eq 0"
    @fact string(@QuadConstraint(1 ≥ q)) --> "-2.5 y*z - 7.1 x - 1.5 $geq 0"
    # 1-6 Number--SOCExpr
    @fact string(1.5 + socexpr) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 0.5"
    @fact string(1.5 - socexpr) --> "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + 3.5"
    @fact string(1.5 * socexpr) --> "2.25 $Vert[w,-w + 1]$Vert$sub2 - 1.5 w - 3"
    @fact_throws 1.5 / socexpr
    @fact_throws @SOCConstraint(2.1 ≤ socexpr)
    @fact_throws @SOCConstraint(2.1 == socexpr)
    @fact string(@SOCConstraint(2.1 ≥ socexpr)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 4.1"
    end

    # 2. Variable tests
    context("Variable--???") do
    # 2-0 Variable unary
    @fact (+x) --> exactly(x)
    @fact string(-x) --> "-x"
    # 2-1 Variable--Number
    @fact string(w + 4.13) --> "w + 4.13"
    @fact string(w - 4.13) --> "w - 4.13"
    @fact string(w * 4.13) --> "4.13 w"
    @fact string(w / 2.00) --> "0.5 w"
    @fact w == w --> true
    @fact string(@LinearConstraint(w ≤ 1)) --> "w $leq 1"
    @fact string(@LinearConstraint(w == 1)) --> "w $eq 1"
    @fact string(@LinearConstraint(w ≥ 1)) --> "w $geq 1"
    @fact string(@QuadConstraint(x*y ≤ 1)) --> "x*y - 1 $leq 0"
    @fact string(@QuadConstraint(x*y == 1)) --> "x*y - 1 $eq 0"
    @fact string(@QuadConstraint(x*y ≥ 1)) --> "x*y - 1 $geq 0"
    # 2-2 Variable--Variable
    @fact string(w + x) --> "w + x"
    @fact string(w - x) --> "w - x"
    @fact string(w * x) --> "w*x"
    @fact string(x - x) --> "0"
    @fact_throws  w / x
    @fact string(@LinearConstraint(w ≤ x)) --> "w - x $leq 0"
    @fact string(@LinearConstraint(w == x)) --> "w - x $eq 0"
    @fact string(@LinearConstraint(w ≥ x)) --> "w - x $geq 0"
    @fact string(@QuadConstraint(y*z ≤ x)) --> "y*z - x $leq 0"
    @fact string(@QuadConstraint(y*z == x)) --> "y*z - x $eq 0"
    @fact string(@QuadConstraint(y*z ≥ x)) --> "y*z - x $geq 0"
    @fact string(@LinearConstraint(x ≤ x)) --> "0 $leq 0"
    @fact string(@LinearConstraint(x == x)) --> "0 $eq 0"
    @fact string(@LinearConstraint(x ≥ x)) --> "0 $geq 0"
    # 2-3 Variable--Norm
    @fact string(w + nrm) --> "$Vert[w,-w + 1]$Vert$sub2 + w"
    @fact string(w - nrm) --> "-1.0 $Vert[w,-w + 1]$Vert$sub2 + w"
    @fact_throws w * nrm
    @fact_throws w / nrm
    @fact_throws @SOCConstraint(w ≤ nrm)
    @fact_throws @SOCConstraint(w == nrm)
    @fact string(@SOCConstraint(w ≥ nrm)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq w"
    # 2-4 Variable--AffExpr
    @fact string(z + aff) --> "7.1 x + z + 2.5"
    @fact string(z - aff) --> "-7.1 x + z - 2.5"
    @fact string(z * aff) --> "7.1 x*z + 2.5 z"
    @fact_throws  z / aff
    @fact_throws  z ≤ aff
    @fact string(@LinearConstraint(z ≤ aff)) --> "z - 7.1 x $leq 2.5"
    @fact string(@LinearConstraint(z == aff)) --> "z - 7.1 x $eq 2.5"
    @fact string(@LinearConstraint(z ≥ aff)) --> "z - 7.1 x $geq 2.5"
    @fact string(@LinearConstraint(7.1 * x - aff ≤ 0)) --> "0 $leq 2.5"
    @fact string(@LinearConstraint(7.1 * x - aff == 0)) --> "0 $eq 2.5"
    @fact string(@LinearConstraint(7.1 * x - aff ≥ 0)) --> "0 $geq 2.5"
    # 2-5 Variable--QuadExpr
    @fact string(w + q) --> "2.5 y*z + 7.1 x + w + 2.5"
    @fact string(w - q) --> "-2.5 y*z - 7.1 x + w - 2.5"
    @fact_throws  w*q
    @fact_throws  w/q
    @fact string(@QuadConstraint(w ≤ q)) --> "-2.5 y*z + w - 7.1 x - 2.5 $leq 0"
    @fact string(@QuadConstraint(w == q)) --> "-2.5 y*z + w - 7.1 x - 2.5 $eq 0"
    @fact string(@QuadConstraint(w ≥ q)) --> "-2.5 y*z + w - 7.1 x - 2.5 $geq 0"
    # 2-6 Variable--SOCExpr
    @fact string(y + socexpr) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + y - 2"
    @fact string(y - socexpr) --> "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + y + 2"
    @fact_throws y * socexpr
    @fact_throws y / socexpr
    @fact_throws @SOCConstraint(y ≤ socexpr)
    @fact_throws @SOCConstraint(y == socexpr)
    @fact string(@SOCConstraint(y ≥ socexpr)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq y + w + 2"
    end

    # 3. Norm tests
    context("Norm--???") do
    # 3-0 Norm unary
    @fact string(+nrm) --> "$Vert[w,-w + 1]$Vert$sub2"
    @fact string(-nrm) --> "-1.0 $Vert[w,-w + 1]$Vert$sub2"
    # 3-1 Norm--Number
    @fact string(nrm + 1.5) --> "$Vert[w,-w + 1]$Vert$sub2 + 1.5"
    @fact string(nrm - 1.5) --> "$Vert[w,-w + 1]$Vert$sub2 - 1.5"
    @fact string(nrm * 1.5) --> "1.5 $Vert[w,-w + 1]$Vert$sub2"
    @fact string(nrm / 1.5) --> "0.6666666666666666 $Vert[w,-w + 1]$Vert$sub2"
    @fact string(@SOCConstraint(nrm ≤ 1.5)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq 1.5"
    @fact_throws @SOCConstraint(nrm == 1.5)
    @fact_throws @SOCConstraint(nrm ≥ 1.5)
    # 3-2 Norm--Variable
    @fact string(nrm + w) --> "$Vert[w,-w + 1]$Vert$sub2 + w"
    @fact string(nrm - w) --> "$Vert[w,-w + 1]$Vert$sub2 - w"
    @fact_throws nrm * w
    @fact_throws nrm / w
    @fact string(@SOCConstraint(nrm ≤ w)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq w"
    @fact_throws @SOCConstraint(nrm == w)
    @fact_throws @SOCConstraint(nrm ≥ w)
    # 3-3 Norm--Norm
    @fact_throws nrm + nrm
    @fact_throws nrm - nrm
    @fact_throws nrm * nrm
    @fact_throws nrm / nrm
    @fact_throws @SOCConstraint(nrm ≤ nrm)
    @fact_throws @SOCConstraint(nrm == nrm)
    @fact_throws @SOCConstraint(nrm ≥ nrm)
    # 3-4 Norm--AffExpr
    @fact string(nrm + aff) --> "$Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
    @fact string(nrm - aff) --> "$Vert[w,-w + 1]$Vert$sub2 - 7.1 x - 2.5"
    @fact_throws nrm * aff
    @fact_throws nrm / aff
    @fact string(@SOCConstraint(nrm ≤ aff)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + 2.5"
    @fact_throws @SOCConstraint(nrm == aff)
    @fact_throws @SOCConstraint(nrm ≥ aff)
    # 3-5 Norm--QuadExpr
    @fact_throws nrm + q
    @fact_throws nrm - q
    @fact_throws nrm * q
    @fact_throws nrm / q
    @fact_throws @SOCConstraint(nrm ≤ q)
    @fact_throws @SOCConstraint(nrm == q)
    @fact_throws @SOCConstraint(nrm ≥ q)
    # 3-6 Norm--SOCExpr
    @fact_throws nrm + socexpr
    @fact_throws nrm - socexpr
    @fact_throws nrm * socexpr
    @fact_throws nrm / socexpr
    @fact_throws @SOCConstraint(nrm ≤ socexpr)
    @fact_throws @SOCConstraint(nrm == socexpr)
    @fact_throws @SOCConstraint(nrm ≥ socexpr)
    end

    # 4. AffExpr tests
    context("AffExpr--???") do
    # 4-0 AffExpr unary
    @fact string(+aff) --> "7.1 x + 2.5"
    @fact string(-aff) --> "-7.1 x - 2.5"
    # 4-1 AffExpr--Number
    @fact string(aff + 1.5) --> "7.1 x + 4"
    @fact string(aff - 1.5) --> "7.1 x + 1"
    @fact string(aff * 2) --> "14.2 x + 5"
    @fact string(aff / 2) --> "3.55 x + 1.25"
    @fact_throws aff ≤ 1
    @fact aff == aff --> true
    @fact_throws aff ≥ 1
    @fact string(@LinearConstraint(aff ≤ 1)) --> "7.1 x $leq -1.5"
    @fact string(@LinearConstraint(aff == 1)) --> "7.1 x $eq -1.5"
    @fact string(@LinearConstraint(aff ≥ 1)) --> "7.1 x $geq -1.5"
    # 4-2 AffExpr--Variable
    @fact string(aff + z) --> "7.1 x + z + 2.5"
    @fact string(aff - z) --> "7.1 x - z + 2.5"
    @fact string(aff * z) --> "7.1 x*z + 2.5 z"
    @fact_throws  aff/z
    @fact string(@LinearConstraint(aff ≤ z)) --> "7.1 x - z $leq -2.5"
    @fact string(@LinearConstraint(aff == z)) --> "7.1 x - z $eq -2.5"
    @fact string(@LinearConstraint(aff ≥ z)) --> "7.1 x - z $geq -2.5"
    @fact string(@LinearConstraint(aff - 7.1 * x ≤ 0)) --> "0 $leq -2.5"
    @fact string(@LinearConstraint(aff - 7.1 * x == 0)) --> "0 $eq -2.5"
    @fact string(@LinearConstraint(aff - 7.1 * x ≥ 0)) --> "0 $geq -2.5"
    # 4-3 AffExpr--Norm
    @fact string(aff + nrm) --> "$Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
    @fact string(aff - nrm) --> "-1.0 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x + 2.5"
    @fact_throws aff * nrm
    @fact_throws aff / nrm
    @fact_throws @SOCConstraint(aff ≤ nrm)
    @fact_throws @SOCConstraint(aff == nrm)
    @fact string(@SOCConstraint(aff ≥ nrm)) --> "$Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + 2.5"
    # 4-4 AffExpr--AffExpr
    @fact string(aff + aff2) --> "7.1 x + 1.2 y + 3.7"
    @fact string(aff - aff2) --> "7.1 x - 1.2 y + 1.3"
    @fact string(aff * aff2) --> "8.52 x*y + 3 y + 8.52 x + 3"
    @fact string((x+x)*(x+3)) --> string((x+3)*(x+x))  # Issue #288
    @fact_throws  aff/aff2
    @fact string(@LinearConstraint(aff ≤ aff2)) --> "7.1 x - 1.2 y $leq -1.3"
    @fact string(@LinearConstraint(aff == aff2)) --> "7.1 x - 1.2 y $eq -1.3"
    @fact string(@LinearConstraint(aff ≥ aff2)) --> "7.1 x - 1.2 y $geq -1.3"
    @fact string(@LinearConstraint(aff-aff ≤ 0)) --> "0 $leq 0"
    @fact string(@LinearConstraint(aff-aff == 0)) --> "0 $eq 0"
    @fact string(@LinearConstraint(aff-aff ≥ 0)) --> "0 $geq 0"
    # 4-5 AffExpr--QuadExpr
    @fact string(aff2 + q) --> "2.5 y*z + 1.2 y + 7.1 x + 3.7"
    @fact string(aff2 - q) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3"
    @fact_throws  aff2 * q
    @fact_throws  aff2 / q
    @fact string(@QuadConstraint(aff2 ≤ q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $leq 0"
    @fact string(@QuadConstraint(aff2 == q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $eq 0"
    @fact string(@QuadConstraint(aff2 ≥ q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $geq 0"
    # 4-6 AffExpr--SOCExpr
    @fact string(aff + socexpr) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x - w + 0.5"
    @fact string(aff - socexpr) --> "-1.5 $Vert[w,-w + 1]$Vert$sub2 + 7.1 x + w + 4.5"
    @fact_throws aff * socexpr
    @fact_throws aff / socexpr
    @fact_throws @SOCConstraint(aff ≤ socexpr)
    @fact_throws @SOCConstraint(aff == socexpr)
    @fact string(@SOCConstraint(aff ≥ socexpr)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq 7.1 x + w + 4.5"
    end

    # 5. QuadExpr
    context("QuadExpr--???") do
    # 5-0 QuadExpr unary
    @fact string(+q) --> "2.5 y*z + 7.1 x + 2.5"
    @fact string(-q) --> "-2.5 y*z - 7.1 x - 2.5"
    # 5-1 QuadExpr--Number
    @fact string(q + 1.5) --> "2.5 y*z + 7.1 x + 4"
    @fact string(q - 1.5) --> "2.5 y*z + 7.1 x + 1"
    @fact string(q * 2) --> "5 y*z + 14.2 x + 5"
    @fact string(q / 2) --> "1.25 y*z + 3.55 x + 1.25"
    @fact q == q --> true
    @fact string(@QuadConstraint(aff2 ≤ q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $leq 0"
    @fact string(@QuadConstraint(aff2 == q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $eq 0"
    @fact string(@QuadConstraint(aff2 ≥ q)) --> "-2.5 y*z + 1.2 y - 7.1 x - 1.3 $geq 0"
    # 5-2 QuadExpr--Variable
    @fact string(q + w) --> "2.5 y*z + 7.1 x + w + 2.5"
    @fact string(q - w) --> "2.5 y*z + 7.1 x - w + 2.5"
    @fact_throws q*w
    @fact_throws q/w
    @fact string(@QuadConstraint(q ≤ w)) --> "2.5 y*z + 7.1 x - w + 2.5 $leq 0"
    @fact string(@QuadConstraint(q == w)) --> "2.5 y*z + 7.1 x - w + 2.5 $eq 0"
    @fact string(@QuadConstraint(q ≥ w)) --> "2.5 y*z + 7.1 x - w + 2.5 $geq 0"
    # 5-3 QuadExpr--Norm
    @fact_throws q + nrm
    @fact_throws q - nrm
    @fact_throws q * nrm
    @fact_throws q / nrm
    @fact_throws @SOCConstraint(q ≤ nrm)
    @fact_throws @SOCConstraint(q == nrm)
    @fact_throws @SOCConstraint(q ≥ nrm)
    # 5-4 QuadExpr--AffExpr
    @fact string(q + aff2) --> "2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @fact string(q - aff2) --> "2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @fact_throws  q * aff2
    @fact_throws  q / aff2
    @fact string(@QuadConstraint(q ≤ aff2)) --> "2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
    @fact string(@QuadConstraint(q == aff2)) --> "2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
    @fact string(@QuadConstraint(q ≥ aff2)) --> "2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
    # 5-5 QuadExpr--QuadExpr
    @fact string(q + q2) --> "8 x*z + 2.5 y*z + 7.1 x + 1.2 y + 3.7"
    @fact string(q - q2) --> "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3"
    @fact_throws  q * q2
    @fact_throws  q / q2
    @fact string(@QuadConstraint(q ≤ q2)) --> "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $leq 0"
    @fact string(@QuadConstraint(q == q2)) --> "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $eq 0"
    @fact string(@QuadConstraint(q ≥ q2)) --> "-8 x*z + 2.5 y*z + 7.1 x - 1.2 y + 1.3 $geq 0"
    # 4-6 QuadExpr--SOCExpr
    @fact_throws q + socexpr
    @fact_throws q - socexpr
    @fact_throws q * socexpr
    @fact_throws q / socexpr
    @fact_throws @SOCConstraint(q ≤ socexpr)
    @fact_throws @SOCConstraint(q == socexpr)
    @fact_throws @SOCConstraint(q ≥ socexpr)
    end

    # 6. SOCExpr tests
    context("SOCExpr--???") do
    # 6-0 SOCExpr unary
    @fact string(+socexpr) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 2"
    @fact string(-socexpr) --> "-1.5 $Vert[w,-w + 1]$Vert$sub2 + w + 2"
    # 6-1 SOCExpr--Number
    @fact string(socexpr + 1.5) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 0.5"
    @fact string(socexpr - 1.5) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 3.5"
    @fact string(socexpr * 1.5) --> "2.25 $Vert[w,-w + 1]$Vert$sub2 - 1.5 w - 3"
    @fact string(socexpr / 1.5) --> "$Vert[w,-w + 1]$Vert$sub2 - 0.6666666666666666 w - 1.3333333333333333"
    @fact string(@SOCConstraint(socexpr ≤ 1.5)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 3.5"
    @fact_throws @SOCConstraint(socexpr == 1.5)
    @fact_throws @SOCConstraint(socexpr ≥ 1.5)
    # 6-2 SOCExpr--Variable
    @fact string(socexpr + y) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + y - 2"
    @fact string(socexpr - y) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - y - 2"
    @fact_throws socexpr * y
    @fact_throws socexpr / y
    @fact string(@SOCConstraint(socexpr ≤ y)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + y + 2"
    @fact_throws @SOCConstraint(socexpr == y)
    @fact_throws @SOCConstraint(socexpr ≥ y)
    # 6-3 SOCExpr--Norm
    @fact_throws socexpr + nrm
    @fact_throws socexpr - nrm
    @fact_throws socexpr * nrm
    @fact_throws socexpr / nrm
    @fact_throws @SOCConstraint(socexpr ≤ nrm)
    @fact_throws @SOCConstraint(socexpr == nrm)
    @fact_throws @SOCConstraint(socexpr ≥ nrm)
    # 6-4 SOCExpr--AffExpr
    @fact string(socexpr + aff) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w + 7.1 x + 0.5"
    @fact string(socexpr - aff) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 - w - 7.1 x - 4.5"
    @fact_throws socexpr * aff
    @fact_throws socexpr / aff
    @fact string(@SOCConstraint(socexpr ≤ aff)) --> "1.5 $Vert[w,-w + 1]$Vert$sub2 $leq w + 7.1 x + 4.5"
    @fact_throws @SOCConstraint(socexpr == aff)
    @fact_throws @SOCConstraint(socexpr ≥ aff)
    # 6-5 SOCExpr--QuadExpr
    @fact_throws socexpr + q
    @fact_throws socexpr - q
    @fact_throws socexpr * q
    @fact_throws socexpr / q
    @fact_throws @SOCConstraint(socexpr ≤ q)
    @fact_throws @SOCConstraint(socexpr == q)
    @fact_throws @SOCConstraint(socexpr ≥ q)
    # 6-6 SOCExpr--SOCExpr
    @fact_throws socexpr + socexpr
    @fact_throws socexpr - socexpr
    @fact_throws socexpr * socexpr
    @fact_throws socexpr / socexpr
    @fact_throws @SOCConstraint(socexpr ≤ socexpr)
    @fact_throws @SOCConstraint(socexpr == socexpr)
    @fact_throws @SOCConstraint(socexpr ≥ socexpr)
    end
end

facts("[operator] Higher-level operators") do
context("sum") do
    sum_m = Model()
    @variable(sum_m, 0 ≤ matrix[1:3,1:3] ≤ 1, start = 1)
    # sum(j::JuMPArray{Variable})
    @fact string(sum(matrix)) --> "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"
    # sum(j::JuMPArray{Variable}) in a macro
    @objective(sum_m, Max, sum(matrix))
    @fact string(sum_m.obj) --> "matrix[1,1] + matrix[2,1] + matrix[3,1] + matrix[1,2] + matrix[2,2] + matrix[3,2] + matrix[1,3] + matrix[2,3] + matrix[3,3]"

    # sum{T<:Real}(j::JuMPArray{T})
    @fact sum(getvalue(matrix)) --> roughly(9, 1e-6)
    # sum(j::Array{Variable})
    @fact string(sum(matrix[1:3,1:3])) --> string(sum(matrix))
    # sum(affs::Array{AffExpr})
    @fact string(sum([2*matrix[i,j] for i in 1:3, j in 1:3])) --> "2 matrix[1,1] + 2 matrix[2,1] + 2 matrix[3,1] + 2 matrix[1,2] + 2 matrix[2,2] + 2 matrix[3,2] + 2 matrix[1,3] + 2 matrix[2,3] + 2 matrix[3,3]"

    S = [1,3]
    @variable(sum_m, x[S], start=1)
    # sum(j::JuMPDict{Variable})
    @fact length(string(sum(x))) --> 11 # order depends on hashing
    @fact contains(string(sum(x)),"x[1]") --> true
    @fact contains(string(sum(x)),"x[3]") --> true
    # sum{T<:Real}(j::JuMPDict{T})
    @fact sum(getvalue(x)) --> 2
end

context("dot") do
    dot_m = Model()
    @variable(dot_m, 0 ≤ x[1:3] ≤ 1)
    c = vcat(1:3)
    @fact string(dot(c,x)) --> "x[1] + 2 x[2] + 3 x[3]"
    @fact string(dot(x,c)) --> "x[1] + 2 x[2] + 3 x[3]"

    A = [1 3 ; 2 4]
    @variable(dot_m, 1 ≤ y[1:2,1:2] ≤ 1)
    @fact string(vecdot(A,y)) --> "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"
    @fact string(vecdot(y,A)) --> "y[1,1] + 2 y[2,1] + 3 y[1,2] + 4 y[2,2]"

    B = ones(2,2,2)
    @variable(dot_m, 0 ≤ z[1:2,1:2,1:2] ≤ 1)
    @fact string(vecdot(B,z)) --> "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"
    @fact string(vecdot(z,B)) --> "z[1,1,1] + z[2,1,1] + z[1,2,1] + z[2,2,1] + z[1,1,2] + z[2,1,2] + z[1,2,2] + z[2,2,2]"

    @objective(dot_m, Max, dot(x, ones(3)) - vecdot(y, ones(2,2)) )
    #solve(dot_m)
    for i in 1:3
        setvalue(x[i], 1)
    end
    for i in 1:2, j in 1:2
        setvalue(y[i,j], 1)
    end
    @fact dot(c, getvalue(x)) --> roughly( 6, 1e-6)
    @fact vecdot(A, getvalue(y)) --> roughly(10, 1e-6)

    # https://github.com/JuliaOpt/JuMP.jl/issues/656
    issue656 = Model()
    @variable(issue656, x)
    floats = Float64[i for i in 1:2]
    anys   = Array(Any, 2)
    anys[1] = 10
    anys[2] = 20 + x
    @fact dot(floats, anys) --> 10 + 40 + 2x
end
end

module TestHelper # weird scoping behavior with FactCheck...
    using JuMP
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
end

facts("Vectorized operations") do

context("Transpose") do
    m = Model()
    @variable(m, x[1:3])
    @variable(m, y[1:2,1:3])
    @variable(m, z[2:5])
    @fact TestHelper.vec_eq(x', [x[1] x[2] x[3]]) --> true
    @fact TestHelper.vec_eq(transpose(x), [x[1] x[2] x[3]]) --> true
    @fact TestHelper.vec_eq(y', [y[1,1] y[2,1]
                      y[1,2] y[2,2]
                      y[1,3] y[2,3]]) --> true
    @fact TestHelper.vec_eq(transpose(y),
                [y[1,1] y[2,1]
                 y[1,2] y[2,2]
                 y[1,3] y[2,3]]) --> true
    @fact_throws z'
    @fact_throws transpose(z)
end

context("Vectorized arithmetic") do
    m = Model()
    @variable(m, x[1:3])
    A = [2 1 0
         1 2 1
         0 1 2]
    B = sparse(A)
    @fact TestHelper.vec_eq(A*x, [2x[1] +  x[2]
                        2x[2] + x[1] + x[3]
                        x[2] + 2x[3]]) --> true
    @fact TestHelper.vec_eq(A*x, B*x) --> true
    @fact TestHelper.vec_eq(A*x, @JuMP.Expression(B*x)) --> true
    @fact TestHelper.vec_eq(@JuMP.Expression(A*x), @JuMP.Expression(B*x)) --> true
    @fact TestHelper.vec_eq(x'*A, [2x[1]+x[2]; 2x[2]+x[1]+x[3]; x[2]+2x[3]]') --> true
    @fact TestHelper.vec_eq(x'*A, x'*B) --> true
    @fact TestHelper.vec_eq(x'*A, @JuMP.Expression(x'*B)) --> true
    @fact TestHelper.vec_eq(@JuMP.Expression(x'*A), @JuMP.Expression(x'*B)) --> true
    @fact TestHelper.vec_eq(x'*A*x, [2x[1]*x[1] + 2x[1]*x[2] + 2x[2]*x[2] + 2x[2]*x[3] + 2x[3]*x[3]]) --> true
    @fact TestHelper.vec_eq(x'A*x, x'*B*x) --> true
    @fact TestHelper.vec_eq(x'*A*x, @JuMP.Expression(x'*B*x)) --> true
    @fact TestHelper.vec_eq(@JuMP.Expression(x'*A*x), @JuMP.Expression(x'*B*x)) --> true

    y = A*x
    @fact TestHelper.vec_eq(-x, [-x[1], -x[2], -x[3]]) --> true
    @fact TestHelper.vec_eq(-y, [-2x[1] -  x[2]
                       -x[1] - 2x[2] -  x[3]
                               -x[2] - 2x[3]]) --> true
    @fact TestHelper.vec_eq(y + 1, [2x[1] +  x[2]         + 1
                          x[1] + 2x[2] +  x[3] + 1
                                  x[2] + 2x[3] + 1]) --> true
    @fact TestHelper.vec_eq(y - 1, [2x[1] +  x[2]         - 1
                          x[1] + 2x[2] +  x[3] - 1
                                  x[2] + 2x[3] - 1]) --> true
    @fact TestHelper.vec_eq(y + 2ones(3), [2x[1] +  x[2]         + 2
                                 x[1] + 2x[2] +  x[3] + 2
                                         x[2] + 2x[3] + 2]) --> true
    @fact TestHelper.vec_eq(y - 2ones(3), [2x[1] +  x[2]         - 2
                                 x[1] + 2x[2] +  x[3] - 2
                                         x[2] + 2x[3] - 2]) --> true
    @fact TestHelper.vec_eq(2ones(3) + y, [2x[1] +  x[2]         + 2
                                 x[1] + 2x[2] +  x[3] + 2
                                         x[2] + 2x[3] + 2]) --> true
    @fact TestHelper.vec_eq(2ones(3) - y, [-2x[1] -  x[2]         + 2
                                 -x[1] - 2x[2] -  x[3] + 2
                                         -x[2] - 2x[3] + 2]) --> true
    @fact TestHelper.vec_eq(y + x, [3x[1] +  x[2]
                          x[1] + 3x[2] +  x[3]
                                  x[2] + 3x[3]]) --> true
    @fact TestHelper.vec_eq(x + y, [3x[1] +  x[2]
                          x[1] + 3x[2] +  x[3]
                                  x[2] + 3x[3]]) --> true
    @fact TestHelper.vec_eq(2y + 2x, [6x[1] + 2x[2]
                           2x[1] + 6x[2] + 2x[3]
                                   2x[2] + 6x[3]]) --> true
    @fact TestHelper.vec_eq(y - x, [ x[1] + x[2]
                          x[1] + x[2] + x[3]
                                 x[2] + x[3]]) --> true
    @fact TestHelper.vec_eq(x - y, [-x[1] - x[2]
                         -x[1] - x[2] - x[3]
                                -x[2] - x[3]]) --> true
    @fact TestHelper.vec_eq(y + x[:], [3x[1] +  x[2]
                             x[1] + 3x[2] +  x[3]
                                     x[2] + 3x[3]]) --> true
    @fact TestHelper.vec_eq(x[:] + y, [3x[1] +  x[2]
                             x[1] + 3x[2] +  x[3]
                                     x[2] + 3x[3]]) --> true

    @fact TestHelper.vec_eq(@JuMP.Expression(A*x/2), A*x/2) --> true
end

context("Dot-ops") do
    m = Model()
    @variable(m, x[1:2,1:2])
    A = [1 2;
         3 4]
    B = sparse(A)
    y = SparseMatrixCSC(2, 2, copy(B.colptr), copy(B.rowval), vec(x))
    @fact TestHelper.vec_eq(A.+x, [1+x[1,1]  2+x[1,2];
                                   3+x[2,1]  4+x[2,2]]) --> true
    @fact TestHelper.vec_eq(A.+x, B.+x) --> true
    # @fact TestHelper.vec_eq(A.+x, A.+y) --> true
    # @fact TestHelper.vec_eq(A.+y, B.+y) --> true
    @fact TestHelper.vec_eq(x.+A, [1+x[1,1]  2+x[1,2];
                                   3+x[2,1]  4+x[2,2]]) --> true
    @fact TestHelper.vec_eq(x.+A, x.+B) --> true
    @fact TestHelper.vec_eq(x.+A, y.+A) --> true
    @fact TestHelper.vec_eq(x .+ x, [2x[1,1] 2x[1,2]; 2x[2,1] 2x[2,2]]) --> true
    # @fact TestHelper.vec_eq(y.+A, y.+B) --> true
    @fact TestHelper.vec_eq(A.-x, [1-x[1,1]  2-x[1,2];
                                   3-x[2,1]  4-x[2,2]]) --> true
    @fact TestHelper.vec_eq(A.-x, B.-x) --> true
    @fact TestHelper.vec_eq(A.-x, A.-y) --> true
    @fact TestHelper.vec_eq(x .- x, [zero(AffExpr) for _1 in 1:2, _2 in 1:2]) --> true
    # @fact TestHelper.vec_eq(A.-y, B.-y) --> true
    @fact TestHelper.vec_eq(x.-A, [-1+x[1,1]  -2+x[1,2];
                                   -3+x[2,1]  -4+x[2,2]]) --> true
    @fact TestHelper.vec_eq(x.-A, x.-B) --> true
    @fact TestHelper.vec_eq(x.-A, y.-A) --> true
    # @fact TestHelper.vec_eq(y.-A, y.-B) --> true
    @fact TestHelper.vec_eq(A.*x, [1*x[1,1]  2*x[1,2];
                                   3*x[2,1]  4*x[2,2]]) --> true
    @fact TestHelper.vec_eq(A.*x, B.*x) --> true
    @fact TestHelper.vec_eq(A.*x, A.*y) --> true
    # @fact TestHelper.vec_eq(A.*y, B.*y) --> true

    @fact TestHelper.vec_eq(x.*A, [1*x[1,1]  2*x[1,2];
                                   3*x[2,1]  4*x[2,2]]) --> true
    @fact TestHelper.vec_eq(x.*A, x.*B) --> true
    @fact TestHelper.vec_eq(x.*A, y.*A) --> true
    # @fact TestHelper.vec_eq(y.*A, y.*B) --> true

    @fact TestHelper.vec_eq(x .* x, [x[1,1]^2 x[1,2]^2; x[2,1]^2 x[2,2]^2]) --> true
    @fact_throws TestHelper.vec_eq(A./x, [1*x[1,1]  2*x[1,2];
                                          3*x[2,1]  4*x[2,2]])
    @fact TestHelper.vec_eq(x./A, [1/1*x[1,1]  1/2*x[1,2];
                                   1/3*x[2,1]  1/4*x[2,2]]) --> true
    @fact TestHelper.vec_eq(x./A, x./B) --> true
    @fact TestHelper.vec_eq(x./A, y./A) --> true
    # @fact TestHelper.vec_eq(A./y, B./y) --> true

    @fact TestHelper.vec_eq((2*x) / 3, full((2*y) / 3)) --> true
    @fact TestHelper.vec_eq(2 * (x/3), full(2 * (y/3))) --> true
    @fact TestHelper.vec_eq(x[1,1] * A, full(x[1,1] * B)) --> true
end

context("Vectorized comparisons") do
    m = Model()
    @variable(m, x[1:3])
    A = [1 2 3
         0 4 5
         6 0 7]
    B = sparse(A)
    @constraint(m, x'*A*x .>= 1)
    @fact TestHelper.vec_eq(m.quadconstr[1].terms, [x[1]*x[1] + 2x[1]*x[2] + 4x[2]*x[2] + 9x[1]*x[3] + 5x[2]*x[3] + 7x[3]*x[3] - 1]) --> true
    @fact m.quadconstr[1].sense --> :(>=)
    @constraint(m, x'*A*x .>= 1)
    @fact TestHelper.vec_eq(m.quadconstr[1].terms, m.quadconstr[2].terms) --> true

    mat = [ 3x[1] + 12x[3] +  4x[2]
            2x[1] + 12x[2] + 10x[3]
           15x[1] +  5x[2] + 21x[3]]

    @constraint(m, (x'A)' + 2A*x .<= 1)
    terms = map(v->v.terms, m.linconstr[1:3])
    lbs   = map(v->v.lb,    m.linconstr[1:3])
    ubs   = map(v->v.ub,    m.linconstr[1:3])
    @fact TestHelper.vec_eq(terms, mat) --> true
    @fact lbs --> fill(-Inf, 3)
    @fact ubs --> fill(   1, 3)
    @fact TestHelper.vec_eq((x'A)' + 2A*x, (x'A)' + 2B*x) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, (x'B)' + 2A*x) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, (x'B)' + 2B*x) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'A)' + 2A*x)) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'B)' + 2A*x)) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'A)' + 2B*x)) --> true
    @fact TestHelper.vec_eq((x'A)' + 2A*x, @JuMP.Expression((x'B)' + 2B*x)) --> true

    @constraint(m, -1 .<= (x'A)' + 2A*x .<= 1)
    terms = map(v->v.terms, m.linconstr[4:6])
    lbs   = map(v->v.lb,    m.linconstr[4:6])
    ubs   = map(v->v.ub,    m.linconstr[4:6])
    @fact TestHelper.vec_eq(terms, mat) --> true
    @fact lbs --> fill(-1, 3)
    @fact ubs --> fill( 1, 3)

    @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= 1)
    terms = map(v->v.terms, m.linconstr[7:9])
    lbs   = map(v->v.lb,    m.linconstr[7:9])
    ubs   = map(v->v.ub,    m.linconstr[7:9])
    @fact TestHelper.vec_eq(terms, mat) --> true
    @fact lbs --> -[1:3;]
    @fact ubs --> fill( 1, 3)

    @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= [3:-1:1;])
    terms = map(v->v.terms, m.linconstr[10:12])
    lbs   = map(v->v.lb,    m.linconstr[10:12])
    ubs   = map(v->v.ub,    m.linconstr[10:12])
    @fact TestHelper.vec_eq(terms, mat) --> true
    @fact lbs --> -[1:3;]
    @fact ubs --> [3:-1:1;]

    @constraint(m, -[1:3;] .<= (x'A)' + 2A*x .<= 3)
    terms = map(v->v.terms, m.linconstr[13:15])
    lbs   = map(v->v.lb,    m.linconstr[13:15])
    ubs   = map(v->v.ub,    m.linconstr[13:15])
    @fact TestHelper.vec_eq(terms, mat) --> true
    @fact lbs --> -[1:3;]
    @fact ubs --> fill(3,3)
end

end

facts("[operator] JuMPArray concatenation") do
    m = Model()
    @variable(m, x[1:3])
    @variable(m, y[1:3,1:3])
    @variable(m, z[1:1])
    @variable(m, w[1:1,1:3])

    @fact TestHelper.vec_eq([x y], [x[1] y[1,1] y[1,2] y[1,3]
                                    x[2] y[2,1] y[2,2] y[2,3]
                                    x[3] y[3,1] y[3,2] y[3,3]]) --> true

    @fact TestHelper.vec_eq([x 2y+1], [x[1] 2y[1,1]+1 2y[1,2]+1 2y[1,3]+1
                                       x[2] 2y[2,1]+1 2y[2,2]+1 2y[2,3]+1
                                       x[3] 2y[3,1]+1 2y[3,2]+1 2y[3,3]+1]) --> true

    @fact TestHelper.vec_eq([1 x'], [1 x[1] x[2] x[3]]) --> true
    @fact TestHelper.vec_eq([2x;x], [2x[1],2x[2],2x[3],x[1],x[2],x[3]]) --> true
    # vcat on JuMPArray
    @fact TestHelper.vec_eq([x;x], [x[1],x[2],x[3],x[1],x[2],x[3]]) --> true
    # hcat on JuMPArray
    @fact TestHelper.vec_eq([x x], [x[1] x[1]
                                    x[2] x[2]
                                    x[3] x[3]]) --> true
    # hvcat on JuMPArray
    tmp1 = [z w; x y]
    tmp2 = [z[1] w[1,1] w[1,2] w[1,3]
            x[1] y[1,1] y[1,2] y[1,3]
            x[2] y[2,1] y[2,2] y[2,3]
            x[3] y[3,1] y[3,2] y[3,3]]
    @fact TestHelper.vec_eq(tmp1, tmp2) --> true
    tmp3 = [1 2x'
            x 2y-x*x']
    tmp4 = [1    2x[1]               2x[2]               2x[3]
            x[1] -x[1]*x[1]+2y[1,1]  -x[1]*x[2]+2y[1,2]  -x[1]*x[3] + 2y[1,3]
            x[2] -x[1]*x[2]+2y[2,1]  -x[2]*x[2]+2y[2,2]  -x[2]*x[3] + 2y[2,3]
            x[3] -x[1]*x[3]+2y[3,1]  -x[2]*x[3]+2y[3,2]  -x[3]*x[3] + 2y[3,3]]
    @fact TestHelper.vec_eq(tmp3, tmp4) --> true

    A = sprand(3, 3, 0.2)
    B = full(A)
    @fact TestHelper.vec_eq([A y], [B y]) --> true
end


# The behavior in this test is no longer well-defined
#let
#    dot_m = Model()
#    @variable(dot_m, 0 <= v[1:3,1:4] <= 1)
#    @variable(dot_m, 0 <= w[3:2:7,7:-2:1] <= 1)
#    @fact string(dot(v,w)) --> "v[1,1]*w[3,7] + v[1,2]*w[3,5] + v[1,3]*w[3,3] + v[1,4]*w[3,1] + v[2,1]*w[5,7] + v[2,2]*w[5,5] + v[2,3]*w[5,3] + v[2,4]*w[5,1] + v[3,1]*w[7,7] + v[3,2]*w[7,5] + v[3,3]*w[7,3] + v[3,4]*w[7,1]"
#
#    @constraint(dot_m, sum(v) + sum(w) --> 2)
#    @objective(dot_m, :Max, v[2,3] + w[5,3])
#    solve(dot_m)
#    @fact dot(getvalue(v),getvalue(w)) --> 1.0
#end

# Ditto
#let
#    slice_m = Model()
#    C = [:cat,:dog]
#    @variable(slice_m, x[-1:1,C])
#    catcoef = [1,2,3]
#    dogcoef = [3,4,5]
#
#    @fact string(dot(catcoef, x[:,:cat])) --> "x[-1,cat] + 2 x[0,cat] + 3 x[1,cat]"
#    @fact string(dot(dogcoef, x[:,:dog])) --> "3 x[-1,dog] + 4 x[0,dog] + 5 x[1,dog]"
#end
