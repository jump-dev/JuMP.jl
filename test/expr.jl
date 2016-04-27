#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/expr.jl
# Testing for AffExpr and QuadExpr
#############################################################################
using JuMP, FactCheck

facts("[expr] Test expression construction") do
    maff = Model()
    @variable(maff, 0 <= x[1:5] <= 1)
    @variable(maff, 0 <= LongName <= 99)

    context("AffExpr") do
        # Test affstr
        a1 = x[1] + LongName + 5
        @fact string(a1) --> "x[1] + LongName + 5"
        # Test like term collection
        a2 = 2*(x[2] + LongName + x[2]) + 0
        @fact string(a2) --> "4 x[2] + 2 LongName"
        # Test appending functionality
        push!(a1, 5.0, x[2])
        @fact string(a1) --> "x[1] + LongName + 5 x[2] + 5"
        append!(a1, a2)
        @fact string(a1) --> "x[1] + 3 LongName + 9 x[2] + 5"
        append!(a1, 2.0)
        @fact string(a1) --> "x[1] + 3 LongName + 9 x[2] + 7"
        append!(a1, LongName)
        @fact string(a1) --> "x[1] + 4 LongName + 9 x[2] + 7"
    end

    context("QuadExpr") do
        # Test string
        q1 = x[1]*x[2] + 27.2*LongName + 5
        @fact string(q1) --> "x[1]*x[2] + 27.2 LongName + 5"
        # Test like term collection
        q2 = x[1]*x[2] + x[2]*x[1]
        @fact string(q2) --> "2 x[1]*x[2]"
    end
end

facts("[expr] Test getvalue(expr)") do
    m = Model()
    @variable(m, 1 <= x[1:3] <= 2)
    setvalue(x[3], 2)
    setvalue(x[2], 2)
    setvalue(x[1], 1)
    @fact getvalue(x[1]-x[2]+2x[3]-1.0) --> roughly(2.0)
    @fact getvalue(x[1]*x[1]-2x[2]*x[1]+3x[2]+1) --> roughly(4.0)
end

facts("[expr] Test expression iterators") do
    m = Model()
    @variable(m, x[1:10])

    a1 = 1*x[1] + 2*x[2]
    k = 1
    for (coeff,var) in linearterms(a1)
        if k == 1
            @fact coeff --> 1
            @fact var --> exactly(x[1])
        elseif k == 2
            @fact coeff --> 2
            @fact var --> exactly(x[2])
        end
        k += 1
    end

    a2 = zero(AffExpr)
    for (coeff, var) in linearterms(a2)
        @fact coeff --> 0.0  # Shouldn't be called!
    end
end
