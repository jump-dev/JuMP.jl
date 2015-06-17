#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
    @defVar(maff, 0 <= x[1:5] <= 1)
    @defVar(maff, 0 <= LongName <= 99)

    context("AffExpr") do
        # Test affToStr
        a1 = x[1] + LongName + 5
        @fact affToStr(a1) => "x[1] + LongName + 5"
        # Test like term collection
        a2 = 2*(x[2] + LongName + x[2]) + 0
        @fact affToStr(a2) => "4 x[2] + 2 LongName"
        # Test appending functionality
        push!(a1, 5.0, x[2])
        @fact affToStr(a1) => "x[1] + LongName + 5 x[2] + 5"
        append!(a1, a2)
        @fact affToStr(a1) => "x[1] + 3 LongName + 9 x[2] + 5"
    end

    context("QuadExpr") do
        # Test quadToStr
        q1 = x[1]*x[2] + 27.2*LongName + 5
        @fact quadToStr(q1) => "x[1]*x[2] + 27.2 LongName + 5"
        # Test like term collection
        q2 = x[1]*x[2] + x[2]*x[1]
        @fact quadToStr(q2) => "2 x[1]*x[2]"
    end
end

facts("[expr] Test getValue(expr)") do
    m = Model()
    @defVar(m, 1 <= x[1:3] <= 2)
    setValue(x[3], 2)
    setValue(x[2], 2)
    setValue(x[1], 1)
    @fact getValue(x[1]-x[2]+2x[3]-1.0) => roughly(2.0)
    @fact getValue(x[1]*x[1]-2x[2]*x[1]+3x[2]+1) => roughly(4.0)
end

facts("[expr] Test expression iterators") do
    m = Model()
    @defVar(m, x[1:10])

    a1 = 1*x[1] + 2*x[2]
    k = 1
    for (coeff,var) in a1
        if k == 1
            @fact coeff => 1
            @fact var => exactly(x[1])
        elseif k == 2
            @fact coeff => 2
            @fact var => exactly(x[2])
        end
        k += 1
    end

    a2 = zero(AffExpr)
    for (coeff, var) in a2
        @fact coeff => 0.0  # Shouldn't be called!
    end
end

# Test ``in(::Variable, AffExpr)``
facts("[expr] Test in(::Variable, ::AffExpr)") do
    m = Model()
    @defVar(m, x[1:3])
    @defVar(m, y)
    @fact x[2] in 2x[2] + x[1] => true
    @fact x[3] in x[1] + 2x[2] => false
    @fact y in @defExpr(sum{i*x[i],i=1:3}) => false
    @fact x[2] in x[1] + 2x[2] - x[2] + x[3] - x[2] => false
end
