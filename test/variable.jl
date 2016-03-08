#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/variable.jl
# Testing for Variable
#############################################################################
using JuMP, FactCheck

facts("[variable] constructors") do
    # Constructors
    mcon = Model()
    @defVar(mcon, nobounds)
    @defVar(mcon, lbonly >= 0)
    @defVar(mcon, ubonly <= 1)
    @defVar(mcon, 0 <= bothb <= 1)
    @defVar(mcon, 0 <= onerange[-5:5] <= 10)
    @defVar(mcon, onerangeub[-7:1] <= 10, Int)
    @defVar(mcon, manyrangelb[0:1,10:20,1:1] >= 2)
    @fact getLower(manyrangelb[0,15,1]) --> 2
    s = ["Green","Blue"]
    @defVar(mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
    @fact getUpper(x[-4,"Green"]) --> 5.5
    @fact getValue(x[-3,"Blue"]) --> -2
    @fact isequal(getVar(mcon, :lbonly),lbonly) --> true
    @fact isequal(getVar(mcon, :ubonly),ubonly) --> true
    @fact isequal(getVar(mcon, :onerangeub)[-7],onerangeub[-7]) --> true
    @defVar(mcon, lbonly)
    @fact_throws ErrorException getVar(mcon, :lbonly)
    @fact_throws ErrorException getVar(mcon, :foo)
    d = Dict()
    @defVar(mcon, d["bar"][1:10] == 1)
    @fact getValue(d["bar"][1]) --> 1
    @fact typeof(zero(nobounds)) --> AffExpr
    @fact typeof(one(nobounds)) --> AffExpr
end

facts("[variable] get and set bounds") do
    m = Model()
    @defVar(m, 0 <= x <= 2)
    @fact getLower(x) --> 0
    @fact getUpper(x) --> 2
    setLower(x, 1)
    @fact getLower(x) --> 1
    setUpper(x, 3)
    @fact getUpper(x) --> 3
    @defVar(m, y, Bin)
    @fact getLower(y) --> 0
    @fact getUpper(y) --> 1
    @defVar(m, 0 <= y <= 1, Bin)
    @fact getLower(y) --> 0
    @fact getUpper(y) --> 1
    @defVar(m, fixedvar == 2)
    @fact getValue(fixedvar) --> 2
    @fact getLower(fixedvar) --> 2
    @fact getUpper(fixedvar) --> 2
    setValue(fixedvar, 5)
    @fact getValue(fixedvar) --> 5
    @fact getLower(fixedvar) --> 5
    @fact getUpper(fixedvar) --> 5
end

facts("[variable] get and set values") do
    m = Model()
    @defVar(m, x[1:3])
    x0 = collect(1:3)
    setValue(x, x0)
    @fact getValue(x) --> x0
    @fact getValue([x[1],x[2],x[3]]) --> x0

    @defVar(m, y[1:3,1:2])
    @fact_throws DimensionMismatch setValue(y, collect(1:6))
end

facts("[variable] get and set category") do
    m = Model()
    @defVar(m, x[1:3])
    setCategory(x[2], :Int)
    @fact getCategory(x[3]) --> :Cont
    @fact getCategory(x[2]) --> :Int
end

facts("[variable] repeated elements in index set (issue #199)") do
    repeatmod = Model()
    s = [:x,:x,:y]
    @defVar(repeatmod, x[s])
    @fact MathProgBase.numvar(repeatmod) --> 3
end

# Test conditions in variable definition
if VERSION >= v"0.4-"
    fa = repl[:for_all]
    inset, dots = repl[:in], repl[:dots]

    facts("[variable] condition in indexing") do
        condmod = Model()
        @defVar(condmod, x[i=1:10; iseven(i)])
        @defVar(condmod, y[j=1:10,k=3:2:9; isodd(j+k) && k <= 8])
        @fact length(x.tupledict) --> 5
        @fact length(y.tupledict) --> 15
        @fact string(condmod) --> "Min 0\nSubject to\n x[i] free $fa i $inset {1,2,$dots,9,10} s.t. iseven(i)\n y[j,k] free $fa j $inset {1,2,$dots,9,10}, k $inset {3,5,7,9} s.t. isodd(j + k) and k <= 8\n"
    end
end

facts("[variable] @defVar returning Array{Variable}") do
    m = Model()
    @defVar(m, x[1:3,1:4,1:2])
    @defVar(m, y[1:0])
    @defVar(m, z[1:4])

    @fact typeof(x) --> Array{Variable,3}
    @fact typeof(y) --> Array{Variable,1}
    @fact typeof(z) --> Array{Variable,1}

    @fact typeof(getValue(x)) --> Array{Float64,3}
    @fact typeof(getValue(y)) --> Array{Float64,1}
    @fact typeof(getValue(z)) --> Array{Float64,1}
end

facts("[variable] getValue on empty things") do
    m = Model()
    @defVar(m, x[1:4,  1:0,1:3])   # Array{Variable}
    @defVar(m, y[1:4,  2:1,1:3]) # JuMPArray
    @defVar(m, z[1:4,Set(),1:3]) # JuMPDict

    @fact getValue(x) --> Array(Float64, 4, 0, 3)
    @fact typeof(getValue(y)) <: JuMP.JuMPArray{Float64} --> true
    @fact size(getValue(y)) --> (4,0,3)
    @fact typeof(getValue(z)) --> JuMP.JuMPArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
    @fact length(getValue(z)) --> 0
end

# Slices three-dimensional JuMPContainer x[I+oi,J+oj,K+ok]
# I,J,K can be singletons, ranges, colons, etc.
# oi, oj, ok are scalar offsets. The last can be omitted, as
# that's useful to test more exotic indexing ([:red,"blue"])
function sliceof(x, I, J, K, oi, oj, ok=nothing)
    if length(K) == 1
        if length(J) == 1
            y = Array(Variable, length(I))
        else
            y = Array(Variable, length(I), length(J))
        end
    else
        y = Array(Variable, length(I), length(J), length(K))
    end
    ii = 1
    jj = 1
    kk = 1
    for i in I
        for j in J
            for k in K
                if ok == nothing
                    y[ii,jj,kk] = x[i+oi,j+oj,k]
                else
                    y[ii,jj,kk] = x[i+oi,j+oj,k+ok]
                end
                kk += 1
            end
            jj += 1
            kk = 1
        end
        ii += 1
        jj = 1
    end
    y
end

facts("[variable] Slices of JuMPArray (#684)") do
    m = Model()
    @defVar(m, x[1:3, 1:4,1:2])
    @defVar(m, y[1:3,-1:2,3:4])
    @defVar(m, z[1:3,-1:2:4,3:4])
    @defVar(m, w[1:3,-1:2,[:red,"blue"]])

    @fact x[:] --> vec(sliceof(x, 1:3, 1:4, 1:2, 0, 0, 0))
    @fact x[:,:,:] --> sliceof(x, 1:3, 1:4, 1:2, 0, 0, 0)
    @fact x[1,:,:] --> sliceof(x, 1, 1:4, 1:2, 0, 0, 0)
    @fact x[1,:,2] --> sliceof(x, 1, 1:4, 2, 0, 0, 0)
    @fact_throws x[1,:,3]
    @fact x[1:2,:,:] --> sliceof(x, 1:2, 1:4, 1:2, 0, 0, 0)
    @fact x[1:2,:,2] --> sliceof(x, 1:2, 1:4, 2, 0, 0, 0)
    @fact x[1:2,:,1:2] --> sliceof(x, 1:2, 1:4, 1:2, 0, 0, 0)
    @fact_throws x[1:2,:,1:3]

    @fact y[:] --> vec(sliceof(y, 1:3, 1:4, 1:2, 0, -2, 2))
    @fact y[:,:,:] --> sliceof(y, 1:3, 1:4, 1:2, 0, -2, 2)
    @fact y[1,:,:] --> sliceof(y, 1, 1:4, 1:2, 0, -2, 2)
    @fact y[1,:,4] --> sliceof(y, 1, 1:4, 2, 0, -2, 2)
    @fact_throws y[1,:,5]
    @fact y[1:2,:,:] --> sliceof(y, 1:2, 1:4, 1:2, 0, -2, 2)
    @fact y[1:2,:,4] --> sliceof(y, 1:2, 1:4, 2, 0, -2, 2)
    @fact y[1:2,:,3:4] --> sliceof(y, 1:2, 1:4, 1:2, 0, -2, 2)
    @fact_throws y[1:2,:,1:3]

    @fact z[:] --> vec(sliceof(z, 1:3, 1:2:6, 1:2, 0, -2, 2))
    @fact z[:,1,:] --> sliceof(z, 1:3, 3, 1:2, 0, -2, 2)
    @fact z[1,1,:] --> sliceof(z, 1, 3, 1:2, 0, -2, 2)
    @fact_throws z[:,5,3]
    @fact z[1:2,1,:] --> sliceof(z, 1:2, 3, 1:2, 0, -2, 2)
    @fact z[1:2,1,4] --> sliceof(z, 1:2, 3, 2, 0, -2, 2)
    @fact z[1:2,1,3:4] --> sliceof(z, 1:2, 3, 1:2, 0, -2, 2)
    @fact_throws z[1:2,1,1:3]

    @fact w[:] --> vec(sliceof(w, 1:3, 1:4, [:red,"blue"], 0, -2))
    @fact_throws w[:,:,:]
    @fact w[1,:,"blue"] --> sliceof(w, 1, 1:4, ["blue"], 0, -2)
    @fact w[1,:,:red] --> sliceof(w, 1, 1:4, [:red], 0, -2)
    @fact_throws w[1,:,"green"]
    @fact w[1:2,:,"blue"] --> sliceof(w, 1:2, 1:4, ["blue"], 0, -2)
    @fact_throws w[1:2,:,[:red,"blue"]]
end
