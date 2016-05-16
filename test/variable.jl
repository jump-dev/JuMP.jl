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
    @variable(mcon, nobounds)
    @variable(mcon, lbonly >= 0)
    @variable(mcon, ubonly <= 1)
    @variable(mcon, 0 <= bothb <= 1)
    @variable(mcon, 0 <= onerange[-5:5] <= 10)
    @variable(mcon, onerangeub[-7:1] <= 10, Int)
    @variable(mcon, manyrangelb[0:1,10:20,1:1] >= 2)
    @fact getlowerbound(manyrangelb[0,15,1]) --> 2
    s = ["Green","Blue"]
    @variable(mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
    @fact getupperbound(x[-4,"Green"]) --> 5.5
    @fact getvalue(x[-3,"Blue"]) --> -2
    @fact isequal(getvariable(mcon, :lbonly),lbonly) --> true
    @fact isequal(getvariable(mcon, :ubonly),ubonly) --> true
    @fact isequal(getvariable(mcon, :onerangeub)[-7],onerangeub[-7]) --> true
    @variable(mcon, lbonly)
    @fact_throws ErrorException getvariable(mcon, :lbonly)
    @fact_throws ErrorException getvariable(mcon, :foo)
    d = Dict()
    @variable(mcon, d["bar"][1:10] == 1)
    @fact getvalue(d["bar"][1]) --> 1
    @fact typeof(zero(nobounds)) --> AffExpr
    @fact typeof(one(nobounds)) --> AffExpr
end

facts("[variable] get and set bounds") do
    m = Model()
    @variable(m, 0 <= x <= 2)
    @fact getlowerbound(x) --> 0
    @fact getupperbound(x) --> 2
    setlowerbound(x, 1)
    @fact getlowerbound(x) --> 1
    setupperbound(x, 3)
    @fact getupperbound(x) --> 3
    @variable(m, y, Bin)
    @fact getlowerbound(y) --> 0
    @fact getupperbound(y) --> 1
    @variable(m, 0 <= y <= 1, Bin)
    @fact getlowerbound(y) --> 0
    @fact getupperbound(y) --> 1
    @variable(m, fixedvar == 2)
    @fact getvalue(fixedvar) --> 2
    @fact getlowerbound(fixedvar) --> 2
    @fact getupperbound(fixedvar) --> 2
    setvalue(fixedvar, 5)
    @fact getvalue(fixedvar) --> 5
    @fact getlowerbound(fixedvar) --> 5
    @fact getupperbound(fixedvar) --> 5
end

facts("[variable] get and set values") do
    m = Model()
    @variable(m, x[1:3])
    x0 = collect(1:3)
    setvalue(x, x0)
    @fact getvalue(x) --> x0
    @fact getvalue([x[1],x[2],x[3]]) --> x0

    @variable(m, y[1:3,1:2])
    @fact_throws DimensionMismatch setvalue(y, collect(1:6))
end

facts("[variable] get and set category") do
    m = Model()
    @variable(m, x[1:3])
    setcategory(x[2], :Int)
    @fact getcategory(x[3]) --> :Cont
    @fact getcategory(x[2]) --> :Int
end

facts("[variable] repeated elements in index set (issue #199)") do
    repeatmod = Model()
    s = [:x,:x,:y]
    @variable(repeatmod, x[s])
    @fact MathProgBase.numvar(repeatmod) --> 3
end

facts("[variable] condition in indexing") do
    fa = repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    condmod = Model()
    @variable(condmod, x[i=1:10; iseven(i)])
    @variable(condmod, y[j=1:10,k=3:2:9; isodd(j+k) && k <= 8])
    @fact length(x.tupledict) --> 5
    @fact length(y.tupledict) --> 15
    @fact string(condmod) --> "Min 0\nSubject to\n x[i] free $fa i $inset {1,2,$dots,9,10} s.t. iseven(i)\n y[j,k] free $fa j $inset {1,2,$dots,9,10}, k $inset {3,5,7,9} s.t. isodd(j + k) and k <= 8\n"
end

facts("[variable] @variable returning Array{Variable}") do
    m = Model()
    @variable(m, x[1:3,1:4,1:2])
    @variable(m, y[1:0])
    @variable(m, z[1:4])

    @fact typeof(x) --> Array{Variable,3}
    @fact typeof(y) --> Array{Variable,1}
    @fact typeof(z) --> Array{Variable,1}

    @fact typeof(getvalue(x)) --> Array{Float64,3}
    @fact typeof(getvalue(y)) --> Array{Float64,1}
    @fact typeof(getvalue(z)) --> Array{Float64,1}
end

facts("[variable] getvalue on empty things") do
    m = Model()
    @variable(m, x[1:4,  1:0,1:3])   # Array{Variable}
    @variable(m, y[1:4,  2:1,1:3]) # JuMPArray
    @variable(m, z[1:4,Set(),1:3]) # JuMPDict

    @fact getvalue(x) --> Array(Float64, 4, 0, 3)
    @fact typeof(getvalue(y)) <: JuMP.JuMPArray{Float64} --> true
    @fact JuMP.size(getvalue(y)) --> (4,0,3)
    @fact typeof(getvalue(z)) --> JuMP.JuMPArray{Float64,3,Tuple{UnitRange{Int},Set{Any},UnitRange{Int}}}
    @fact length(getvalue(z)) --> 0
end

# Slices three-dimensional JuMPContainer x[I,J,K]
# I,J,K can be singletons, ranges, colons, etc.
function sliceof(x, I, J, K)
    y = Array(Variable, length(I), length(J), length(K))

    ii = 1
    jj = 1
    kk = 1
    for i in I
        for j in J
            for k in K
                y[ii,jj,kk] = x[i,j,k]
                kk += 1
            end
            jj += 1
            kk = 1
        end
        ii += 1
        jj = 1
    end
    if VERSION >= v"0.5-"
        idx = [length(I)==1, length(J)==1, length(K)==1]
    else
        idx = [false, false, length(K)==1]
        if idx[3] == true && length(J) == 1
            idx[2] = true
        end
    end
    squeeze(y, tuple(find(idx)...))
end

facts("[variable] Slices of JuMPArray (#684)") do
    m = Model()
    @variable(m, x[1:3, 1:4,1:2])
    @variable(m, y[1:3,-1:2,3:4])
    @variable(m, z[1:3,-1:2:4,3:4])
    @variable(m, w[1:3,-1:2,[:red,"blue"]])

    @fact x[:] --> vec(sliceof(x, 1:3, 1:4, 1:2))
    @fact x[:,:,:] --> sliceof(x, 1:3, 1:4, 1:2)
    @fact x[1,:,:] --> sliceof(x, 1, 1:4, 1:2)
    @fact x[1,:,2] --> sliceof(x, 1, 1:4, 2)
    @fact_throws x[1,:,3]
    @fact x[1:2,:,:] --> sliceof(x, 1:2, 1:4, 1:2)
    @fact x[1:2,:,2] --> sliceof(x, 1:2, 1:4, 2)
    @fact x[1:2,:,1:2] --> sliceof(x, 1:2, 1:4, 1:2)
    @fact_throws x[1:2,:,1:3]

    @fact y[:] --> vec(sliceof(y, 1:3, -1:2, 3:4))
    @fact y[:,:,:] --> sliceof(y, 1:3, -1:2, 3:4)
    @fact y[1,:,:] --> sliceof(y, 1, -1:2, 3:4)
    @fact y[1,:,4] --> sliceof(y, 1, -1:2, 4)
    @fact_throws y[1,:,5]
    @fact y[1:2,:,:] --> sliceof(y, 1:2, -1:2, 3:4)
    @fact y[1:2,:,4] --> sliceof(y, 1:2, -1:2, 4)
    @fact y[1:2,:,3:4] --> sliceof(y, 1:2, -1:2, 3:4)
    @fact_throws y[1:2,:,1:3]

    @fact z[:] --> vec(sliceof(z, 1:3, -1:2:4, 3:4))
    @fact z[:,1,:] --> sliceof(z, 1:3, 1, 3:4)
    @fact z[1,1,:] --> sliceof(z, 1, 1, 3:4)
    @fact_throws z[:,5,3]
    @fact z[1:2,1,:] --> sliceof(z, 1:2, 1, 3:4)
    @fact z[1:2,1,4] --> sliceof(z, 1:2, 1, 4)
    @fact z[1:2,1,3:4] --> sliceof(z, 1:2, 1, 3:4)
    @fact_throws z[1:2,1,1:3]

    @fact w[:] --> vec(sliceof(w, 1:3, -1:2, [:red,"blue"]))
    @fact_throws w[:,:,:]
    @fact w[1,:,"blue"] --> sliceof(w, 1, -1:2, ["blue"])
    @fact w[1,:,:red] --> sliceof(w, 1, -1:2, [:red])
    @fact_throws w[1,:,"green"]
    @fact w[1:2,:,"blue"] --> sliceof(w, 1:2, -1:2, ["blue"])
    @fact_throws w[1:2,:,[:red,"blue"]]
end

facts("[variable] Can't use end for indexing a JuMPContainer") do
    m = Model()
    @variable(m, x[0:2,1:4])
    @variable(m, y[i=1:4,j=1:4;true])
    @variable(m, z[0:2])
    @fact_throws x[end,1]
    @fact_throws x[end-1]
    @fact_throws x[0,end-1]
    @fact_throws y[end,end-1]
    @fact_throws y[end,1]
    @fact_throws z[end]
end
