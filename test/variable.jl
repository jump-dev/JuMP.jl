#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
    @fact getLower(manyrangelb[0,15,1]) => 2
    s = ["Green","Blue"]
    @defVar(mcon, x[i=-10:10,s] <= 5.5, Int, start=i+1)
    @fact getUpper(x[-4,"Green"]) => 5.5
    @fact getValue(x[-3,"Blue"]) => -2
    @fact isequal(getVar(mcon, :lbonly),lbonly) => true
    @fact isequal(getVar(mcon, :ubonly),ubonly) => true
    @fact isequal(getVar(mcon, :onerangeub)[-7],onerangeub[-7]) => true
    @defVar(mcon, lbonly)
    @fact_throws ErrorException getVar(mcon, :lbonly)
    @fact_throws ErrorException getVar(mcon, :foo)
    d = Dict()
    @defVar(mcon, d["bar"][1:10] == 1)
    @fact getValue(d["bar"][1]) => 1
end

facts("[variable] get and set bounds") do
    m = Model()
    @defVar(m, 0 <= x <= 2)
    @fact getLower(x) => 0
    @fact getUpper(x) => 2
    setLower(x, 1)
    @fact getLower(x) => 1
    setUpper(x, 3)
    @fact getUpper(x) => 3
    @defVar(m, y, Bin)
    @fact getLower(y) => 0
    @fact getUpper(y) => 1
    @defVar(m, 0 <= y <= 1, Bin)
    @fact getLower(y) => 0
    @fact getUpper(y) => 1
    @defVar(m, fixedvar == 2)
    @fact getValue(fixedvar) => 2
    @fact getLower(fixedvar) => 2
    @fact getUpper(fixedvar) => 2
    setValue(fixedvar, 5)
    @fact getValue(fixedvar) => 5
    @fact getLower(fixedvar) => 5
    @fact getUpper(fixedvar) => 5
end

facts("[variable] get and set category") do
    m = Model()
    @defVar(m, x[1:3])
    setCategory(x[2], :Int)
    @fact getCategory(x[3]) => :Cont
    @fact getCategory(x[2]) => :Int
end

facts("[variable] repeated elements in index set (issue #199)") do
    repeatmod = Model()
    s = [:x,:x,:y]
    @defVar(repeatmod, x[s])
    @fact MathProgBase.numvar(repeatmod) => 3
end

# Test conditions in variable definition
if VERSION >= v"0.4-"
    facts("[variable] condition in indexing") do
        condmod = Model()
        @defVar(condmod, x[i=1:10; iseven(i)])
        @defVar(condmod, y[j=1:10,k=3:2:9; isodd(j+k) && k <= 8])
        @fact length(x.tupledict) => 5
        @fact length(y.tupledict) => 15
        @fact string(condmod) => "Min 0\nSubject to\n x[i] free for all i in {1,2..9,10} s.t. iseven(i)\n y[j,k] free for all j in {1,2..9,10}, k in {3,5,7,9} s.t. isodd(j + k) and k <= 8\n"
    end
end
