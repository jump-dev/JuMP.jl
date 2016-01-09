#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/callback.jl
# Testing callbacks
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, MathProgBase, FactCheck

facts("[callback] Test lazy constraints") do
for lazysolver in lazy_solvers
context("With solver $(typeof(lazysolver))") do
    entered = [false,false]

    mod = Model(solver=lazysolver)
    @defVar(mod, 0 <= x <= 2, Int)
    @defVar(mod, 0 <= y <= 2, Int)
    @setObjective(mod, Max, y + 0.5x)
    function corners(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        TOL = 1e-6
        # Check top right
        if y_val + x_val > 3 + TOL
            @addLazyConstraint(cb, y + 0.5x + 0.5x <= 3)
        end
        entered[1] = true
        @fact_throws ErrorException @defVar(cb, z)
        @fact_throws ErrorException @addLazyConstraint(cb, x^2 <= 1)
    end
    addLazyCallback(mod, corners)
    addLazyCallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact getValue(x) --> roughly(1.0, 1e-6)
    @fact getValue(y) --> roughly(2.0, 1e-6)
end; end; end


facts("[callback] Test user cuts") do
for cutsolver in cut_solvers
context("With solver $(typeof(cutsolver))") do
    entered = [false,false]

    N = 1000
    # Include explicit data from srand(234) so that we can reproduce across platforms
    include(joinpath("data","usercut.jl"))
    mod = Model(solver=cutsolver)
    @defVar(mod, x[1:N], Bin)
    @setObjective(mod, Max, dot(r1,x))
    @addConstraint(mod, c[i=1:10], dot(r2[i],x) <= rhs[i]*N/10)
    function mycutgenerator(cb)
        # add a trivially valid cut
        @addUserCut(cb, sum{x[i], i=1:N} <= N)
        entered[1] = true
    end
    addCutCallback(mod, mycutgenerator)
    addCutCallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getValue(x)[:]) --> [35,38,283,305,359,397,419,426,442,453,526,553,659,751,840,865,878,978]
end; end; end


facts("[callback] Test heuristics") do
for heursolver in heur_solvers
context("With solver $(typeof(heursolver))") do
    entered = [false,false]

    N = 100
    # Include explicit data from srand(250) so that we can reproduce across platforms
    include(joinpath("data","heuristic.jl"))
    mod = Model(solver=heursolver)
    @defVar(mod, x[1:N], Bin)
    @setObjective(mod, Max, dot(r1,x))
    @addConstraint(mod, dot(ones(N),x) <= rhs*N)
    function myheuristic1(cb)
        entered[1] == true && return
        entered[1] = true
        for i in 1:100
            if i in [9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100]
                setSolutionValue!(cb, x[i], 0)
            else
                setSolutionValue!(cb, x[i], 1)
            end
        end
        addSolution(cb)
    end
    addHeuristicCallback(mod, myheuristic1)
    addHeuristicCallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getValue(x)[:]) --> setdiff(1:N,[9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100])

    empty!(mod.callbacks)
    entered[1] = false
    # Test that solver rejects infeasible partial solutions...
    # ...the second solution has higher objective value, but is infeasible
    function myheuristic2(cb)
        entered[1] == true && return
        entered[1] = true
        for i in 1:90 # not every component, but close
            setSolutionValue!(cb, x[i], 1)
        end
        addSolution(cb)
    end
    addHeuristicCallback(mod, myheuristic2)
    addHeuristicCallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getValue(x)[:]) --> setdiff(1:N,[9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100])
end; end; end

facts("[callback] Test informational callback") do
for infosolver in info_solvers
context("With solver $(typeof(infosolver))") do
    nodes      = Int[]
    objs       = Float64[]
    bestbounds = Float64[]
    entered = [false,false]

    N = 10000
    include(joinpath("data","informational.jl"))
    mod = Model(solver=infosolver)
    @defVar(mod, x[1:N], Bin)
    @setObjective(mod, Max, dot(r1,x))
    @addConstraint(mod, c[i=1:10], dot(r2[i],x) <= rhs[i]*N/10)
    # Test that solver fills solution correctly
    function myinfo(cb)
        entered[1] = true
        push!(nodes,      MathProgBase.cbgetexplorednodes(cb))
        push!(objs,       MathProgBase.cbgetobj(cb))
        push!(bestbounds, MathProgBase.cbgetbestbound(cb))
    end
    addInfoCallback(mod, myinfo)
    addInfoCallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    mono_node, mono_obj, mono_bestbound = true, true, true
    for n in 2:length(nodes)
        mono_node &= (nodes[n-1] <= nodes[n] + 1e-8)
        if nodes[n] > 0 # all bets are off at monotonicity at root node
            mono_obj &= (objs[n-1] <= objs[n] + 1e-8)
            mono_bestbound &= (bestbounds[n-1] >= bestbounds[n] - 1e-8)
        end
    end
    @fact mono_node      --> true
    @fact mono_obj       --> true
    @fact mono_bestbound --> true
end; end; end

facts("[callback] Callback exit on CallbackAbort") do
for solver in lazy_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @defVar(mod, 0 <= x <= 2, Int)
    @defVar(mod, 0 <= y <= 2, Int)
    @setObjective(mod, Max, x + 2y)
    @addConstraint(mod, y + x <= 3.5)

    mycallback = _ -> throw(CallbackAbort())
    addLazyCallback(mod, mycallback)
    @fact_throws solve(mod)
end; end; end
