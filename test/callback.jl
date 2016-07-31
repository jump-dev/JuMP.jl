#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
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
    @variable(mod, 0 <= x <= 2, Int)
    @variable(mod, 0 <= y <= 2, Int)
    @objective(mod, Max, y + 0.5x)
    function corners(cb)
        x_val = getvalue(x)
        y_val = getvalue(y)
        TOL = 1e-6
        # Check top right
        if y_val + x_val > 3 + TOL
            @lazyconstraint(cb, y + 0.5x + 0.5x <= 3)
        end
        entered[1] = true
        @fact_throws ErrorException @variable(cb, z)
        @fact_throws ErrorException @lazyconstraint(cb, x^2 <= 1)
    end
    addlazycallback(mod, corners)
    addlazycallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact getvalue(x) --> roughly(1.0, 1e-6)
    @fact getvalue(y) --> roughly(2.0, 1e-6)
end; end; end

facts("[callback] Test vectorized lazy constraints") do
for lazysolver in lazy_solvers
context("With solver $(typeof(lazysolver))") do
    entered = [false,false]

    mod = Model(solver=lazysolver)
    @variable(mod, 0 <= x <= 2, Int)
    @variable(mod, 0 <= y <= 2, Int)
    @objective(mod, Max, y + 0.5x)
    function corners(cb)
        x_val = getvalue(x)
        y_val = getvalue(y)
        TOL = 1e-6
        # Check top right
        if y_val + x_val > 3 + TOL
            A = [1.0 1.0
                 1.0 0.5]
            b = [3, 4]
            @lazyconstraint(cb, A*[x,y] .<= b)
        end
        entered[1] = true
        @fact_throws ErrorException @variable(cb, z)
        @fact_throws ErrorException @lazyconstraint(cb, x^2 <= 1)
    end
    addlazycallback(mod, corners)
    addlazycallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact getvalue(x) --> roughly(1.0, 1e-6)
    @fact getvalue(y) --> roughly(2.0, 1e-6)
end; end; end

facts("[callback] Test local lazy constraints") do
for lazylocalsolver in lazylocal_solvers
context("With solver $(typeof(lazylocalsolver))") do
    entered = [false,false, false]

    weights = [24714888; 21118272; 5063487; 23450813; 5598179; 4049178; 13516450; 8385365; 9684076; 31317634; 14084148; 21750211; 29261668; 17996589; 12115244]
    values = weights
    W = floor(sum(weights)/2) # knsapsack instance generated according to criteria in 'Hard knapsack Instances' (Chvatal, Op. Res., 1980), with weights composed of random integers in [1:10^(nbItems / 2)]
    mod = Model(solver=lazylocalsolver)
    @variable(mod, 0 <= x[1:length(weights)] <= 1, Int)
    @objective(mod, Max, dot(x, values))
    @constraint(mod, dot(x, weights) <=  W)

    function mycb_localzero(cb)
        nodesexpl = MathProgBase.cbgetexplorednodes(cb)
        if  entered[1] == false && nodesexpl >= 1
            # the following lazy cut  constrains all x[i] to be zero, but applies only locally at the node of the first feasible solution found: it doesn't preclude the existence of "optimal" non-trival solutions
            @lazyconstraint(cb, sum{x[i], i=1:length(weights)} <= 0, localcut=true)
            # @lazyconstraint(cb, sum{x[i], i=1:length(weights)} <= 0) # applying the cut globally would lead the solver to x=0 as the optimal solution
            @fact macroexpand(:(@lazyconstraint(cb, sum{x[i], i=1:length(weights)} <= 0, badkwarg=true))).head --> :error
            entered[1] = true
        end
        entered[2] = true
    end
    addlazycallback(mod, mycb_localzero)
    addlazycallback(mod, cb -> (entered[3] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true,true]
    @fact sum(getvalue(x)) --> greater_than(0)
end; end; end


facts("[callback] Test user cuts") do
for cutsolver in cut_solvers
context("With solver $(typeof(cutsolver))") do
    entered = [false,false]

    N = 1000
    # Include explicit data from srand(234) so that we can reproduce across platforms
    include(joinpath("data","usercut.jl"))
    mod = Model(solver=cutsolver)
    @variable(mod, x[1:N], Bin)
    @objective(mod, Max, dot(r1,x))
    @constraint(mod, c[i=1:10], dot(r2[i],x) <= rhs[i]*N/10)
    function mycutgenerator(cb)
        # add a trivially valid cut
        @usercut(cb, sum{x[i], i=1:N} <= N)
        entered[1] = true
    end
    addcutcallback(mod, mycutgenerator)
    addcutcallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getvalue(x)[:]) --> [35,38,283,305,359,397,419,426,442,453,526,553,659,751,840,865,878,978]
end; end; end

facts("[callback] Test local user cuts") do
for cutlocalsolver in cutlocal_solvers
context("With solver $(typeof(cutlocalsolver))") do
    entered = [false,false, false]

    weights = [24714888; 21118272; 5063487; 23450813; 5598179; 4049178; 13516450; 8385365; 9684076; 31317634; 14084148; 21750211; 29261668; 17996589; 12115244]
    values = weights
    W = floor(sum(weights)/2) # knsapsack instance generated according to criteria in 'Hard knapsack Instances' (Chvatal, Op. Res., 1980), with weights composed of random integers in [1:10^(nbItems / 2)]
    mod = Model(solver=cutlocalsolver)
    @variable(mod, 0 <= x[1:length(weights)] <= 1, Int)
    @objective(mod, Max, dot(x, values))
    @constraint(mod, dot(x, weights) <=  W)

    function mycb_localzero(cb)
        nodesexpl = MathProgBase.cbgetexplorednodes(cb)
        if  entered[1] == false && nodesexpl >= 1
            # the following user cut  constrains all x[i] to be zero, but applies only locally at the first node after the root node, and doesn't preclude the existence non-trival "optimal" solutions
            @usercut(cb, sum{x[i], i=1:length(weights)} <= 0, localcut=true)
            # @usercut(cb, sum{x[i], i=1:length(weights)} <= 0) # applying the cut globally would lead the solver to x=0 as the optimal solution
            @fact macroexpand(:(@usercut(cb, sum{x[i], i=1:length(weights)} <= 0, badkwarg=true))).head --> :error
            entered[1] = true
        end
        entered[2] = true
    end
    addcutcallback(mod, mycb_localzero)
    addcutcallback(mod, cb -> (entered[3] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true,true]
    @fact sum(getvalue(x)) --> greater_than(0)
end; end; end


facts("[callback] Test heuristics") do
for heursolver in heur_solvers
context("With solver $(typeof(heursolver))") do
    entered = [false,false]

    N = 100
    # Include explicit data from srand(250) so that we can reproduce across platforms
    include(joinpath("data","heuristic.jl"))
    mod = Model(solver=heursolver)
    @variable(mod, x[1:N], Bin)
    @objective(mod, Max, dot(r1,x))
    @constraint(mod, dot(ones(N),x) <= rhs*N)
    function myheuristic1(cb)
        entered[1] == true && return
        entered[1] = true
        for i in 1:100
            if i in [9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100]
                setsolutionvalue(cb, x[i], 0)
            else
                setsolutionvalue(cb, x[i], 1)
            end
        end
        addsolution(cb)
    end
    addheuristiccallback(mod, myheuristic1)
    addheuristiccallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getvalue(x)[:]) --> setdiff(1:N,[9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100])

    empty!(mod.callbacks)
    entered[1] = false
    # Test that solver rejects infeasible partial solutions...
    # ...the second solution has higher objective value, but is infeasible
    function myheuristic2(cb)
        entered[1] == true && return
        entered[1] = true
        for i in 1:90 # not every component, but close
            setsolutionvalue(cb, x[i], 1)
        end
        addsolution(cb)
    end
    addheuristiccallback(mod, myheuristic2)
    addheuristiccallback(mod, cb -> (entered[2] = true))
    @fact solve(mod) --> :Optimal
    @fact entered --> [true,true]
    @fact find(getvalue(x)[:]) --> setdiff(1:N,[9,10,11,14,15,16,25,30,32,41,44,49,50,53,54,98,100])
end; end; end

facts("[callback] Test informational callback") do
for infosolver in info_solvers
context("With solver $(typeof(infosolver))") do
    nodes      = Int[]
    objs       = Float64[]
    objvals    = Float64[]
    bestbounds = Float64[]
    entered = fill(false, 2)

    N = 10000
    include(joinpath("data","informational.jl"))
    mod = Model(solver=infosolver)
    @variable(mod, x[1:N], Bin)
    @objective(mod, Max, dot(r1,x))
    @constraint(mod, c[i=1:10], dot(r2[i],x) <= rhs[i]*N/10)
    # Test that solver fills solution correctly
    function myinfo1(cb)
        @assert MathProgBase.cbgetstate(cb) == :Intermediate
        entered[1] = true
        push!(nodes,      MathProgBase.cbgetexplorednodes(cb))
        push!(objs,       MathProgBase.cbgetobj(cb))
        push!(bestbounds, MathProgBase.cbgetbestbound(cb))
    end
    function myinfo2(cb)
        @assert MathProgBase.cbgetstate(cb) == :MIPSol
        push!(objvals, dot(r1,JuMP.getvalue(x[:])))
    end
    addinfocallback(mod, myinfo1, when = :Intermediate)
    addinfocallback(mod, myinfo2, when = :MIPSol)
    addinfocallback(mod, cb -> (entered[2] = true), when = :Intermediate)
    @fact solve(mod) --> :Optimal
    @fact entered --> fill(true, 2)
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
    for (s1,s2) in zip(objvals[1:end-1], objvals[2:end])
        @fact  s1 <= s2 --> true
    end
end; end; end

facts("[callback] Callback exit on CallbackAbort") do
for solver in lazy_solvers
context("With solver $(typeof(solver))") do
    mod = Model(solver=solver)
    @variable(mod, 0 <= x <= 2, Int)
    @variable(mod, 0 <= y <= 2, Int)
    @objective(mod, Max, x + 2y)
    @constraint(mod, y + x <= 3.5)

    mycallback = _ -> throw(CallbackAbort())
    addlazycallback(mod, mycallback)
    @fact solve(mod) --> :UserLimit
end; end; end

cbc && facts("[callback] Solver doesn't support callbacks") do
    mycb(cb) = nothing
    mod = Model(solver=Cbc.CbcSolver())
    addlazycallback(mod, mycb)
    @fact_throws ErrorException solve(mod)
    mod = Model(solver=Cbc.CbcSolver())
    addcutcallback(mod, mycb)
    @fact_throws ErrorException solve(mod)
    mod = Model(solver=Cbc.CbcSolver())
    addheuristiccallback(mod, mycb)
    @fact_throws ErrorException solve(mod)
    mod = Model(solver=Cbc.CbcSolver())
    addinfocallback(mod, mycb, when = :Intermediate)
    @fact_throws ErrorException solve(mod)
end
