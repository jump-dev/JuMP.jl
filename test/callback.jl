#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/callback.jl
# Testing callbacks
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, MathProgBase, FactCheck, Compat

facts("[callback] Test lazy constraints") do
for lazysolver in lazy_solvers
context("With solver $(typeof(lazysolver))") do

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
    end
    addLazyCallback(mod, corners)
    @fact solve(mod) => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)
end; end; end


facts("[callback] Test user cuts") do
for cutsolver in cut_solvers
context("With solver $(typeof(cutsolver))") do

    mod = Model(solver=cutsolver)
    @defVar(mod, 0 <= x <= 2, Int)
    @defVar(mod, 0 <= y <= 2, Int)
    @setObjective(mod, Max, x + 2y)
    @addConstraint(mod, y + x <= 3.5)
    function mycutgenerator(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        TOL = 1e-6
        # Check top right
        if y_val + x_val > 3 + TOL
            @addUserCut(cb, y + x <= 3)
        end
    end
    addCutCallback(mod, mycutgenerator)
    @fact solve(mod) => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)
end; end; end


facts("[callback] Test heuristics") do
for heursolver in heur_solvers
context("With solver $(typeof(heursolver))") do

    mod = Model(solver=heursolver)
    @defVar(mod, 0 <= x <= 2, Int)
    @defVar(mod, 0 <= y <= 2, Int)
    @setObjective(mod, Max, x + 2y)
    @addConstraint(mod, y + x <= 3.5)
    # Test that solver fills solution correctly
    function myheuristic1(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        # Heuristic is to round solution down
        setSolutionValue!(cb, x, floor(x_val))
        # Leave y undefined - solver should handle as it sees fit
        # In case of Gurobi - try to figure out what it should be
        addSolution(cb)
    end
    addHeuristicCallback(mod, myheuristic1)
    @fact solve(mod) => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)

    empty!(mod.callbacks)
    # Test that solver rejects infeasible partial solutions...
    function myheuristic2(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        setSolutionValue!(cb, x, 3)
        addSolution(cb)
    end
    addHeuristicCallback(mod, myheuristic2)
    @fact solve(mod) => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)
end; end; end

facts("[callback] Test informational callback") do
for infosolver in info_solvers
context("With solver $(typeof(infosolver))") do
    nodes      = Int[]
    objs       = Float64[]
    bestbounds = Float64[]

    srand(100)
    N = 10000
    mod = Model(solver=infosolver)
    @defVar(mod, x[1:N], Bin)
    @setObjective(mod, Max, dot(rand(N),x))
    @addConstraint(mod, c[1:10], dot(rand(N),x) <= rand()*N/10)
    # Test that solver fills solution correctly
    function myinfo(cb)
        push!(nodes,      MathProgBase.cbgetexplorednodes(cb))
        push!(objs,       MathProgBase.cbgetobj(cb))
        push!(bestbounds, MathProgBase.cbgetbestbound(cb))
    end
    addInfoCallback(mod, myinfo)
    @fact solve(mod) => :Optimal
    mono_node, mono_obj, mono_bestbound = true, true, true
    for n in 2:length(nodes)
        mono_node &= (nodes[n-1] <= nodes[n] + 1e-8)
        if nodes[n] > 0 # all bets are off at monotonicity at root node
            mono_obj &= (objs[n-1] <= objs[n] + 1e-8)
            mono_bestbound &= (bestbounds[n-1] >= bestbounds[n] - 1e-8)
        end
    end
    @fact mono_node      => true
    @fact mono_obj       => true
    @fact mono_bestbound => true
end; end; end

facts("[callback] Test multiple callbacks") do
for solver in intersect(lazy_solvers,cut_solvers,heur_solvers)
context("With solver $(typeof(solver))") do

    mod = Model(solver=solver)
    @defVar(mod, 0 <= x <= 2, Int)
    @defVar(mod, 0 <= y <= 2, Int)
    @setObjective(mod, Max, x + 2y)
    @addConstraint(mod, y + x <= 3.5)
    cb_tracker = @compat Dict(
        :l_1 => false,
        :l_2 => false,
        :c_1 => false,
        :c_2 => false,
        :h_1 => false,
        :h_2 => false
    )
    addLazyCallback(mod, cb -> (cb_tracker[:l_1] = true))
    addLazyCallback(mod, cb -> (cb_tracker[:l_2] = true))
    addCutCallback(mod, cb -> (cb_tracker[:c_1] = true))
    addCutCallback(mod, cb -> (cb_tracker[:c_2] = true))
    addHeuristicCallback(mod, cb -> (cb_tracker[:h_1] = true))
    addHeuristicCallback(mod, cb -> (cb_tracker[:h_2] = true))

    @fact solve(mod) => :Optimal
    @fact collect(values(cb_tracker)) => fill(true,6)
end; end; end

