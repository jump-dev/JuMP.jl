#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/callback.jl
# Testing callbacks
# Must be run as part of runtests.jl, as it needs a list of solvers.
#############################################################################
using JuMP, FactCheck

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
    setLazyCallback(mod, corners)
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
    setCutCallback(mod, mycutgenerator)
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
    function myheuristic(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        # Heuristic is to round solution down
        setSolutionValue!(cb, x, floor(x_val))
        # Leave y undefined - solver should handle as it sees fit
        # In case of Gurobi - try to figure out what it should be
        addSolution(cb)
        # Check that solvers ignore infeasible solutions
        setSolutionValue!(cb, x, 3)
        addSolution(cb)
        setSolutionValue!(cb, x, 3)
        setSolutionValue!(cb, y, 5)
        addSolution(cb)
    end
    setHeuristicCallback(mod, myheuristic)
    @fact solve(mod) => :Optimal
    @fact getValue(x) => roughly(1.0, 1e-6)
    @fact getValue(y) => roughly(2.0, 1e-6)
end; end; end