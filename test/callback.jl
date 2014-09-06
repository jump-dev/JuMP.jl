function callback_test(solvername, lazysolver, cutsolver, heursolver)
    # Lazy Constraints
    println(string("  Running ", solvername, " lazy"))
    let
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
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
    end  # lazy let

    # User cuts
    println(string("  Running ", solvername, " cuts"))
    let
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
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
    end  # cut let


    # User heuristic
    println(string("  Running ", solvername, " heuristic"))
    let
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
        end
        setHeuristicCallback(mod, myheuristic)
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
    end  # cut let
end

if Pkg.installed("Gurobi") != nothing  
    using Gurobi
    callback_test("Gurobi", GurobiSolver(OutputFlag=0), 
                            GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0),
                            GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0))
end
if Pkg.installed("CPLEX") != nothing
    using CPLEX
    callback_test("CPLEX", CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0),
                           CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0),
                           CplexSolver(CPX_PARAM_PRELINEAR=0, CPX_PARAM_PREIND=0, CPX_PARAM_ADVIND=0, CPX_PARAM_MIPSEARCH=1,CPX_PARAM_MIPCBREDLP=0)) 
end
if Pkg.installed("GLPKMathProgInterface") != nothing
    using GLPKMathProgInterface
    callback_test("GLPK", GLPKSolverMIP(), GLPKSolverMIP(), GLPKSolverMIP())
end
