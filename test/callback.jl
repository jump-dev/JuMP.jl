

function callback_test(solvername, lazysolver, cutsolver)
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
                @addLazyConstraint(cb, y + x <= 3)
            end
        end
        setLazyCallback(mod, corners)
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
    end  # lazy let
    # User cuts
    println(string("  Running ", solvername, " user"))
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


end

if Pkg.installed("Gurobi") != nothing  
    using Gurobi
    callback_test("Gurobi", GurobiSolver(LazyConstraints=1, OutputFlag=0), GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0)) 
end
if Pkg.installed("CPLEXLink") != nothing
    using CPLEXLink
    callback_test("CPLEXLink", CplexSolver(), CplexSolver()) 
end
if Pkg.installed("GLPKMathProgInterface") != nothing
    using GLPKMathProgInterface
    callback_test("GLPK", GLPKSolverMIP(), GLPKSolverMIP())
end
