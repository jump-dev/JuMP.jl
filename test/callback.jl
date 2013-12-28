# Pull packages in that are available
# TODO: Put this inside macro
if Pkg.installed("Gurobi") != nothing  
  using Gurobi
end
if Pkg.installed("CPLEXLink") != nothing
  using CPLEXLink
end
if Pkg.installed("GLPKMathProgInterface") != nothing
  using GLPKMathProgInterface
end

# Only write test code once, use macros to copy-paste it
macro callback_test(solvername)
  if solvername == :Gurobi
    pkg_install         = :(Pkg.installed("Gurobi"))
    mod_creation_lazy   = :(mod = Model(solver=GurobiSolver(LazyConstraints=1, OutputFlag=0)))
    mod_creation_cut    = :(mod = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0, OutputFlag=0)))
  elseif solvername == :CPLEXLink
    pkg_install         = :(Pkg.installed("CPLEXLink"))
    mod_creation_lazy   = :(mod = Model(solver=CplexSolver()))
    mod_creation_cut    = :(mod = Model(solver=CplexSolver()))
  elseif solvername == :GLPKMathProgInterface
    pkg_install         = :(Pkg.installed("GLPKMathProgInterface"))
    mod_creation_lazy   = :(mod = Model(solver=GLPKSolverMIP()))
    mod_creation_cut    = :(mod = Model(solver=GLPKSolverMIP()))
  end

  quote
    if $pkg_install != nothing  
      # Lazy Constraints
      println(string("  Running ", $solvername, " lazy"))
      let
        $mod_creation_lazy
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
        setlazycallback(mod, corners)
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
      end  # lazy let
      # User cuts
      println(string("  Running ", $solvername, " user"))
      let
        $mod_creation_cut
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
        setcutcallback(mod, mycutgenerator)
        solve(mod)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
        @test_approx_eq_eps getValue(y) 2.0 1e-6
      end  # cut let
    end  # pkg test
  end  # quote
end  # macro

# Generate the three tests
@callback_test Gurobi
@callback_test CPLEXLink
@callback_test GLPKMathProgInterface
