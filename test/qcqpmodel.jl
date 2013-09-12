# This whole test only executes if we have Gurobi 
if Pkg.installed("Gurobi") != nothing
  
  cursolver = MathProgBase.lpsolver
  MathProgBase.setlpsolver(:Gurobi)
  let
    modQ = Model(:Min)
    
    @defVar(modQ, -2 <= x <= 2 )
    @defVar(modQ, -2 <= y <= 2 )
    
    @setObjective(modQ, x - y )
    addConstraint(modQ, x + x*x + x*y + y*y <= 1 )
    
    status = solve(modQ)
    @test status == :Optimal
    @test_approx_eq_eps modQ.objVal -1-4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
  end

  MathProgBase.setlpsolver(cursolver)

  cursolver = MathProgBase.mipsolver
  MathProgBase.setmipsolver(:Gurobi)
  let
    modQ = Model(:Min)
    
    @defVar(modQ, -2 <= x <= 2, Int )
    @defVar(modQ, -2 <= y <= 2, Int )
    
    @setObjective(modQ, x - y )
    addConstraint(modQ, x + x*x + x*y + y*y <= 1 )
    
    status = solve(modQ)
    @test status == :Optimal
    @test_approx_eq_eps modQ.objVal -3 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1 1e-6
  end

  MathProgBase.setmipsolver(cursolver)
else 
  println("WARNING: Gurobi not installed, cannot execute QCQP test")
end
