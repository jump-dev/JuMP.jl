# This whole test only executes if we have Gurobi 
if Pkg.installed("Gurobi") != nothing
  
  let
    modQ = Model(:Min,lpsolver=LPSolver(:Gurobi))
    
    @defVar(modQ, -2 <= x <= 2 )
    @defVar(modQ, -2 <= y <= 2 )
    
    @setObjective(modQ, x - y )
    addConstraint(modQ, x + x*x + x*y + y*y <= 1 )
    
    status = solve(modQ)
    @test status == :Optimal
    @test_approx_eq_eps modQ.objVal -1-4/sqrt(3) 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1/3 1e-3
  end

  let
    modQ = Model(:Min,mipsolver=MIPSolver(:Gurobi))
    
    @defVar(modQ, -2 <= x <= 2, Int )
    @defVar(modQ, -2 <= y <= 2, Int )
    
    @setObjective(modQ, x - y )
    addConstraint(modQ, x + x*x + x*y + y*y <= 1 )
    
    status = solve(modQ)
    @test status == :Optimal
    @test_approx_eq_eps modQ.objVal -3 1e-6
    @test_approx_eq_eps (getValue(x) + getValue(y)) -1 1e-6
  end

else 
  println("WARNING: Gurobi not installed, cannot execute QCQP test")
end
