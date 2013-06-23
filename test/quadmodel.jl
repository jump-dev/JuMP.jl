# This whole test only executes if we have Gurobi 
if Pkg.installed("Gurobi") != nothing
  MathProgBase.setmipsolver(Gurobi)
  modQ = Model(:Min)

  @defVar(modQ, 1.1*i <= x[i=1:3] <= 2.5*i, Int)
  
  setObjective(modQ, 3*x[1]*x[2] + x[2]*x[2] + 9*x[3]*x[1])
 
  @addConstraint(modQ, x[2] <= 1.7*x[3])
  @addConstraint(modQ, x[2] >= 0.5*x[1])
  
  status = solve(modQ)
  @test status == :Optimal
  @test_approx_eq modQ.objVal 86.4
  println(getValue(x))
  #@test_approx_eq getValue(x) 0.0
  #@test_approx_eq getValue(y) 1.0
  println(modQ.colCat)


  MathProgBase.setmipsolver(CoinMP)
end

