#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# urbanplan.jl
#
# An "urban planning" problem.
# Based on 
#  http://puzzlor.editme.com/urbanplanning
#############################################################################

using JuMP

function SolveUrban()

  m = Model(:Max)

  # x is indexed by row and column
  @defVar(m, 0 <= x[1:5,1:5] <= 1, Int)

  # y is indexed by R or C, and the points
  # JuMP allows indexing on arbitrary sets
  rowcol = ["R","C"]
  points = [+5,+4,+3,-3,-4,-5]
  @defVar(m, 0 <= y[rowcol,points,i=1:5] <= 1, Int)

  # Objective - combine the positive and negative parts
  @setObjective(m, sum{
    3*(y["R", 3,i] + y["C", 3,i]) 
  + 1*(y["R", 4,i] + y["C", 4,i]) 
  + 1*(y["R", 5,i] + y["C", 5,i]) 
  - 3*(y["R",-3,i] + y["C",-3,i]) 
  - 1*(y["R",-4,i] + y["C",-4,i])
  - 1*(y["R",-5,i] + y["C",-5,i]), i=1:5})

  # Constrain the number of residential lots
  @addConstraint(m, sum{x[i,j], i=1:5, j=1:5} == 12)

  # Add the constraints that link the auxiliary y variables
  # to the x variables
  # Rows
  for i = 1:5
    @addConstraint(m, y["R", 5,i] <=   1/5*sum{x[i,j], j=1:5}) # sum = 5
    @addConstraint(m, y["R", 4,i] <=   1/4*sum{x[i,j], j=1:5}) # sum = 4
    @addConstraint(m, y["R", 3,i] <=   1/3*sum{x[i,j], j=1:5}) # sum = 3
    @addConstraint(m, y["R",-3,i] >= 1-1/3*sum{x[i,j], j=1:5}) # sum = 2
    @addConstraint(m, y["R",-4,i] >= 1-1/2*sum{x[i,j], j=1:5}) # sum = 1
    @addConstraint(m, y["R",-5,i] >= 1-1/1*sum{x[i,j], j=1:5}) # sum = 0
  end
  # Columns
  for j = 1:5
    @addConstraint(m, y["C", 5,j] <=   1/5*sum{x[i,j], i=1:5}) # sum = 5
    @addConstraint(m, y["C", 4,j] <=   1/4*sum{x[i,j], i=1:5}) # sum = 4
    @addConstraint(m, y["C", 3,j] <=   1/3*sum{x[i,j], i=1:5}) # sum = 3
    @addConstraint(m, y["C",-3,j] >= 1-1/3*sum{x[i,j], i=1:5}) # sum = 2
    @addConstraint(m, y["C",-4,j] >= 1-1/2*sum{x[i,j], i=1:5}) # sum = 1
    @addConstraint(m, y["C",-5,j] >= 1-1/1*sum{x[i,j], i=1:5}) # sum = 0
  end

  # Solve it with the default solver (CBC)
  status = solve(m)
  if status == :Infeasible
    error("Solver couldn't find solution!")
  end

  # Print results
  println("Best objective: $(round(getObjectiveValue(m)))")
  println(getValue(x))
end

SolveUrban()
