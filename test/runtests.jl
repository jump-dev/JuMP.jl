require("JuMP")
using JuMP
using Base.Test

tests = ["expr.jl",
         "variable.jl",
         "operator.jl",
         "macros.jl",
         "model.jl"]

println("Running tests:")

for curtest in tests
  println(" Test: $(curtest)")
  include(curtest)
end

if Pkg.installed("Gurobi") != nothing
  gurobitests = ["qcqpmodel.jl", "quadmodel.jl"]
  for curtest in gurobitests
    println(" Test: $(curtest)")
    include(curtest)
  end
else
  println("WARNING: Gurobi not installed, cannot execute corresponding tests")
end

# hygiene.jl should be run separately
