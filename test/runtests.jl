require("MathProg")
using MathProg
using Base.Test

tests = ["quadmodel.jl",
         "expr.jl",
         "variable.jl",
         "operator.jl",
         "macros.jl",
         "model.jl"]

println("Running tests:")

for curtest in tests
  println(" Test: $(curtest)")
  include(curtest)
end

# hygiene.jl should be run separately
