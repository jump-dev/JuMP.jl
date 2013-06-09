require("MathProg")
using MathProg
using Base.Test

tests = ["affexpr.jl",
         "quadexpr.jl",
         "variable.jl"]

println("Running tests:")

for curtest in tests
  println(" Test: $(curtest)")
  include(curtest)
end
