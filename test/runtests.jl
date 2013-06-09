require("MathProg")
using MathProg

tests = ["affexpr.jl",
         "variable.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end
