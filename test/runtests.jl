require("MathProg")
using MathProg

tests = ["model.jl",
         "variable.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(my_test)
end
