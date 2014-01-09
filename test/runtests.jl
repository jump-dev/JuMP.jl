require("JuMP")
using JuMP
using Base.Test

tests = ["expr.jl",
         "variable.jl",
         "operator.jl",
         "macros.jl",
         "model.jl",
         "probmod.jl",
         "callback.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end

if Pkg.installed("Gurobi") != nothing || Pkg.installed("CPLEXLink") != nothing
    quadtests = ["qcqpmodel.jl", "quadmodel.jl"]
    for curtest in quadtests
        println(" Test: $(curtest)")
        include(curtest)
    end
else
    println("WARNING: Gurobi or CPLEXLink not installed, cannot execute corresponding tests")
end

# hygiene.jl should be run separately
