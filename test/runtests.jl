require("JuMP")
using JuMP
using Base.Test

tests = ["expr.jl",
         "variable.jl",
         "operator.jl",
         "macros.jl",
         "model.jl",
         "probmod.jl",
         "callback.jl",
         "sosmodel.jl",
         "vector_validity.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end

if Pkg.installed("Gurobi") != nothing || 
   Pkg.installed("CPLEX") != nothing ||
   Pkg.installed("Mosek") != nothing
    quadtests = ["qcqpmodel.jl", "quadmodel.jl"]
    for curtest in quadtests
        println(" Test: $(curtest)")
        include(curtest)
    end
else
    println("WARNING: Neither Gurobi nor CPLEX nor Mosek installed, cannot execute corresponding tests")
end

if Pkg.installed("Ipopt") != nothing
    println(" Test: nonlinear.jl")
    include("nonlinear.jl")
end

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
