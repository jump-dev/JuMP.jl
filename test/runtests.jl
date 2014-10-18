#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using JuMP
using Base.Test

tests =["print.jl",
        "expr.jl",
        "variable.jl",
        "operator.jl",
        "macros.jl",
        "model.jl",
        "probmod.jl",
        "callback.jl",
        "sosmodel.jl"]

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

#############################################################################
println(" Test: nonlinear.jl")
include("nonlinear.jl")
run_nl_tests(load_nl_solvers())

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
