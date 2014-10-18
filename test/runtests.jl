using JuMP
using Base.Test

tests =["print.jl",
        "expr.jl",
        "variable.jl",
        "operator.jl",
        "macros.jl",
        #"model.jl",
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
nl_solvers = Any[]
if Pkg.installed("Ipopt") != nothing
    eval(Expr(:import,:Ipopt))
    push!(nl_solvers, ("Ipopt",Ipopt.IpoptSolver(print_level=0)))
end
if Pkg.installed("NLopt") != nothing
    eval(Expr(:import,:NLopt))
    push!(nl_solvers, ("NLopt",NLopt.NLoptSolver(algorithm=:LD_SLSQP)))
end
if Pkg.installed("KNITRO") != nothing
    eval(Expr(:import,:KNITRO))
    push!(nl_solvers, ("KNITRO",KNITRO.KnitroSolver()))
end
include("nonlinear.jl")
run_nl_tests(nl_solvers)

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
