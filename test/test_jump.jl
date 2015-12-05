using ReverseDiffSparse2
using JuMP, Ipopt

# need to support ifelse() for nlp_solvers
nlp_solvers = [] #[RDSSolver(IpoptSolver())]

# All tests pass except the unboundedness test
convex_nlp_solvers = [RDSSolver(IpoptSolver(print_level=0))]
minlp_solvers = []

include(Pkg.dir("JuMP","test","nonlinear.jl"))
