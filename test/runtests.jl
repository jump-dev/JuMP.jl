#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using JuMP
using FactCheck
#FactCheck.setstyle(:compact)

# Static tests - don't require a solver
include("print.jl")
include("variable.jl")
include("expr.jl")
include("operator.jl")
include("macros.jl")

# Load solvers
include("solvers.jl")

# Solver-dependent tests
include("model.jl");        length(   lp_solvers) == 0 && warn("Model tests not run!")
include("probmod.jl");      length(   lp_solvers) == 0 && warn("Prob. mod. tests not run!")
include("callback.jl");     length( lazy_solvers) == 0 && warn("Callback tests not run!")
include("qcqpmodel.jl");    length( quad_solvers) == 0 && warn("Quadratic tests not run!")
include("nonlinear.jl");    length(  nlp_solvers) == 0 && warn("Nonlinear tests not run!")
                            length(minlp_solvers) == 0 && warn("Mixed-integer Nonlinear tests not run!")

# Throw an error if anything failed
FactCheck.exitstatus()

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
