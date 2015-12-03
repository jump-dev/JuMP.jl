#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using JuMP
using FactCheck
using Base.Test
#FactCheck.setstyle(:compact)

# Static tests - don't require a solver
include("print.jl")
include("variable.jl")
include("expr.jl")
include("operator.jl")
include("macros.jl")

# Fuzzer of macros to build expressions
include("fuzzer.jl")

# Load solvers
include("solvers.jl")

# Solver-dependent tests
include("model.jl");        length(   lp_solvers) == 0 && warn("Model tests not run!")
include("probmod.jl");      length(   lp_solvers) == 0 && warn("Prob. mod. tests not run!")
include("callback.jl");     length( lazy_solvers) == 0 && warn("Callback tests not run!")
include("qcqpmodel.jl");    length( quad_solvers) == 0 && warn("Quadratic tests not run!")
include("nonlinear.jl");    length(  nlp_solvers) == 0 && warn("Nonlinear tests not run!")
                            length(minlp_solvers) == 0 && warn("Mixed-integer Nonlinear tests not run!")
include("sdp.jl");          length(  sdp_solvers) == 0 && warn("Semidefinite tests not run!")
include("socduals.jl");     length(conic_solvers_with_duals) == 0 && warn("Conic solvers with duals tests not run!")
# Throw an error if anything failed
FactCheck.exitstatus()

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
