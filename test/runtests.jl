#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using JuMP

using Compat
using Compat.LinearAlgebra  # for dot and tr
using Compat.SparseArrays # for sparse
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

# TODO: This should be defined in JuMP source with the other versions of
# isequal_canonical.
function JuMP.isequal_canonical(x::AbstractArray{<:JuMP.AbstractJuMPScalar}, y::AbstractArray{<:JuMP.AbstractJuMPScalar})
    size(x) == size(y) && all(JuMP.isequal_canonical.(x, y))
end

macro test_expression(expr)
    esc(quote
            @test JuMP.isequal_canonical(@expression(m, $expr), $expr)
    end)
end

macro test_expression_with_string(expr, str)
    esc(quote
            @test string(@inferred $expr) == $str
            @test_expression $expr
    end)
end

# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(errortype, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    :(@test_throws $errortype try @eval $m catch err; throw(err.error) end)
end

# include("JuMPExtension.jl")

# include("derivatives.jl")
# include("derivatives_coloring.jl")

# include("containers.jl")
# include("model.jl")
# include("variable.jl")
# include("expr.jl")
# include("objective.jl")
# include("constraint.jl")
# include("nlp.jl")
# include("generate_and_solve.jl")
# include("print.jl")
include("operator.jl")
# include("macros.jl")

# Fuzzer of macros to build expressions
#include("fuzzer.jl")

# Solver-dependent tests
#include("model.jl");        length(   lp_solvers) == 0 && warn("Model tests not run!")
#include("probmod.jl");      length(   lp_solvers) == 0 && warn("Prob. mod. tests not run!")
#include("callback.jl");     length( lazy_solvers) == 0 && warn("Callback tests not run!")
#include("qcqpmodel.jl");    length( quad_solvers) == 0 && warn("Quadratic tests not run!")
#include("nonlinear.jl");    length(  nlp_solvers) == 0 && warn("Nonlinear tests not run!")
#                            length(minlp_solvers) == 0 && warn("Mixed-integer Nonlinear tests not run!")
#include("sdp.jl");          length(  sdp_solvers) == 0 && warn("Semidefinite tests not run!")
#include("socduals.jl");     length(conic_solvers_with_duals) == 0 && warn("Conic solvers with duals tests not run!")

# hygiene.jl should be run separately
# hockschittkowski/runhs.jl has additional nonlinear tests
