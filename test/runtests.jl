#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using JuMP

using LinearAlgebra  # for dot and tr
using SparseArrays # for sparse
using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

include("Containers/Containers.jl")

include("utilities.jl")
include("JuMPExtension.jl")

include("derivatives.jl")
include("derivatives_coloring.jl")
include("model.jl")
include("variable.jl")
include("expr.jl")
include("objective.jl")
include("constraint.jl")
include("nlp.jl")
include("generate_and_solve.jl")
include("print.jl")
include("operator.jl")
include("macros.jl")
include("lp_sensitivity.jl")
# TODO: The hygiene test should run in a separate Julia instance where JuMP hasn't been loaded via `using`.
include("hygiene.jl")
