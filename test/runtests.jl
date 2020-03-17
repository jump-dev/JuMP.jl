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

include("Containers/Containers.jl")

include("utilities.jl")
include("JuMPExtension.jl")

@testset "$(file)" for file in filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
    if file in [
        "runtests.jl",
        "utilities.jl",
        "JuMPExtension.jl",
        "nlp_solver.jl",
        "hygiene.jl",
    ]
        continue
    end
    include(file)
end

# TODO: The hygiene test should run in a separate Julia instance where JuMP
# hasn't been loaded via `using`.
include("hygiene.jl")
