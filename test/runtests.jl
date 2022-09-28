#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# test/runtests.jl
#############################################################################

using Test
import JuMP

@testset "Ambiguities" begin
    # It is important to run these first, before the tests start adding methods.
    @test isempty(Test.detect_ambiguities(JuMP))
    # TODO(odow): there are still some ambiguities in Containers
    # @test isempty(Test.detect_ambiguities(JuMP.Containers))
end

t = time()
include("Containers/Containers.jl")
println("Containers.jl took $(round(time() - t; digits = 1)) seconds.")

for file in filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
    if file in [
        "runtests.jl",
        "utilities.jl",
        "JuMPExtension.jl",
        "nlp_solver.jl",
        "hygiene.jl",
    ]
        continue
    end

    @testset "$(file)" begin
        t = time()
        include(file)
        println("$(file) took $(round(time() - t; digits = 1)) seconds.")
    end
end

# TODO: The hygiene test should run in a separate Julia instance where JuMP
# hasn't been loaded via `using`.
include("hygiene.jl")
