#  Copyright 2018, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

#=
    This script runs all the JuMP examples.
    julia --project=examples run_examples.jl
=#

using Test

const EXAMPLES = [
    "basic.jl",
    "cannery.jl",
    "minellipse.jl",
    "optcontrol.jl",
    "qcp.jl",
    "rosenbrock.jl"
]

@testset "run_examples.jl" begin
    for example in EXAMPLES
        @testset "$(example)" begin
            example_function = include(example)
            example_function(verbose = false)
        end
    end
end
