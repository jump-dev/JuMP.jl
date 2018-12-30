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
    "corr_sdp.jl",
    "mindistortion.jl",
    "knapsack.jl",
    "maxcut_sdp.jl",
    "minellipse.jl",
    "mle.jl",
    "optcontrol.jl",
    "qcp.jl",
    "robust_uncertainty.jl",
    "rosenbrock.jl",
    "sudoku.jl",
    "urbanplan.jl"
]

@testset "run_examples.jl" begin
    for example in EXAMPLES
        @testset "$(example)" begin
            include(example)
        end
    end
end
