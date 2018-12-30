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
    "clnlbeam.jl",
    "cluster.jl",
    "corr_sdp.jl",
    "diet.jl",
    "knapsack.jl",
    "max_cut_sdp.jl",
    "min_distortion.jl",
    "min_ellipse.jl",
    "mle.jl",
    "multi.jl",
    "prod.jl",
    "qcp.jl",
    "robust_uncertainty.jl",
    "rosenbrock.jl",
    "steelT3.jl",
    "sudoku.jl",
    "transp.jl",
    "urban_plan.jl"
]

@testset "run_examples.jl" begin
    @testset "$(example)" for example in EXAMPLES
        include(example)
    end
end
