#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

using Test

import PackageCompiler

PackageCompiler.create_app(
    joinpath(@__DIR__, "MyApp"),
    joinpath(@__DIR__, "build");
    force = true,
)

@testset "test --solver=ipopt" begin
    app = joinpath(@__DIR__, "build", "bin", "MyApp")
    output = sprint(io -> run(pipeline(`$app --solver=highs`; stdout = io)))
    @test occursin("HiGHS", output)
    @test occursin("OPTIMAL", output)
end

@testset "test --solver=ipopt" begin
    app = joinpath(@__DIR__, "build", "bin", "MyApp")
    output = sprint(io -> run(pipeline(`$app --solver=ipopt`; stdout = io)))
    @test occursin("Ipopt", output)
    @test occursin("LOCALLY_SOLVED", output)
end
