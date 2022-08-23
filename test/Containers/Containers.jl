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

@testset "Containers" begin
    @testset "$(file)" for file in
                           filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
        if file in ["Containers.jl"]
            continue
        end
        include(joinpath(@__DIR__, file))
    end
end
