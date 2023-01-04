#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestContainersNestedIterator

using JuMP.Containers
using Test

function test_NestedIterator()
    iterators = (() -> 1:3, i -> 1:i)
    condition = (i, j) -> j > i
    @test isempty(Containers.nested(iterators...; condition = condition))
    @test isempty(
        collect(Containers.nested(iterators...; condition = condition)),
    )
    condition = (i, j) -> isodd(i) || isodd(j)
    @test collect(Containers.nested(iterators...; condition = condition)) ==
          [(1, 1), (2, 1), (3, 1), (3, 2), (3, 3)]
    return
end

function test_StackOverflow_2335()
    iterator = Containers.nested(() -> 1:100_000; condition = _ -> false)
    @test iterate(iterator) === nothing
    return
end

end  # module
