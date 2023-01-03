#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestContainersVectorizedProductIterator

using JuMP.Containers
using Test

function test_VectorizedProductIterator()
    I = [1 2; 3 4]
    @test isempty(Containers.vectorized_product(2, I, 1:0))
    @test isempty(collect(Containers.vectorized_product(2, I, 1:0)))
    @test collect(Containers.vectorized_product(2, I)) ==
          [(2, 1) (2, 3) (2, 2) (2, 4)]
    return
end

function test_unknown_size()
    f = Iterators.filter(k -> isodd(k), 1:10)
    v = Containers.vectorized_product(f)
    @test axes(v) == (Base.OneTo(5),)
    return
end

function test_infinite_size()
    f = Iterators.repeated(1)
    @test_throws ErrorException Containers.vectorized_product(f)
    return
end

end  # module
