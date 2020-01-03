#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Test
using JuMP
using JuMP.Containers

macro dummycontainer(expr, requestedtype)
    name = gensym()
    indexvars, indexsets, condition = JuMP.Containers._parse_ref_sets(error, expr)
    if condition == :()
        return Containers.generate_container(Bool, indexvars, indexsets,
                                             requestedtype)[1]
    else
        if requestedtype != :Auto && requestedtype != :SparseAxisArray
            return :(error(""))
        end
        return Containers.generate_container(Bool, indexvars, indexsets,
                                             :SparseAxisArray)[1]
    end
end

function containermatches(c1::AbstractArray,c2::AbstractArray)
    return typeof(c1) == typeof(c2) && size(c1) == size(c2)
end

function containermatches(c1::DenseAxisArray,c2::DenseAxisArray)
    return typeof(c1) == typeof(c2) && axes(c1) == axes(c2)
end

function containermatches(c1::SparseAxisArray,
                          c2::SparseAxisArray)
    return eltype(c1) == eltype(c2)
end
containermatches(c1, c2) = false

@testset "Container syntax" begin
    @test containermatches(@dummycontainer([i=1:10], Auto), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=1:10], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=1:10], DenseAxisArray), DenseAxisArray(Vector{Bool}(undef,10), 1:10))
    @test containermatches(@dummycontainer([i=1:10], SparseAxisArray),
                           SparseAxisArray(Dict{Tuple{Any},Bool}()))

    @test containermatches(@dummycontainer([i=1:10,1:2], Auto), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=1:10,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=1:10,n=1:2], DenseAxisArray), DenseAxisArray(Matrix{Bool}(undef,10,2), 1:10, 1:2))
    @test containermatches(@dummycontainer([i=1:10,1:2], SparseAxisArray),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))

    @test containermatches(@dummycontainer([i=1:10,n=2:3], Auto), DenseAxisArray(Matrix{Bool}(undef,10,2), 1:10, 2:3))
    @test_throws ErrorException @dummycontainer([i=1:10,2:3], Array)
    @test containermatches(@dummycontainer([i=1:10,n=2:3], DenseAxisArray), DenseAxisArray(Matrix{Bool}(undef,10,2), 1:10, 2:3))
    @test containermatches(@dummycontainer([i=1:10,n=2:3], SparseAxisArray),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))


    S = Base.OneTo(10)
    @test containermatches(@dummycontainer([i=S], Auto), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], DenseAxisArray), DenseAxisArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], SparseAxisArray),
                           SparseAxisArray(Dict{Tuple{Any},Bool}()))

    @test containermatches(@dummycontainer([i=S,1:2], Auto), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,n=1:2], DenseAxisArray), DenseAxisArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], SparseAxisArray),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))

    S = 1:10
    # Not type stable to return an Array by default even when S is one-based interval
    @test containermatches(@dummycontainer([i=S], Auto), DenseAxisArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], DenseAxisArray), DenseAxisArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], SparseAxisArray),
                           SparseAxisArray(Dict{Tuple{Any},Bool}()))

    @test containermatches(@dummycontainer([i=S,n=1:2], Auto), DenseAxisArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,n=1:2], DenseAxisArray), DenseAxisArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], SparseAxisArray),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))

    # TODO: test case where S is index set not supported by DenseAxisArrays (does this exist?)

    # Conditions
    @test containermatches(@dummycontainer([i=1:10; iseven(i)], Auto),
                           SparseAxisArray(Dict{Tuple{Any},Bool}()))
    @test_throws ErrorException @dummycontainer([i=1:10; iseven(i)], Array)
    @test_throws ErrorException @dummycontainer([i=1:10; iseven(i)], DenseAxisArray)
    @test containermatches(@dummycontainer([i=1:10; iseven(i)], SparseAxisArray),
                           SparseAxisArray(Dict{Tuple{Any},Bool}()))

    # Dependent axes
    @test containermatches(@dummycontainer([i=1:10, j=1:i], Auto),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))
    @test_throws ErrorException @dummycontainer([i=1:10, j=1:i], Array)
    @test_throws ErrorException @dummycontainer([i=1:10, j=1:i], DenseAxisArray)
    @test containermatches(@dummycontainer([i=1:10, j=1:i], SparseAxisArray),
                           SparseAxisArray(Dict{NTuple{2,Any},Bool}()))

end
