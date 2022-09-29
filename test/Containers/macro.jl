#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

using JuMP
using JuMP.Containers
using Test

@testset "Macro" begin
    @testset "Array" begin
        Containers.@container(x[i = 1:3], i^2)
        @test x isa Vector{Int}
        x = Containers.@container([i = 1:3, j = 1:3], i^2)
        @test x isa Matrix{Int}
        # alternative syntax
        Containers.@container(x[i in 1:3], i^2)
        @test x isa Vector{Int}
        Containers.@container(x[i ∈ 1:3], i^2)
        @test x isa Vector{Int}
    end
    @testset "Forced array" begin
        set = 1:2
        x = Containers.@container([i = set, j = 1:2], i + j, container = Array)
        @test x == [2 3; 3 4]
        set_2 = [1, 2]
        y = Containers.@container(
            [i = set_2, j = 1:2],
            i + j,
            container = Array
        )
        @test y == [2 3; 3 4]
        @test_throws(
            ErrorException,
            Containers.@container(
                [i = [:a, :b], j = 1:2],
                i + j,
                container = Array
            )
        )
        @test_throws(
            ErrorException,
            Containers.@container(
                [i = 1:2:4, j = 1:2],
                i + j,
                container = Array
            )
        )
    end
    @testset "DenseAxisArray" begin
        Containers.@container(x[i = 2:3], i^2)
        @test x isa Containers.DenseAxisArray{Int,1}
        Containers.@container(x[i = 2:3, j = 1:2], i + j)
        @test x isa Containers.DenseAxisArray{Int,2}
        Containers.@container(x[4], 0.0)
        @test x isa Containers.DenseAxisArray{Float64,1}
        Containers.@container(x[4, 5], 0)
        @test x isa Containers.DenseAxisArray{Int,2}
        Containers.@container(x[4, 1:3, 5], 0)
        @test x isa Containers.DenseAxisArray{Int,3}
    end
    @testset "SparseAxisArray" begin
        Containers.@container(x[i = 1:3, j = 1:i], i + j)
        @test x isa Containers.SparseAxisArray{Int,2,Tuple{Int,Int}}
        Containers.@container(x[i = 1:10; iseven(i)], i)
        @test x isa Containers.SparseAxisArray{Int,1,Tuple{Int}}
        # Return types are not the same across Julia versions. Check for a
        # variety of plausible results.
        T = Union{Tuple{Any,Any},Tuple{Int,Any},Tuple{Int,Int}}
        # Here the iterators are empty, and have a linked dependence, so we
        # can't always infer the key or value types.
        Containers.@container(x[i = 1:0, j = i:0], i)
        @test(
            x isa SparseAxisArray{Int,2,<:T} ||
            x isa SparseAxisArray{Any,2,<:T}
        )
        # This one is better, we can infer the value type, but the keys are
        # difficult to infer, and depend on the Julia version you are running.
        Containers.@container(x[i = 1:2, j = 1:2; false], i)
        @test x isa SparseAxisArray{Int,2,<:T}
        Containers.@container(
            x[i = 1:0, j = 2:1],
            i,
            container = SparseAxisArray
        )
        @test x isa SparseAxisArray{Int,2,Tuple{Int,Int}}
        Containers.@container(
            x[i = 1:0, j = 1:0],
            i,
            container = SparseAxisArray
        )
        @test x isa SparseAxisArray{Int,2,Tuple{Int,Int}}
    end
    @testset "duplicate_indices" begin
        expr = :(Containers.@container(x[i = 1:2, i = 1:2], i + i))
        @test_throws(LoadError, try
            @eval $expr
        catch err
            throw(err)
        end,)
    end
    @testset "double_filter_typed_vcat" begin
        expr =
            :(Containers.@container(x[i = 1:2; isodd(i); iseven(i + 1)], i + i))
        @test_throws(LoadError, try
            @eval $expr
        catch err
            throw(err)
        end,)
    end
    @testset "double_filter_vect" begin
        expr = :(Containers.@container(
            x[i = 1:2, 1:2; isodd(i); iseven(i + 1)],
            i + i,
        ))
        @test_throws(LoadError, try
            @eval $expr
        catch err
            throw(err)
        end,)
    end
    @testset "Dict" begin
        Containers.@container(v[i = 1:3], sin(i), container = Dict)
        @test v isa Dict{Int,Float64}
        @test length(v) == 3
        @test v[2] ≈ sin(2)
        Containers.@container(w[i = 1:3, j = 1:3], i + j, container = Dict)
        @test w isa Dict{Tuple{Int,Int},Int}
        @test length(w) == 9
        @test w[2, 3] == 5
        Containers.@container(
            x[i = 1:3, j = [:a, :b]],
            (j, i),
            container = Dict
        )
        @test x isa Dict{Tuple{Int,Symbol},Tuple{Symbol,Int}}
        @test length(x) == 6
        @test x[2, :a] == (:a, 2)
        Containers.@container(y[i = 1:3, j = 1:i], i + j, container = Dict)
        @test y isa Dict{Tuple{Int,Int},Int}
        @test length(y) == 6
        @test y[2, 1] == 3
        Containers.@container(
            z[i = 1:3, j = 1:3; isodd(i + j)],
            i + j,
            container = Dict
        )
        @test z isa Dict{Tuple{Int,Int},Int}
        @test length(z) == 4
        @test z[1, 2] == 3
    end
    @testset "Invalid container" begin
        err = ErrorException(
            "Unable to build a container with the provided type $(Int). " *
            "Implement `Containers.container(::Function, indices, ::Type{$Int})`.",
        )
        @test_throws err Containers.@container(
            x[i = 1:2, j = 1:2],
            i + j,
            container = Int
        )
    end
    @testset "Compound indexing expressions" begin
        Containers.@container(
            x[(i, j) in [(1, 1), (2, 2)], k in i:3],
            i + j + k,
        )
        @test x isa Containers.SparseAxisArray
        @test length(x) == 5
        @test x[(2, 2), 3] == 7
    end
end

struct _MyContainer end

function Containers.container(f::Function, indices, ::Type{_MyContainer})
    key(i::Tuple) = i
    key(i::Tuple{T}) where {T} = i[1]
    return Dict(key(i) => f(i...) for i in indices)
end

@testset "_MyContainer" begin
    Containers.@container(v[i = 1:3], sin(i), container = _MyContainer)
    @test v isa Dict{Int,Float64}
    @test length(v) == 3
    @test v[2] ≈ sin(2)
    Containers.@container(w[i = 1:3, j = 1:3], i + j, container = _MyContainer)
    @test w isa Dict{Tuple{Int,Int},Int}
    @test length(w) == 9
    @test w[2, 3] == 5
    Containers.@container(
        x[i = 1:3, j = [:a, :b]],
        (j, i),
        container = _MyContainer
    )
    @test x isa Dict{Tuple{Int,Symbol},Tuple{Symbol,Int}}
    @test length(x) == 6
    @test x[2, :a] == (:a, 2)
    Containers.@container(y[i = 1:3, j = 1:i], i + j, container = _MyContainer)
    @test y isa Dict{Tuple{Int,Int},Int}
    @test length(y) == 6
    @test y[2, 1] == 3
    Containers.@container(
        z[i = 1:3, j = 1:3; isodd(i + j)],
        i + j,
        container = _MyContainer
    )
    @test z isa Dict{Tuple{Int,Int},Int}
    @test length(z) == 4
    @test z[1, 2] == 3
end

# Test containers that use subindex names
struct _MyContainer2
    names::Any
    d::Any
end

function Containers.container(
    f::Function,
    indices,
    ::Type{_MyContainer2},
    names,
)
    key(i::Tuple) = i
    key(i::Tuple{T}) where {T} = i[1]
    return _MyContainer2(names, Dict(key(i) => f(i...) for i in indices))
end

@testset "_MyContainer2" begin
    Containers.@container(v[i = 1:3], sin(i), container = _MyContainer2)
    @test v.d isa Dict{Int,Float64}
    @test v.names == [:i]
end
