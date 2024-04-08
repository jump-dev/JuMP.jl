#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestContainersSparseAxisArray

using JuMP.Containers
using Test

import LinearAlgebra
import OrderedCollections

function _util_sparse_test(d, sum_d, d2, d3, dsqr, d_bads)
    sqr(x) = x^2
    # map
    @test d == @inferred map(identity, d)
    @test dsqr == @inferred map(sqr, d)
    @test d3 == @inferred map(x -> x * 3, d)
    @test d3 == @inferred map(x -> 3 * x, d)
    # reduce
    @test sum_d == @inferred sum(d)
    # broadcasting
    @test dsqr == @inferred d .* d
    @test d2 == @inferred d .+ d
    @test d3 == @inferred d .* 3
    @test d3 == @inferred 3 .* d
    @test d == identity.(d)
    @test dsqr == sqr.(d)
    err = ArgumentError(
        "Cannot broadcast" *
        " Containers.SparseAxisArray with" *
        " another array of different type",
    )
    @test_throws err [1, 2] .+ d
    @test_throws err d .* [1, 2]
    err = ArgumentError(
        "Cannot broadcast" *
        " Containers.SparseAxisArray with " *
        "different indices",
    )
    for d_bad in d_bads
        @test_throws err d_bad .+ d
        @test_throws err d .+ d_bad
    end
    return
end

function test_1_dimensional()
    d = @inferred SparseAxisArray(
        OrderedCollections.OrderedDict((:a,) => 1, (:b,) => 2),
    )
    @test sprint(summary, d) == """
$(SparseAxisArray{Int,1,Tuple{Symbol}}) with 2 entries"""
    @test sprint(show, "text/plain", d) == """
$(SparseAxisArray{Int,1,Tuple{Symbol}}) with 2 entries:
  [a]  =  1
  [b]  =  2"""
    @test d isa SparseAxisArray{Int,1,Tuple{Symbol}}
    d2 = @inferred SparseAxisArray(Dict((:a,) => 2, (:b,) => 4))
    d3 = @inferred SparseAxisArray(Dict((:a,) => 3, (:b,) => 6))
    dsqr = @inferred SparseAxisArray(Dict((:a,) => 1, (:b,) => 4))
    da = @inferred SparseAxisArray(Dict((:b,) => 2))
    dc = @inferred SparseAxisArray(Dict((:a,) => 1, (:b,) => 2, (:c,) => 3))
    _util_sparse_test(d, 3, d2, d3, dsqr, [da, dc])
    # f(::Int) has return type Union{Int, Float64}
    f(x) = x == 1 ? 1 : 1 / x
    fd = f.(d)
    @test fd isa SparseAxisArray{Real,1,Tuple{Symbol}}
    @test fd == SparseAxisArray(Dict((:a,) => 1, (:b,) => 0.5))
    g(x, y) = f(x) + f(y)
    fd = g.(d, 1)
    @test fd isa SparseAxisArray{Real,1,Tuple{Symbol}}
    @test fd == SparseAxisArray(Dict((:a,) => 2, (:b,) => 1.5))
    fd = g.(1, d)
    @test fd isa SparseAxisArray{Real,1,Tuple{Symbol}}
    @test fd == SparseAxisArray(Dict((:a,) => 2, (:b,) => 1.5))
    @test 3 * (d2 / 2) == d3
    @test (d2 / 2) * 3 == d3
    return
end

function test_2_dimensional()
    d = @inferred SparseAxisArray(
        OrderedCollections.OrderedDict((:a, 'u') => 2.0, (:b, 'v') => 0.5),
    )
    @test d isa SparseAxisArray{Float64,2,Tuple{Symbol,Char}}
    @test_throws BoundsError(d, (:a,)) d[:a]
    @test sprint(summary, d) == """
$(SparseAxisArray{Float64,2,Tuple{Symbol,Char}}) with 2 entries"""
    @test sprint(show, "text/plain", d) == """
    $(SparseAxisArray{Float64,2,Tuple{Symbol,Char}}) with 2 entries:
      [a, u]  =  2.0
      [b, v]  =  0.5"""
    d2 = @inferred SparseAxisArray(Dict((:b, 'v') => 1.0, (:a, 'u') => 4.0))
    d3 = @inferred SparseAxisArray(Dict((:a, 'u') => 6.0, (:b, 'v') => 1.5))
    dsqr = @inferred SparseAxisArray(Dict((:a, 'u') => 4.0, (:b, 'v') => 0.25))
    da = @inferred SparseAxisArray(Dict((:b, 'v') => 2.0))
    db = @inferred SparseAxisArray(Dict((:a, 'u') => 1.0, (:b, 'u') => 2.0))
    dc = @inferred SparseAxisArray(
        Dict((:a, 'u') => 1.0, (:b, 'v') => 2.0, (:c, 'w') => 3.0),
    )
    _util_sparse_test(d, 2.5, d2, d3, dsqr, [da, db, dc])
    return
end

function test_3_dimensional()
    d = @inferred SparseAxisArray(
        Dict((:a, 'u', 2) => 2.0, (:b, 'v', 3) => 0.5),
    )
    @test d isa SparseAxisArray{Float64,3,Tuple{Symbol,Char,Int}}
    @test_throws BoundsError(d, (:a,)) d[:a]
    @test_throws BoundsError(d, (:a, 'u')) d[:a, 'u']
    d2 = @inferred SparseAxisArray(
        Dict((:b, 'v', 3) => 1.0, (:a, 'u', 2) => 4.0),
    )
    d3 = @inferred SparseAxisArray(
        Dict((:a, 'u', 2) => 6.0, (:b, 'v', 3) => 1.5),
    )
    dsqr = @inferred SparseAxisArray(
        Dict((:a, 'u', 2) => 4.0, (:b, 'v', 3) => 0.25),
    )
    da = @inferred SparseAxisArray(Dict((:b, 'v', 3) => 2.0))
    db = @inferred SparseAxisArray(
        Dict((:a, 'u', 3) => 1.0, (:b, 'u', 2) => 2.0),
    )
    dc = @inferred SparseAxisArray(
        Dict((:a, 'u', 2) => 1.0, (:b, 'v', 3) => 2.0, (:c, 'w', 4) => 3.0),
    )
    _util_sparse_test(d, 2.5, d2, d3, dsqr, [da, db, dc])
    return
end

function test_empty_array()
    a = Containers.@container([i = 1:3; i > 5], sqrt(i))
    @test a isa SparseAxisArray{Float64,1,Tuple{Int}}
    @test length(a) == 0
    S = [["a"], [:b]]
    b = Containers.@container([i = 1:2, j = S[i]; i > 3], fill(i, j))
    # The compiler doesn't always return the same thing for
    # `@default_eltype`. It gets Tuple{Int,Any} if run from the REPL, but
    # `Tuple` if run from within this testset. Just test for either.
    @test(
        b isa SparseAxisArray{Any,2,Tuple{Int,Any}} ||
        b isa SparseAxisArray{Any,2,Tuple{Any,Any}}
    )
    @test length(b) == 0
    c = Containers.@container(
        [i = 1:0, j = Any[]],
        i,
        container = SparseAxisArray
    )
    @test c isa SparseAxisArray{Int,2,Tuple{Int,Any}}
    @test length(c) == 0
    d = Containers.@container([i = Any[], j = Any[]; isodd(i)], i)
    @test d isa SparseAxisArray{Any,2,Tuple{Any,Any}}
    @test length(d) == 0
    return
end

function test_half_screen_printing()
    d = SparseAxisArray(Dict((i,) => 2 * i for i in 1:100))
    io = IOBuffer()
    show(IOContext(io, :limit => true, :compact => true), "text/plain", d)
    seekstart(io)
    @test occursin("\u22ee", read(io, String))
    return
end

function test_hash()
    a = Containers.@container([i = 1:3; i > 5], sqrt(i))
    @test hash(a) isa UInt
    s = Set{Any}()
    push!(s, a)
    @test length(s) == 1
    return
end

function test_size()
    err = ErrorException(
        "`Base.size` is not implemented for `SparseAxisArray` because " *
        "although it is a subtype of `AbstractArray`, it is conceptually " *
        "closer to a dictionary with `N`-dimensional keys. If you encounter " *
        "this error and you didn't call `size` explicitly, it is because " *
        "you called a method that is unsupported for `SparseAxisArray`s. " *
        "Consult the JuMP documentation for a list of supported operations.",
    )
    x = Containers.@container([i = 1:3, j = i:3], i + j)
    @test_throws err size(x)
    return
end

function test_empty_broadcasting()
    S = Any[]
    x = Containers.@container([S, 1:2], 0, container = SparseAxisArray)
    f(x) = 2x
    y = f.(x)
    @test y isa SparseAxisArray{Any,2,Tuple{Any,Int}}
    @test isempty(y)
    return
end

function test_slicing()
    Containers.@container(x[i = 1:4, j = 1:2; isodd(i + j)], i + j)
    @test x[:, :] == x
    @test x[1, :] == Containers.@container(y[j = 1:2; isodd(1 + j)], 1 + j)
    @test x[:, 1] == Containers.@container(z[i = 1:4; isodd(i + 1)], i + 1)
    @test isempty(x[[1, 3], [1, 3]])
    @test typeof(x[[1, 3], [1, 3]]) == typeof(x)
    @test typeof(x[[1, 3], 1]) == Containers.SparseAxisArray{Int,1,Tuple{Int}}
    @test isempty(x[[1, 3], 1])
    Containers.@container(y[i = 1:4; isodd(i)], i)
    @test y[:] == y
    Containers.@container(y[i = 1:4; isodd(i)], i)
    @test y[[1, 3]] == y
    z = Containers.@container([i = 1:3, j = [:A, :B]; i > 1], (i, j))
    @test z[2, :] == Containers.@container([j = [:A, :B]; true], (2, j))
    @test z[:, :A] == Containers.@container([i = 2:3; true], (i, :A))
    @test z[:, :] == z
    @test z[1:2, :A] == Containers.@container([i = 2:2; true], (i, :A))
    @test z[2, [:A, :B]] == Containers.@container([j = [:A, :B]; true], (2, j))
    @test z[1:2, [:A, :B]] ==
          Containers.@container([i = 2:2, j = [:A, :B]; true], (i, j))
    return
end

function test_slicing_on_set()
    Containers.@container(x[i = 1:4, j = 1:2; isodd(i + j)], i + j)
    err = ArgumentError(
        "Slicing is not supported when calling `setindex!` on a" *
        " SparseAxisArray",
    )
    @test_throws(err, x[:, :] = 1)
    @test_throws(err, x[1, :] = 1)
    @test_throws(err, x[1, 1:2] = 1)
    return
end

function test_ambiguity_broadcast_preserving_zero_d()
    Containers.@container(x[i = 1:2, j = i:3], i + j)
    @test Broadcast.broadcast_preserving_zero_d(*, x, x) == x .* x
    return
end

function test_ambuguity_BroadcastStyleUnknown()
    Containers.@container(x[i = 1:2, j = i:3], i + j)
    style = Base.BroadcastStyle(typeof(x))
    @test_throws ArgumentError Base.BroadcastStyle(style, Broadcast.Unknown())
    return
end

function test_containers_sparseaxisarray_kwarg_indexing()
    Containers.@container(
        x[i = 2:3, j = 1:2],
        i + j,
        container = SparseAxisArray,
    )
    for i in (2, 3, 2:2, 2:3, :), j in (1, 2, 1:2, 1:1, 2:2, :)
        @test x[i = i, j = j] == x[i, j]
        @test_throws ErrorException x[j = j, i = i]
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j = 1, i = 2],
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i = 2, k = 2],
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i = 2, 2],
    )
    Containers.@container(y[i = 2:3, 1:2], i, container = SparseAxisArray,)
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        y[i = 2, 2],
    )
    @test_throws(BoundsError, y[i = 2] = 1)
    @test_throws(BoundsError, y[2] = 1)
    return
end

function test_containers_sparseaxisarray_kwarg_indexing_slicing()
    Containers.@container(
        x[i = 2:3, j = 1:2],
        i + j,
        container = SparseAxisArray,
    )
    y = x[i = 2, j = :]
    @test y[j = 2] == 4
    y = x[i = :, j = 1]
    @test y[i = 3] == 4
    y = x[i = :, j = :]
    @test y[i = 3, j = 1] == 4
    return
end

function test_containers_sparseaxisarray_kwarg_setindex()
    Containers.@container(
        x[i = 2:3, j = 1:2],
        i + j,
        container = SparseAxisArray,
    )
    for i in 2:3, j in 1:2
        @test x[i = i, j = j] == i + j
        x[i = i, j = j] = i + j + 2
        @test x[i = i, j = j] == i + j + 2
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j = 1, i = 2] = 2,
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i = 2, k = 2] = 2,
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i = 2, 2] = 3,
    )
    @test_throws(BoundsError, x[i = 2] = 3)
    return
end

function test_multi_arg_eachindex()
    Containers.@container(x[i = 2:3], i, container = SparseAxisArray)
    Containers.@container(y[i = 2:3], i, container = SparseAxisArray)
    Containers.@container(
        z[i = 2:4, j = 1:2],
        i + j,
        container = SparseAxisArray,
    )
    @test eachindex(x) == keys(x.data)
    @test eachindex(y) == keys(y.data)
    @test eachindex(z) == keys(z.data)
    @test eachindex(x, y) == eachindex(x)
    @test_throws DimensionMismatch eachindex(x, z)
    return
end

function test_sparseaxisarray_order()
    A = [[1, 2, 10], [2, 3, 30]]
    Containers.@container(
        x[i in 1:2, j in A[i]],
        i + j,
        container = SparseAxisArray,
    )
    Containers.@container(x1[j in A[1]], 1 + j, container = SparseAxisArray)
    Containers.@container(x2[j in A[2]], 2 + j, container = SparseAxisArray)
    @test x[1, :] == x1
    @test x[2, :] == x2
    @test LinearAlgebra.dot(x[1, :], 1:3) == 41
    return
end

end  # module
