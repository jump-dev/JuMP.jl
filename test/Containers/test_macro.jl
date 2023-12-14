#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestContainersMacro

using JuMP.Containers
using Test

function test_Array()
    Containers.@container(x[i = 1:3], i^2)
    @test x isa Vector{Int}
    x = Containers.@container([i = 1:3, j = 1:3], i^2)
    @test x isa Matrix{Int}
    # alternative syntax
    Containers.@container(x[i in 1:3], i^2)
    @test x isa Vector{Int}
    Containers.@container(x[i ∈ 1:3], i^2)
    @test x isa Vector{Int}
    return
end

function test_forced_array()
    set = 1:2
    x = Containers.@container([i = set, j = 1:2], i + j, container = Array)
    @test x == [2 3; 3 4]
    set_2 = [1, 2]
    y = Containers.@container([i = set_2, j = 1:2], i + j, container = Array)
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
        Containers.@container([i = 1:2:4, j = 1:2], i + j, container = Array)
    )
    return
end

function test_DenseAxisArray()
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
    return
end

function test_SparseAxisArray()
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
    @test(x isa SparseAxisArray{Int,2,<:T} || x isa SparseAxisArray{Any,2,<:T})
    # This one is better, we can infer the value type, but the keys are
    # difficult to infer, and depend on the Julia version you are running.
    Containers.@container(x[i = 1:2, j = 1:2; false], i)
    @test x isa SparseAxisArray{Int,2,<:T}
    Containers.@container(x[i = 1:0, j = 2:1], i, container = SparseAxisArray)
    @test x isa SparseAxisArray{Int,2,Tuple{Int,Int}}
    Containers.@container(x[i = 1:0, j = 1:0], i, container = SparseAxisArray)
    @test x isa SparseAxisArray{Int,2,Tuple{Int,Int}}
    return
end

function test_duplicate_indices()
    expr = :(Containers.@container(x[i = 1:2, i = 1:2], i + i))
    @test_throws(LoadError, try
        @eval $expr
    catch err
        throw(err)
    end,)
    return
end

function test_double_filter_typed_vcat()
    expr = :(Containers.@container(x[i = 1:2; isodd(i); iseven(i + 1)], i + i))
    @test_throws(LoadError, try
        @eval $expr
    catch err
        throw(err)
    end,)
    return
end

function test_double_filter_vect()
    expr = :(Containers.@container(
        x[i = 1:2, 1:2; isodd(i); iseven(i + 1)],
        i + i,
    ))
    @test_throws(LoadError, try
        @eval $expr
    catch err
        throw(err)
    end,)
    return
end

function test_Dict()
    Containers.@container(v[i = 1:3], sin(i), container = Dict)
    @test v isa Dict{Int,Float64}
    @test length(v) == 3
    @test v[2] ≈ sin(2)
    Containers.@container(w[i = 1:3, j = 1:3], i + j, container = Dict)
    @test w isa Dict{Tuple{Int,Int},Int}
    @test length(w) == 9
    @test w[2, 3] == 5
    Containers.@container(x[i = 1:3, j = [:a, :b]], (j, i), container = Dict)
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
    return
end

function test_invalid_container()
    err = ErrorException(
        "Unable to build a container with the provided type $(Int). " *
        "Implement `Containers.container(::Function, indices, ::Type{$Int})`.",
    )
    @test_throws err Containers.@container(
        x[i = 1:2, j = 1:2],
        i + j,
        container = Int
    )
    return
end

function test_compound_indexing_expressions()
    Containers.@container(x[(i, j) in [(1, 1), (2, 2)], k in i:3], i + j + k,)
    @test x isa Containers.SparseAxisArray
    @test length(x) == 5
    @test x[(2, 2), 3] == 7
    return
end

struct _MyContainer end

function Containers.container(f::Function, indices, ::Type{_MyContainer})
    key(i::Tuple) = i
    key(i::Tuple{T}) where {T} = i[1]
    return Dict(key(i) => f(i...) for i in indices)
end

function test__MyContainer()
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
    return
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

function test__MyContainer2()
    Containers.@container(v[i = 1:3], sin(i), container = _MyContainer2)
    @test v.d isa Dict{Int,Float64}
    @test v.names == [:i]
    return
end

function test_parse_macro_arguments()
    args, kwargs = Containers.parse_macro_arguments(error, ())
    @test args == Any[]
    @test isempty(kwargs)
    return
end

function test_add_additional_args()
    call = :(f(1; a = 2))
    kwargs = Dict{Symbol,Any}()
    @test Containers.add_additional_args(call, [:(foo)], kwargs) === nothing
    @test call == :(f(1, $(Expr(:escape, :foo)); a = 2))
    call = :(f(1))
    Containers.add_additional_args(call, [2, 3], kwargs)
    @test call == :(f(1, $(esc(2)), $(esc(3))))
    call = :(f.(1))
    Containers.add_additional_args(call, [2, 3], kwargs)
    @test call == :(f.(1, $(esc(2)), $(esc(3))))
    call = :(f(1; a = 4))
    Containers.add_additional_args(call, [2, 3], kwargs)
    @test call == :(f(1, $(esc(2)), $(esc(3)); a = 4))
    call = :(f.(1; a = 4))
    Containers.add_additional_args(call, [2, 3], kwargs)
    @test call == :(f.(1, $(esc(2)), $(esc(3)); a = 4))
    call = :(f.(1, a = 4))
    kwargs = Dict{Symbol,Any}(:b => 4, :c => false)
    Containers.add_additional_args(call, Any[2], kwargs; kwarg_exclude = [:b])
    @test call == Expr(
        :.,
        :f,
        Expr(:tuple, 1, esc(2), Expr(:kw, :a, 4), esc(Expr(:kw, :c, false))),
    )
    return
end

end  # module
