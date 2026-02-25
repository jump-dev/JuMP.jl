#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

module TestContainersDenseAxisArray

using JuMP.Containers
using Test

function test_undef_constructor()
    A = @inferred DenseAxisArray{Int}(undef, [:a, :b], 1:2)
    @test isassigned(A, :a, 1)  # Because the eltype is Int, isassigned=true.
    @test !isassigned(A, :c, 1)
    @test !isassigned(A, :c, 1, :d)
    A[:a, 1] = 1
    A[:b, 1] = 2
    A[:a, 2] = 3
    A[:b, 2] = 4
    @test A[:a, 1] == 1
    @test A[:b, 1] == 2
    @test A[:a, 2] == 3
    @test A[:b, 2] == 4
    @test isassigned(A, :a, 1)
    @test !isassigned(A, :c, 1)
    @test 10 == @inferred sum(A)
    return
end

function test_undef_constructor_range()
    A = @inferred DenseAxisArray{String}(undef, 1:2)
    @test !isassigned(A, 1)
    @test !isassigned(A, 2)
    @test !isassigned(A, 3)
    A[1] = "abc"
    @test isassigned(A, 1)
    @test !isassigned(A, 2)
    @test !isassigned(A, 3)
    @test !isassigned(A, 2, 2)
    return
end

function test_range_index_set()
    A = @inferred DenseAxisArray([1.0, 2.0], 2:3)
    @test size(A) == (2,)
    @test size(A, 1) == 2
    @test @inferred A[2] == 1.0
    @test A[3] == 2.0
    @test A[2, 1] == 1.0
    @test A[3, 1, 1, 1, 1] == 2.0
    @test isassigned(A, 2)
    @test !isassigned(A, 1)
    @test length.(axes(A)) == (2,)
    @test_throws KeyError A["2"]

    correct_answer = DenseAxisArray([2.0, 3.0], 2:3)
    @test sprint(show, correct_answer) == """
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, 2:3
And data, a 2-element $(Vector{Float64}):
 2.0
 3.0"""

    plus1(x) = x + 1
    @test plus1.(A) == correct_answer
    @test correct_answer == @inferred map(plus1, A)
    @test A .+ 1 == correct_answer
    @test correct_answer == @inferred map(x -> x + 1, A)
    @test 1 .+ A == correct_answer
    @test correct_answer == @inferred map(x -> 1 + x, A)

    correct_answer = DenseAxisArray([2.0, 4.0], 2:3)
    @test 2 * A == correct_answer
    @test correct_answer == @inferred map(x -> 2 * x, A)
    @test A * 2 == correct_answer
    @test correct_answer == @inferred map(x -> x * 2, A)
    @test A / (1 / 2) == correct_answer
    @test correct_answer == @inferred map(x -> x / (1 / 2), A)
    return
end

function test_symbol_index_set()
    A = @inferred DenseAxisArray([1.0, 2.0], [:a, :b])
    @test size(A) == (2,)
    @test size(A, 1) == 2
    @test @inferred A[:a] == 1.0
    @test A[:b] == 2.0
    @test length.(axes(A)) == (2,)
    correct_answer = DenseAxisArray([2.0, 3.0], [:a, :b])
    @test sprint(show, correct_answer) == """
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, $([:a, :b])
And data, a 2-element $(Vector{Float64}):
 2.0
 3.0"""
    plus1(x) = x + 1
    @test plus1.(A) == correct_answer
    @test A .+ 1 == correct_answer
    @test 1 .+ A == correct_answer
    return
end

function test_string_index_set()
    A = @inferred DenseAxisArray([1.0, 2.0], ["a", "b"])
    @test (@inferred A["a"]) == (@inferred A[GenericString("a")]) == 1.0
    @test (@inferred A[["a", "b"]]) ==
          (@inferred A[[GenericString("a"), GenericString("b")]]) ==
          A
    return
end

function test_mixed_range_symbol_index_sets()
    A = @inferred DenseAxisArray([1 2; 3 4], 2:3, [:a, :b])
    @test size(A) == (2, 2)
    @test size(A, 1) == 2
    @test size(A, 2) == 2
    @test_throws BoundsError(A, (2,)) A[2]
    @test length.(axes(A)) == (2, 2)
    @test @inferred A[2, :a] == 1
    @test A[3, :a] == 3
    @test A[2, :b] == 2
    @test A[3, :b] == 4
    @test A[2, :a, 1] == 1
    @test A[2, :a, 1, 1] == 1
    @test A[3, :a, 1, 1, 1] == 3
    @test DenseAxisArray([1, 3], 2:3) == @inferred A[:, :a]
    @test A[2, :] == DenseAxisArray([1, 2], [:a, :b])
    @test sprint(show, A) == """
2-dimensional DenseAxisArray{$Int,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, $([:a, :b])
And data, a 2×2 $(Matrix{Int}):
 1  2
 3  4"""
    return
end

function test_4_dimensional_DenseAxisArray()
    A = DenseAxisArray(zeros(2, 2, 2, 2), 2:3, [:a, :b], -1:0, ["a", "b"])
    @test size(A) == (2, 2, 2, 2)
    @test size(A, 1) == 2
    @test size(A, 2) == 2
    @test size(A, 3) == 2
    @test size(A, 4) == 2
    @test_throws BoundsError(A, (2,)) A[2]
    @test_throws BoundsError(A, (2, :a)) A[2, :a]
    @test_throws BoundsError(A, (2, :a, 0)) A[2, :a, 0]
    A[2, :a, -1, "a"] = 1.0
    f = 0.0
    for I in eachindex(A)
        f += A[I]
    end
    @test f == 1.0
    @test isassigned(A, 2, :a, -1, "a")
    @test A[:, :, -1, "a"] == DenseAxisArray([1.0 0.0; 0.0 0.0], 2:3, [:a, :b])
    @test_throws KeyError A[2, :a, -1, :a]
    @test sprint(summary, A) == """
4-dimensional DenseAxisArray{Float64,4,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, $([:a, :b])
    Dimension 3, -1:0
    Dimension 4, ["a", "b"]
And data, a 2×2×2×2 $(Array{Float64,4})"""
    @test sprint(show, A) == """
4-dimensional DenseAxisArray{Float64,4,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, $([:a, :b])
    Dimension 3, -1:0
    Dimension 4, ["a", "b"]
And data, a 2×2×2×2 $(Array{Float64,4}):
[:, :, -1, "a"] =
 1.0  0.0
 0.0  0.0

[:, :, 0, "a"] =
 0.0  0.0
 0.0  0.0

[:, :, -1, "b"] =
 0.0  0.0
 0.0  0.0

[:, :, 0, "b"] =
 0.0  0.0
 0.0  0.0"""
    return
end

function test_0_dimensional_DenseAxisArray()
    a = Array{Int,0}(undef)
    a[] = 10
    A = DenseAxisArray(a)
    @test size(A) == tuple()
    @test A[] == 10
    A[] = 1
    @test sprint(show, A) == """
0-dimensional DenseAxisArray{$Int,0,...} with index sets:
And data, a 0-dimensional $(Array{Int,0}):
1"""
    return
end

function test_DenseAxisArray_keys()
    A = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, [:a, :b])
    A_keys = collect(keys(A))
    @test A[A_keys[3]] == 6.0
    @test A[A_keys[4]] == 8.0
    @test A_keys[3][1] == 2
    @test A_keys[3][2] == :b
    @test A_keys[4][1] == 3
    @test A_keys[4][2] == :b

    B = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, Set([:a, :b]))
    B_keys = keys(B)
    @test Containers.DenseAxisArrayKey((2, :a)) in B_keys
    @test Containers.DenseAxisArrayKey((2, :b)) in B_keys
    @test Containers.DenseAxisArrayKey((3, :a)) in B_keys
    @test Containers.DenseAxisArrayKey((3, :b)) in B_keys
    return
end

# See https://github.com/jump-dev/JuMP.jl/issues/1988
function test_filter()
    A = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, [:a, :b])
    k = filter(k -> 6 <= A[k] <= 7, keys(A))
    @test k isa Vector{Containers.DenseAxisArrayKey{Tuple{Int,Symbol}}}
    @test k[1] == Containers.DenseAxisArrayKey((3, :a))
    @test k[2] == Containers.DenseAxisArrayKey((2, :b))
    return
end

function test_AxisLookup()
    A = DenseAxisArray([5.0 6.0; 7.0 8.0], [:a, :b], [:a, :b])
    @test A.lookup[1] isa Containers._AxisLookup{Dict{Symbol,Int}}
    @test_throws KeyError A[:c, :a]
    @test_throws KeyError A[1, 1]
    @test_throws KeyError A[:a, :b, 2] == 6.0
    @test isassigned(A, :a, :a)
    @test !isassigned(A, :a, :c)

    @test (@inferred A[:a, :b]) == 6.0
    @test (@inferred A[:a, :b, 1]) == 6.0
    @test (@inferred A[:b, :a]) == 7.0
    @test (@inferred A[[:a, :b], [:a, :b]]) == A
    @test (@inferred A[:a, [:a, :b]]) == DenseAxisArray([5.0, 6.0], [:a, :b])
    @test (@inferred A[[:a, :b], :b]) == DenseAxisArray([6.0, 8.0], [:a, :b])

    B = DenseAxisArray([5.0 6.0; 7.0 8.0], Base.OneTo(2), [:a, :b])
    @test B.lookup[1] isa Containers._AxisLookup{Base.OneTo{Int}}
    @test_throws KeyError B[0, :a]
    @test isassigned(B, 1, :a)
    @test !isassigned(B, 3, :b)

    @test (@inferred B[1, :b]) == 6.0
    @test (@inferred B[2, :a]) == 7.0
    @test (@inferred B[1:2, [:a, :b]]) == B
    @test (@inferred B[1, [:a, :b]]) == DenseAxisArray([5.0, 6.0], [:a, :b])
    @test (@inferred B[1:2, :b]) == DenseAxisArray([6.0, 8.0], 1:2)

    C = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, [:a, :b])
    @test C.lookup[1] isa Containers._AxisLookup{Tuple{Int,Int}}
    @test_throws KeyError C[0, :a]
    @test isassigned(C, 2, :a)
    @test !isassigned(C, 4, :b)
    @test (@inferred C[2, :b]) == 6.0
    @test (@inferred C[3, :a]) == 7.0
    @test (@inferred C[2:3, [:a, :b]]) == C
    @test (@inferred C[2, [:a, :b]]) == DenseAxisArray([5.0, 6.0], [:a, :b])
    @test (@inferred C[2:3, :b]) == DenseAxisArray([6.0, 8.0], 2:3)

    D = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, ["a", "b"])
    @test (@inferred D[2, GenericString("b")]) == 6.0
    @test (@inferred D[2, [GenericString("a"), GenericString("b")]]) ==
          DenseAxisArray([5.0, 6.0], ["a", "b"])
    return
end

function test_BitArray()
    x = DenseAxisArray([0 1; 1 0], [:a, :b], 1:2)
    y = Bool.(x)
    @test y isa DenseAxisArray
    @test x == y
    return
end

function test_Broadcast()
    foo(x, y) = x + y
    foo_b(x, y) = foo.(x, y)
    bar(x, y) = (foo.(x, y) .+ x) .^ 2
    a = [5.0 6.0; 7.0 8.0]
    A = DenseAxisArray(a, [:a, :b], [:a, :b])
    b = a .+ 1
    B = A .+ 1
    @test B == DenseAxisArray(b, [:a, :b], [:a, :b])
    C = @inferred foo_b(A, B)
    @test C == DenseAxisArray(foo_b(a, b), [:a, :b], [:a, :b])
    D = @inferred bar(A, B)
    @test D == DenseAxisArray(bar(a, b), [:a, :b], [:a, :b])
    return
end

function test_Broadcast_errors()
    a = [5.0 6.0; 7.0 8.0]
    A = DenseAxisArray(a, [:a, :b], [:a, :b])
    B = DenseAxisArray(a, [:b, :a], [:a, :b])
    @test_throws ErrorException A .+ B
    b = [5.0 6.0; 7.0 8.0; 9.0 10.0]
    @test_throws DimensionMismatch A .+ b
    return
end

function test_DenseAxisArray_with_Base_OneTo()
    A = @inferred DenseAxisArray([1, 3, 2], Base.OneTo(3))
    B = @inferred map(x -> x^2, A)
    @test B isa DenseAxisArray
    @test B.data == [1, 9, 4]
    @test B.axes == (Base.OneTo(3),)
    C = @inferred DenseAxisArray([1 3; 2 4], Base.OneTo(2), Base.OneTo(2))
    D = @inferred map(x -> x - 1, C)
    @test D isa DenseAxisArray
    @test D.data == [0 2; 1 3]
    @test D.axes == (Base.OneTo(2), Base.OneTo(2))
    return
end

function test_Array()
    A = DenseAxisArray([1, 3, 2], Base.OneTo(3))
    B = @inferred Array(A)
    @test B isa Vector{Int}
    @test B == [1, 3, 2]
    # Test mutating B doesn't mutate A
    B[2] = 4
    @test A[2] == 3
    C = @inferred Array{Float64}(A)
    @test C isa Vector{Float64}
    @test C == [1.0, 3.0, 2.0]
    return
end

function test_hash()
    a = [5.0 6.0; 7.0 8.0]
    A = DenseAxisArray(a, [:a, :b], Base.OneTo(2))
    @test hash(A) isa UInt
    s = Set{Any}()
    push!(s, A)
    @test length(s) == 1
    return
end

function test_non_AbstractArray_axes()
    x = [1.0, 2.0, 3.0]
    d = Dict(:a => "a", :b => "b", :c => "c")
    X = DenseAxisArray(x, d)
    for (k, v) in d
        @test X[(k, v)] in x
        @test X[(k, v)] == X[k=>v]
    end
    @test_throws KeyError X[(:a, "b")]
    @test isassigned(X, (:a, "a"))
    @test !isassigned(X, (:a, "b"))
    @test length(X[[(:a, "a"), (:c, "c")]]) == 2
    return
end

function test_non_AbstractArray_matrix()
    x = [1.0 2.0 3.0; 1.0 2.0 3.0; 1.0 2.0 3.0]
    d = Dict(:a => "a", :b => "b", :c => "c")
    X = DenseAxisArray(x, d, d)
    for (k, v) in d
        @test X[(k, v), (k, v)] in x
        @test X[(k, v), (k, v)] == X[k=>v, k=>v]
    end
    @test_throws BoundsError X[(:a, "b")]
    @test_throws KeyError X[(:a, "b"), (:a, "a")]
    @test_throws KeyError X[(:a, "a"), (:a, "b")]
    @test isassigned(X, (:a, "a"), (:a, "a"))
    @test !isassigned(X, (:a, "b"))
    @test isassigned(X, (:a, "a"), (:b, "b"))
    y = Array(X[:, (:a, "a")])
    @test all(y .== y[1])
    return
end

function test_Singular_axis()
    x = @test_logs (:warn,) DenseAxisArray([1.1 2.2], 1, 1:2)
    @test x[1, 2] == 2.2
    y = @test_logs DenseAxisArray([1.1 2.2], [1], 1:2)
    @test y[1, 2] == 2.2
    return
end

function test_CartesianIndex_error()
    S = CartesianIndex.([2, 4])
    err = ErrorException(
        "Unsupported index type `CartesianIndex` in axis: $S. Cartesian " *
        "indices are restricted for indexing into and iterating over " *
        "multidimensional arrays.",
    )
    @test_throws(err, DenseAxisArray([1.1, 2.2], S))
    return
end

function test_Matrix_indices()
    sources = ["A", "B", "C"]
    sinks = ["D", "E"]
    S = [(source, sink) for source in sources, sink in sinks]
    x = DenseAxisArray(1:6, S)
    @test size(x) == (6,)
    return
end

function test_DenseAxisArray_show_nd()
    S = zeros(Int, 2, 2, 3, 3, 3)
    for i in eachindex(S)
        S[i] = i
    end
    x = DenseAxisArray(S, 1:2, 1:2, 1:3, 1:3, 1:3)
    str = sprint((io, x) -> Base.show_nd(io, x, Base.print_matrix, true), x)
    @test occursin("[:, :, 1, 2, 3] =\n 85  87\n 86  88\n", str)
    str_limit = sprint(x) do io, x
        return Base.show_nd(
            IOContext(io, :limit => true),
            x,
            Base.print_matrix,
            true,
        )
    end
    @test occursin("[:, :, 1, 2, 3] =\n 85  87\n 86  88\n", str_limit)
    return
end

function test_DenseAxisArray_show_nd_limit()
    S = zeros(Int, 2, 2, 3, 3, 20)
    for i in eachindex(S)
        S[i] = i
    end
    x = DenseAxisArray(S, 1:2, 1:2, 1:3, 1:3, 1:20)
    str = sprint((io, x) -> Base.show_nd(io, x, Base.print_matrix, true), x)
    @test occursin("[:, :, 1, 1, 3]", str)
    @test occursin("[:, :, 1, 1, 4]", str)
    @test occursin("[:, :, 1, 1, 17]", str)
    @test occursin("[:, :, 1, 1, 18]", str)
    str_limit = sprint(x) do io, x
        return Base.show_nd(
            IOContext(io, :limit => true),
            x,
            Base.print_matrix,
            true,
        )
    end
    @test occursin("[:, :, 1, 1, 3]", str_limit)
    @test !occursin("[:, :, 1, 1, 4]", str_limit)
    @test !occursin("[:, :, 1, 1, 17]", str_limit)
    @test occursin("[:, :, 1, 1, 18]", str_limit)
    return
end

function test_DenseAxisArray_show_nd_empty()
    x = DenseAxisArray(Int[], 1:0)
    str = sprint((io, x) -> Base.show_nd(io, x, Base.print_matrix, true), x)
    @test isempty(str)
    return
end

function test_DenseAxisArray_vector_keys()
    paths = [[1, 2, 15, 3, 20], [1, 9, 16, 20], [1, 2, 20]]
    x = DenseAxisArray(1:3, paths)
    for i in 1:3
        @test x[paths[i]] == i
        @test isassigned(x, paths[i]) == true
    end
    @test isassigned(x, Int[]) == false
    return
end

function test_containers_denseaxisarray_setindex_vector()
    A = Containers.DenseAxisArray(zeros(3), 1:3)
    A[2:3] .= 1.0
    @test A.data == [0.0, 1.0, 1.0]
    A = Containers.DenseAxisArray(zeros(3), 1:3)
    A[[2, 3]] .= 1.0
    @test A.data == [0.0, 1.0, 1.0]
    A = Containers.DenseAxisArray(zeros(3), 1:3)
    A[[1, 3]] .= 1.0
    @test A.data == [1.0, 0.0, 1.0]
    A = Containers.DenseAxisArray(zeros(3), 1:3)
    A[[2]] .= 1.0
    @test A.data == [0.0, 1.0, 0.0]
    A[2:3] = Containers.DenseAxisArray([2.0, 3.0], 2:3)
    @test A.data == [0.0, 2.0, 3.0]
    A = Containers.DenseAxisArray(zeros(3), 1:3)
    A[:] .= 1.0
    @test A.data == [1.0, 1.0, 1.0]
    return
end

function test_containers_denseaxisarray_setindex_matrix()
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    A[:, [:a, :b]] .= 1.0
    @test A.data == [1.0 1.0 0.0; 1.0 1.0 0.0; 1.0 1.0 0.0]
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    A[2:3, [:a, :b]] .= 1.0
    @test A.data == [0.0 0.0 0.0; 1.0 1.0 0.0; 1.0 1.0 0.0]
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    A[3:3, [:a, :b]] .= 1.0
    @test A.data == [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 1.0 0.0]
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    A[[1, 3], [:a, :b]] .= 1.0
    @test A.data == [1.0 1.0 0.0; 0.0 0.0 0.0; 1.0 1.0 0.0]
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    A[[1, 3], [:a, :c]] .= 1.0
    @test A.data == [1.0 0.0 1.0; 0.0 0.0 0.0; 1.0 0.0 1.0]
    return
end

function test_containers_denseaxisarray_view()
    A = Containers.DenseAxisArray(zeros(3, 3), 1:3, [:a, :b, :c])
    B = view(A, :, [:a, :b])
    @test_throws KeyError view(A, :, [:d])
    @test size(B) == (3, 2)
    @test B[1, :a] == A[1, :a]
    @test B[3, :a] == A[3, :a]
    @test_throws KeyError B[3, :c]
    @test sprint(show, B) == sprint(show, B.data)
    @test sprint(Base.print_array, B) == sprint(show, B.data)
    @test sprint(Base.summary, B) ==
          "view(::DenseAxisArray, 1:3, [:a, :b]), over"
    return
end

function test_containers_denseaxisarray_jump_3151()
    D = Containers.DenseAxisArray(zeros(3), [:a, :b, :c])
    E = Containers.DenseAxisArray(ones(3), [:a, :b, :c])
    I = [:a, :b]
    D[I] = E[I]
    @test D.data == [1.0, 1.0, 0.0]
    D = Containers.DenseAxisArray(zeros(3), [:a, :b, :c])
    I = [:b, :c]
    D[I] = E[I]
    @test D.data == [0.0, 1.0, 1.0]
    D = Containers.DenseAxisArray(zeros(3), [:a, :b, :c])
    I = [:a, :c]
    D[I] = E[I]
    @test D.data == [1.0, 0.0, 1.0]
    return
end

function test_containers_denseaxisarray_view_drops_dimension()
    x = Containers.@container([i = 4:6, j = [:A, :B]], (i, j))
    y = @view x[5:6, :A]
    @test axes(y) == (5:6,)
    @test (@inferred y[5]) == (5, :A)
    @test (@inferred y[6]) == (6, :A)
    x = Containers.@container([i = 4:6, j = [:A, :B], k = 1:3], (i, j, k))
    y = @view x[5:6, :A, 1:2]
    @test axes(y) == (5:6, 1:2)
    @test (@inferred y[5, 1]) == (5, :A, 1)
    @test (@inferred y[6, 2]) == (6, :A, 2)
    return
end

function test_containers_denseaxisarray_view_operations()
    c = Containers.@container([i = 1:4, j = 2:3], i + 2 * j)
    d = view(c, 2:3, :)
    @test sum(c) == 60
    @test sum(d) == 30
    d .= 1
    @test sum(d) == 4
    @test sum(c) == 34
    return
end

function test_containers_denseaxisarray_view_addition()
    c = Containers.@container([i = 1:4, j = 2:3], i + 2 * j)
    d = view(c, 2:3, :)
    @test_throws MethodError d + d
    return
end

function test_containers_denseaxisarray_view_colon()
    c = Containers.@container([i = 1:4, j = 2:3], i + 2 * j)
    d = view(c, 2:3, :)
    @test d[:, 2] == Containers.@container([i = 2:3], i + 2 * 2)
    return
end

function test_containers_denseaxisarray_setindex_invalid()
    c = Containers.@container([i = 1:4, j = 2:3], 0)
    d = Containers.@container([i = 1:4, j = 2:3], i + 2 * j)
    setindex!(c, d, 1:4, 2:3)
    @test c == d
    c .= 0
    setindex!(c, d, 1:4, 2:2)
    @test c == Containers.@container([i = 1:4, j = 2:3], (4 + i) * (j == 2))
    d = Containers.@container([i = 5:6, j = 2:3], i + 2 * j)
    @test_throws KeyError setindex!(c, d, 1:4, 2:3)
    return
end

function test_containers_denseaxisarray_setindex_keys()
    c = Containers.@container([i = 1:4, j = 2:3], 0)
    for (i, k) in enumerate(keys(c))
        c[k] = c[k] + i
    end
    @test c == Containers.@container([i = 1:4, j = 2:3], 4 * (j - 2) + i)
    for (i, k) in enumerate(keys(c))
        c[k] = c[k] + i
    end
    @test c == Containers.@container([i = 1:4, j = 2:3], 2 * (4 * (j - 2) + i))
    return
end

function test_ambiguity_isassigned()
    x = DenseAxisArray([:a, :b, :c], 2:4)
    @test !isassigned(x, 1)
    @test isassigned(x, 2)
    @test isassigned(x, CartesianIndex(1))
    @test isassigned(x, CartesianIndex(2, 1))
    @test !isassigned(x, CartesianIndex(2, 1), 1)
    return
end

function test_containers_denseaxisarray_view_axes_n()
    x = Containers.@container([i = 4:6, j = [:A, :B]], (i, j))
    y = @view x[5:6, :A]
    @test sprint(show, MIME("text/plain"), y) isa String
    @test axes(y) == (5:6,)
    @test axes(y, 1) == 5:6
    @test axes(y, 2) == Base.OneTo(1)
    return
end

function test_containers_denseaxisarray_vector_any()
    key_1 = Any[Any["a", 1], "b"]
    key_2 = Any[Any["a", 2], "c"]
    K = Any[key_1, key_2]
    Containers.@container(x[k=K], k[1][2])
    @test axes(x) == (K,)
    @test x[key_1] == 1
    @test x[key_2] == 2
    @test x[Any[key_2, key_1]] ==
          Containers.DenseAxisArray([2, 1], Any[key_2, key_1])
    return
end

function test_containers_denseaxisarray_ambiguous_slice()
    K = Any[Any["a"], Any["b"], Any[Any["a"], Any["b"]]]
    Containers.@container(x[k=K], length(k))
    @test axes(x) == (K,)
    @test x[Any["a"]] == 1
    @test x[Any["b"]] == 1
    new_key = Any[Any["b"], Any["a"]]
    @test x[new_key] == Containers.DenseAxisArray([1, 1], new_key)
    key = reverse(new_key)
    @test_throws(
        ErrorException(
            "ambiguous use of getindex with key $key. We cannot tell if " *
            "you meant to return the single element corresponding to the " *
            "key, or a slice for each element in the key.",
        ),
        x[key],
    )
    return
end

function test_containers_denseaxisarray_kwarg_indexing()
    Containers.@container(x[i=2:3, j=1:2], i + j)
    for i in (2, 3, 2:2, 2:3, :), j in (1, 2, 1:2, 1:1, 2:2, :)
        @test x[i=i, j=j] == x[i, j]
        @test_throws ErrorException x[j=j, i=i]
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j=1, i=2],
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i=2, k=2],
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i=2, 2],
    )
    Containers.@container(y[i=2:3, 1:2], i)
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        y[i=2, 2],
    )
    return
end

function test_containers_denseaxisarray_kwarg_setindex()
    Containers.@container(x[i=2:3, j=1:2], i + j)
    for i in 2:3, j in 1:2
        @test x[i=i, j=j] == i + j
        x[i=i, j=j] = i + j + 2
        @test x[i=i, j=j] == i + j + 2
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j=1, i=2] = 2,
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i=2, k=2] = 2,
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i=2, 2] = 3,
    )
    return
end

function test_containers_denseaxisarray_kwarg_indexing_slicing()
    Containers.@container(x[i=2:3, j=1:2], i + j)
    y = x[i=2, j = :]
    @test y[j=2] == 4
    y = x[i = :, j=1]
    @test y[i=3] == 4
    y = x[i = :, j = :]
    @test y[i=3, j=1] == 4
    return
end

function test_containers_denseaxisarrayview_kwarg_indexing()
    Containers.@container(a[i=2:3, j=1:2], i + j)
    x = view(a, :, :)
    for i in (2, 3, 2:2, 2:3, :), j in (1, 2, 1:2, 1:1, 2:2, :)
        @test x[i=i, j=j] == x[i, j]
        @test_throws ErrorException x[j=j, i=i]
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j=1, i=2],
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i=2, k=2],
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i=2, 2],
    )
    return
end

function test_containers_denseaxisarrayview_kwarg_indexing_drop_dim()
    Containers.@container(a[i=2:3, j=1:2], i + j)
    x = view(a, 2, 1:2)
    for j in (1, 2, 1:2, 1:1, 2:2, :)
        @test x[j=j] == x[j]
    end
    @test_throws(
        ErrorException(
            "Invalid index i in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i=2],
    )
    return
end

function test_containers_denseaxisarrayview_kwarg_indexing_slicing()
    Containers.@container(a[i=2:3, j=1:2], i + j)
    x = view(a, :, :)
    y = x[i=2, j = :]
    @test y[j=2] == 4
    y = x[i = :, j=1]
    @test y[i=3] == 4
    y = x[i = :, j = :]
    @test y[i=3, j=1] == 4
    return
end

function test_containers_denseaxisarrayview_kwarg_setindex()
    Containers.@container(a[i=2:3, j=1:2], i + j)
    x = view(a, :, :)
    for i in 2:3, j in 1:2
        @test x[i=i, j=j] == i + j
        x[i=i, j=j] = i + j + 2
        @test x[i=i, j=j] == i + j + 2
    end
    @test_throws(
        ErrorException(
            "Invalid index j in position 1. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[j=1, i=2] = 2,
    )
    @test_throws(
        ErrorException(
            "Invalid index k in position 2. When using keyword indexing, the " *
            "indices must match the exact name and order used when creating " *
            "the container.",
        ),
        x[i=2, k=2] = 2,
    )
    @test_throws(
        ErrorException(
            "Cannot index with mix of positional and keyword arguments",
        ),
        x[i=2, 2] = 3,
    )
    return
end

function test_sum_dims()
    Containers.@container(x[i=1:2, j=1:2], i + j, container = DenseAxisArray)
    @test_throws(
        ErrorException(
            "`sum(x::DenseAxisArray; dims)` is not supported. Convert the array " *
            "to an `Array` using `sum(Array(x); dims=2)`, or use an explicit " *
            "for-loop summation instead.",
        ),
        sum(x; dims = 2),
    )
    return
end

function test_multi_arg_eachindex()
    Containers.@container(x[i=2:3], i)
    Containers.@container(y[i=2:3], i)
    Containers.@container(z[i=2:4, j=1:2], i + j)
    @test eachindex(x) == CartesianIndices((2,))
    @test eachindex(y) == CartesianIndices((2,))
    @test eachindex(z) == CartesianIndices((3, 2))
    @test eachindex(x, y) == CartesianIndices((2,))
    @test_throws DimensionMismatch eachindex(x, z)
    return
end

function test_LinearIndices()
    Containers.@container(x[i in 2:3], i)
    @test_throws(
        ErrorException("DenseAxisArray does not support this operation."),
        LinearIndices(x),
    )
    return
end

function test_CartesianIndices()
    Containers.@container(x[i in 2:3], i)
    @test CartesianIndices(x) == CartesianIndices((2,))
    return
end

function test_show_nd()
    Containers.@container(
        x[a in 2:3, b in 2:3, c in 2:14, d in 2:14],
        (a, b, c, d),
    )
    s = sprint(io -> show(IOContext(io, :limit => true), x))
    limit_indices = [2, 3, 4, 12, 13, 14]
    for c in 2:14, d in 2:14
        is_visible = (c in limit_indices && d in limit_indices)
        @test occursin("[:, :, $c, $d]", s) == is_visible
        for a in 2:3, b in 2:3
            @test occursin("($a, $b, $c, $d)", s) == is_visible
        end
    end
    s = sprint(io -> show(IOContext(io, :limit => false), x))
    limit_indices = [2, 3, 4, 12, 13, 14]
    for c in 2:14, d in 2:14
        @test occursin("[:, :, $c, $d]", s)
        for a in 2:3, b in 2:3
            @test occursin("($a, $b, $c, $d)", s)
        end
    end
    return
end

function test_view_DenseAxisArray()
    Containers.@container(x[a in 2:3], a)
    @test_throws KeyError view(x, 3:4)
    @test_throws KeyError view(x, 4)
    return
end

function test_promote_shape()
    Containers.@container(x[a in 2:3], a)
    Containers.@container(y[a in 2:3], 2 * a)
    Containers.@container(z[a in 2:3, b in ["a", "b"]], (a, b))
    @test x + x == y
    @test promote_shape(x, x) == (2:3,)
    @test promote_shape(x, y) == (2:3,)
    @test promote_shape(y, x) == (2:3,)
    @test_throws DimensionMismatch promote_shape(x, z)
    return
end

function test_container_Base_OneTo_Integer()
    Containers.@container(x[i=0x01:0x08], 2 * i)
    @test x[0x01:0x02] == Containers.@container([i = 0x01:0x02], 2 * i)
    return
end

function test_conntainer_AbstractUnitRange_Integer()
    Containers.@container(A[i in 0x01:0x03], i)
    for i in 0x01:0x03
        @test A[i] == i
    end
    Containers.@container(B[i in 0x02:0x03], i)
    @test_throws KeyError B[0x01]
    for i in 0x02:0x03
        @test B[i] == i
    end
    return
end

struct Int4053 <: Integer
    x::Int
end

Int4053(f::Int4053) = f

Base.:<(a::Int4053, b::Int4053) = a.x < b.x

Base.:<=(a::Int4053, b::Int4053) = a.x <= b.x

Base.:+(a::Int4053, b::Int4053) = Int4053(a.x + b.x)

Base.:-(a::Int4053, b::Int4053) = Int4053(a.x - b.x)

Base.:(==)(a::Int4053, b::Int4053) = a.x == b.x

Base.hash(a::Int4053, h::UInt) = hash((Int4053, a.x), h)

Base.:(==)(a::Int4053, b::Integer) = false

Base.:(==)(a::Integer, b::Int4053) = false

Base.promote_rule(T::Type{<:Integer}, ::Type{Int4053}) = T

Base.Int(x::Int4053) = x.x

function test_issue_4053()
    Containers.@container(A[i in Int4053(1):Int4053(3)], i.x)
    @test_throws KeyError A[1]
    @test !isassigned(A, 1)
    @test_throws KeyError A[0x01]
    @test !isassigned(A, 0x01)
    @test A[Int4053(1)] === 1
    @test isassigned(A, Int4053(1))
    Containers.@container(B[i in 2:4], i)
    @test B[2] === 2
    @test B[0x02] === 2
    @test_throws KeyError B[Int4053(2)]
    @test !isassigned(B, Int4053(2))
    return
end

const _REDUCING_ERROR = VERSION < v"1.11" ? MethodError : ArgumentError

function test_sum_init()
    x = Containers.@container([i in Int[]], i)
    @test_throws _REDUCING_ERROR sum(x)
    @test sum(x; init = 1) == 1
    y = Containers.@container([i in BigInt[]], i)
    @test_throws _REDUCING_ERROR sum(y)
    y_2 = sum(y; init = 0)
    @test y_2 === 0
    return
end

function test_sum_init_any()
    x = Containers.@container([i in Any[]], i)
    @test_throws _REDUCING_ERROR sum(x)
    return
end

end  # module
