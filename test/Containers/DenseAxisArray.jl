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

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

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

end  # module

TestContainersDenseAxisArray.runtests()
