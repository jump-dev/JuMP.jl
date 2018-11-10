#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuMP
using Compat # For undef
using Compat.Test

macro dummycontainer(expr, requestedtype)
    name = gensym()
    refcall, indexvars, indexsets, condition = JuMP.buildrefsets(expr, name)
    if condition == :()
        return JuMP.generatecontainer(Bool, indexvars, indexsets, requestedtype)[1]
    else
        if requestedtype != :Auto && requestedtype != :Dict
            return :(error(""))
        end
        return JuMP.generatecontainer(Bool, indexvars, indexsets, :Dict)[1]
    end
end

function containermatches(c1::AbstractArray,c2::AbstractArray)
    return typeof(c1) == typeof(c2) && size(c1) == size(c2)
end

function containermatches(c1::JuMPArray,c2::JuMPArray)
    return typeof(c1) == typeof(c2) && Compat.axes(c1) == Compat.axes(c2)
end

containermatches(c1::Dict, c2::Dict) = (eltype(c1) == eltype(c2))
containermatches(c1, c2) = false

@testset "Container syntax" begin
    @test containermatches(@dummycontainer([i=1:10], Auto), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=1:10], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=1:10], JuMPArray), JuMPArray(Vector{Bool}(undef,10), 1:10))
    @test containermatches(@dummycontainer([i=1:10], Dict), Dict{Any,Bool}())

    @test containermatches(@dummycontainer([i=1:10,1:2], Auto), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=1:10,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=1:10,n=1:2], JuMPArray), JuMPArray(Matrix{Bool}(undef,10,2), 1:10, 1:2))
    @test containermatches(@dummycontainer([i=1:10,1:2], Dict), Dict{Any,Bool}())

    @test containermatches(@dummycontainer([i=1:10,n=2:3], Auto), JuMPArray(Matrix{Bool}(undef,10,2), 1:10, 2:3))
    @test_throws ErrorException @dummycontainer([i=1:10,2:3], Array)
    @test containermatches(@dummycontainer([i=1:10,n=2:3], JuMPArray), JuMPArray(Matrix{Bool}(undef,10,2), 1:10, 2:3))
    @test containermatches(@dummycontainer([i=1:10,n=2:3], Dict), Dict{Any,Bool}())


    S = Base.OneTo(10)
    @test containermatches(@dummycontainer([i=S], Auto), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], JuMPArray), JuMPArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], Dict), Dict{Any,Bool}())

    @test containermatches(@dummycontainer([i=S,1:2], Auto), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,n=1:2], JuMPArray), JuMPArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], Dict), Dict{Any,Bool}())

    S = 1:10
    # Not type stable to return an Array by default even when S is one-based interval
    @test containermatches(@dummycontainer([i=S], Auto), JuMPArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], Array), Vector{Bool}(undef,10))
    @test containermatches(@dummycontainer([i=S], JuMPArray), JuMPArray(Vector{Bool}(undef,10), S))
    @test containermatches(@dummycontainer([i=S], Dict), Dict{Any,Bool}())

    @test containermatches(@dummycontainer([i=S,n=1:2], Auto), JuMPArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], Array), Matrix{Bool}(undef,10,2))
    @test containermatches(@dummycontainer([i=S,n=1:2], JuMPArray), JuMPArray(Matrix{Bool}(undef,10,2), S, 1:2))
    @test containermatches(@dummycontainer([i=S,1:2], Dict), Dict{Any,Bool}())

    # TODO: test case where S is index set not supported by JuMPArrays (does this exist?)

    # Conditions
    @test containermatches(@dummycontainer([i=1:10; iseven(i)], Auto), Dict{Any,Bool}())
    @test_throws ErrorException @dummycontainer([i=1:10; iseven(i)], Array)
    @test_throws ErrorException @dummycontainer([i=1:10; iseven(i)], JuMPArray)
    @test containermatches(@dummycontainer([i=1:10; iseven(i)], Dict), Dict{Any,Bool}())

    # Dependent axes
    @test containermatches(@dummycontainer([i=1:10, j=1:i], Auto), Dict{Any,Bool}())
    @test_throws ErrorException @dummycontainer([i=1:10, j=1:i], Array)
    @test_throws ErrorException @dummycontainer([i=1:10, j=1:i], JuMPArray)
    @test containermatches(@dummycontainer([i=1:10, j=1:i], Dict), Dict{Any,Bool}())

end

@testset "JuMPArray" begin
    @testset "Range index set" begin
        A = @inferred JuMPArray([1.0,2.0], 2:3)
        if VERSION >= v"0.7-"
            @test size(A) == (2,)
            @test size(A, 1) == 2
        end
        @test @inferred A[2] == 1.0
        @test A[3] == 2.0
        @test A[2,1] == 1.0
        @test A[3,1,1,1,1] == 2.0
        @test isassigned(A, 2)
        @test !isassigned(A, 1)
        @test length.(Compat.axes(A)) == (2,)
        plus1(x) = x + 1
        B = plus1.(A)
        @test B[2] == 2.0
        @test B[3] == 3.0
        @test sprint(show, B) == """
1-dimensional JuMPArray{Float64,1,...} with index sets:
    Dimension 1, 2:3
And data, a 2-element Array{Float64,1}:
 2.0
 3.0"""
    end

    @testset "Symbol index set" begin
        A = @inferred JuMPArray([1.0,2.0], [:a, :b])
        if VERSION >= v"0.7-"
            @test size(A) == (2,)
            @test size(A, 1) == 2
        end
        @test @inferred A[:a] == 1.0
        @test A[:b] == 2.0
        @test length.(Compat.axes(A)) == (2,)
        plus1(x) = x + 1
        B = plus1.(A)
        @test B[:a] == 2.0
        @test B[:b] == 3.0
        @test sprint(show, B) == """
1-dimensional JuMPArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
 2.0
 3.0"""
    end

    @testset "Mixed range/symbol index sets" begin
        A = @inferred JuMPArray([1 2; 3 4], 2:3, [:a, :b])
        if VERSION >= v"0.7-"
            @test size(A) == (2, 2)
            @test size(A, 1) == 2
            @test size(A, 2) == 2
        end
        @test length.(Compat.axes(A)) == (2,2)
        @test @inferred A[2,:a] == 1
        @test A[3,:a] == 3
        @test A[2,:b] == 2
        @test A[3,:b] == 4
        @test A[2,:a,1] == 1
        @test A[2,:a,1,1] == 1
        @test A[3,:a,1,1,1] == 3
        @test @inferred A[:,:a] == JuMPArray([1,3], 2:3)
        @test A[2, :] == JuMPArray([1,2], [:a, :b])
        @test sprint(show, A) == """
2-dimensional JuMPArray{$Int,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Symbol[:a, :b]
And data, a 2×2 Array{$Int,2}:
 1  2
 3  4"""
    end

    @testset "4-dimensional JuMPArray" begin
        if VERSION >= v"0.7-"
            # TODO: This inference tests fails on 0.7. Investigate and fix.
            A = JuMPArray(zeros(2,2,2,2), 2:3, [:a, :b], -1:0, ["a","b"])
        else
            A = @inferred JuMPArray(zeros(2,2,2,2), 2:3, [:a, :b], -1:0,
                                    ["a","b"])
        end
        if VERSION >= v"0.7-"
            @test size(A) == (2, 2, 2, 2)
            @test size(A, 1) == 2
            @test size(A, 2) == 2
            @test size(A, 3) == 2
            @test size(A, 4) == 2
        end
        A[2,:a,-1,"a"] = 1.0
        f = 0.0
        for I in eachindex(A)
            f += A[I]
        end
        @test f == 1.0
        @test isassigned(A, 2, :a, -1, "a")
        @test A[:,:,-1,"a"] == JuMPArray([1.0 0.0; 0.0 0.0], 2:3, [:a,:b])
        @test_throws KeyError A[2,:a,-1,:a]
        if VERSION >= v"0.7-"
            @test sprint(show, A) == """
4-dimensional JuMPArray{Float64,4,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Symbol[:a, :b]
    Dimension 3, -1:0
    Dimension 4, ["a", "b"]
And data, a 2×2×2×2 Array{Float64,4}:
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
        else
            @test sprint(show, A) == """
4-dimensional JuMPArray{Float64,4,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Symbol[:a, :b]
    Dimension 3, -1:0
    Dimension 4, String["a", "b"]
And data, a 2×2×2×2 Array{Float64,4}:
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
        end
    end

    @testset "0-dimensional JuMPArray" begin
        a = Array{Int,0}(undef)
        a[] = 10
        A = JuMPArray(a)
        if VERSION >= v"0.7-"
            @test size(A) == tuple()
        end
        @test A[] == 10
        A[] = 1
        @test sprint(show, A) == """
0-dimensional JuMPArray{$Int,0,...} with index sets:
And data, a 0-dimensional Array{$Int,0}:
1"""
    end

    @testset "JuMPArray keys" begin
        A = JuMPArray([5.0 6.0; 7.0 8.0], 2:3, [:a,:b])
        A_keys = keys(A)
        @test A[A_keys[1,2]] == 6.0
        @test A[A_keys[2,2]] == 8.0
        @test A_keys[1,2][1] == 2
        @test A_keys[1,2][2] == :b
        @test A_keys[2,2][1] == 3
        @test A_keys[2,2][2] == :b
    end
end
