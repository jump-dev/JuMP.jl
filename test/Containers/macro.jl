using Test
using JuMP
using JuMP.Containers

@testset "Macro" begin
    @testset "Array" begin
        Containers.@container(x[i = 1:3], i^2)
        @test x isa Vector{Int}
        x = Containers.@container([i = 1:3, j = 1:3], i^2)
        @test x isa Matrix{Int}
    end
    @testset "DenseAxisArray" begin
        Containers.@container(x[i = 2:3], i^2)
        @test x isa Containers.DenseAxisArray{Int, 1}
        Containers.@container(x[i = 2:3, j = 1:2], i + j)
        @test x isa Containers.DenseAxisArray{Int, 2}
        Containers.@container(x[4], 0.0)
        @test x isa Containers.DenseAxisArray{Float64, 1}
        Containers.@container(x[4, 5], 0)
        @test x isa Containers.DenseAxisArray{Int, 2}
        Containers.@container(x[4, 1:3, 5], 0)
        @test x isa Containers.DenseAxisArray{Int, 3}
    end
    @testset "SparseAxisArray" begin
        Containers.@container(x[i = 1:3, j = 1:i], i + j)
        @test x isa Containers.SparseAxisArray{Int, 2}
        Containers.@container(x[i=1:10; iseven(i)], i)
        @test x isa Containers.SparseAxisArray{Int, 1}
        Containers.@container(x[i=1:0, j=i:0], i)
        @test x isa SparseAxisArray{Any,2,Tuple{Any,Any}}
        Containers.@container(x[i=1:2, j=1:2; false], i)
        @test x isa SparseAxisArray{Any,2,Tuple{Any,Any}}
        Containers.@container(x[i=1:0, j=2:1], i, container = SparseAxisArray)
        @test x isa SparseAxisArray{Any,2,Tuple{Int,Int}}
        Containers.@container(x[i=1:0, j=1:0], i, container = SparseAxisArray)
        @test x isa SparseAxisArray{Any,2,Tuple{Int,Int}}
    end
end
