using Compat
using Compat.Test
using JuMP

@testset "Containers" begin
    @testset "SparseAxisArray" begin
        SA = JuMP.Containers.SparseAxisArray
        d = @inferred SA(Dict((:a,) => 1, (:b,) => 2))
        @test d isa SA{Int, 1, Tuple{Symbol}}
        sqr(x) = x^2
        @testset "Colon indexing" begin
            if VERSION < v"0.7-"
                @test_throws ArgumentError d[:, 1]
                @test_throws ArgumentError d[:a, :]
            else
                err = ArgumentError("Indexing with `:` is not supported by" *
                                    " Containers.SparseAxisArray")
                @test_throws err d[:, 1]
                @test_throws err d[:a, :]
            end
        end
        @testset "Map" begin
            @test (@inferred map(x -> x * 3, d)) == SA(Dict((:a,) => 3, (:b,) => 6))
            @test (@inferred map(x -> 3 * x, d)) == SA(Dict((:a,) => 3, (:b,) => 6))
            @test (@inferred map(identity, d)) == d
            @test (@inferred map(sqr, d)) == SA(Dict((:a,) => 1, (:b,) => 4))
        end
        @testset "Reduce" begin
            @test 3 == @inferred sum(d)
        end
        @testset "Broadcasting" begin
            @test (@inferred d .* d) == SA(Dict((:a,) => 1, (:b,) => 4))
            @test (@inferred d .+ d) == SA(Dict((:a,) => 2, (:b,) => 4))
            @test (@inferred d .* 3) == SA(Dict((:a,) => 3, (:b,) => 6))
            @test (@inferred 3 .* d) == SA(Dict((:a,) => 3, (:b,) => 6))
            @test identity.(d) == d
            @test sqr.(d) == SA(Dict((:a,) => 1, (:b,) => 4))
            @testset "Different array" begin
                if VERSION < v"0.7-"
                    @test_throws ArgumentError [1, 2] .+ d
                    @test_throws ArgumentError d .* [1, 2]
                else
                    err = ArgumentError("Cannot broadcast" *
                                        " Containers.SparseAxisArray with" *
                                        " another array of different type")
                    @test_throws err [1, 2] .+ d
                    @test_throws err d .* [1, 2]
                end
            end
            @testset "Different indices" begin
                dc = @inferred SA(Dict((:a,) => 1, (:b,) => 2, (:c,) => 3))
                da = @inferred SA(Dict((:b,) => 2))
                if VERSION < v"0.7-"
                    @test_throws ArgumentError dc .+ d
                    @test_throws ArgumentError d .+ dc
                    @test_throws ArgumentError da .+ d
                    @test_throws ArgumentError d .+ da
                else
                    err = ArgumentError("Cannot broadcast" *
                        " Containers.SparseAxisArray with different indices")
                    @test_throws err dc .+ d
                    @test_throws err d .+ dc
                    @test_throws err da .+ d
                    @test_throws err d .+ da
                end
            end
        end
    end
end
