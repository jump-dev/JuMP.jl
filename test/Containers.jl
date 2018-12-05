using Compat
using Compat.Test
using JuMP

@testset "Containers" begin
    @testset "SparseArray" begin
        SAInt1 = JuMP.Containers.SparseArray{Int, 1}
        d = SAInt1(Dict((:a,) => 1, (:b,) => 2))
        @testset "Reduce" begin
            @test sum(d) == 3
        end
        @testset "Broadcasting" begin
            @test d .* d == SAInt1(Dict((:a,) => 1, (:b,) => 4))
            @test d .+ d == SAInt1(Dict((:a,) => 2, (:b,) => 4))
            @test d .* 3 == SAInt1(Dict((:a,) => 3, (:b,) => 6))
            @test 3 .* d == SAInt1(Dict((:a,) => 3, (:b,) => 6))
            @test identity.(d) == d
            sqr(x) = x^2
            @test sqr.(d) == SAInt1(Dict((:a,) => 1, (:b,) => 4))
            @testset "Different indices" begin
                dc = SAInt1(Dict((:a,) => 1, (:b,) => 2, (:c,) => 3))
                da = SAInt1(Dict((:b,) => 2))
                if VERSION < v"0.7-"
                    @test_throws ArgumentError dc .+ d
                    @test_throws ArgumentError d .+ dc
                    @test_throws ArgumentError da .+ d
                    @test_throws ArgumentError d .+ da
                else
                    err = ArgumentError("Cannot broadcast" *
                        " Containers.SparseArray with different indices")
                    @test_throws err dc .+ d
                    @test_throws err d .+ dc
                    @test_throws err da .+ d
                    @test_throws err d .+ da
                end
            end
        end
    end
end
