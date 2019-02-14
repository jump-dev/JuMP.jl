@testset "SparseAxisArray" begin
    function sparse_test(d, sum_d, d2, d3, dsqr, d_bads)
        sqr(x) = x^2
        @testset "Colon indexing" begin
            err = ArgumentError("Indexing with `:` is not supported by" *
                                " Containers.SparseAxisArray")
            @test_throws err d[:, 1]
            @test_throws err d[:a, :]
        end
        @testset "Map" begin
            @test d == @inferred map(identity, d)
            @test dsqr == @inferred map(sqr, d)
            @test d3 == @inferred map(x -> x * 3, d)
            @test d3 == @inferred map(x -> 3 * x, d)
        end
        @testset "Reduce" begin
            @test sum_d == @inferred sum(d)
        end
        @testset "Broadcasting" begin
            @test dsqr == @inferred d .* d
            @test d2 == @inferred d .+ d
            @test d3 == @inferred d .* 3
            @test d3 == @inferred 3 .* d
            @test d == identity.(d)
            @test dsqr == sqr.(d)
            @testset "Different array" begin
                err = ArgumentError("Cannot broadcast" *
                                    " Containers.SparseAxisArray with" *
                                    " another array of different type")
                @test_throws err [1, 2] .+ d
                @test_throws err d .* [1, 2]
            end
            @testset "Different indices" begin
                err = ArgumentError("Cannot broadcast" *
                                    " Containers.SparseAxisArray with " *
                                    "different indices")
                for d_bad in d_bads
                    @test_throws err d_bad .+ d
                    @test_throws err d .+ d_bad
                end
            end
        end
    end
    @testset "1-dimensional" begin
        SA = SparseAxisArray
        d = @inferred SA(Dict((:a,) => 1, (:b,) => 2))
        @testset "Printing" begin
            @test sprint(show, "text/plain", d) == """
SparseAxisArray{$Int,1,Tuple{Symbol}} with 2 entries:
  [a]  =  1
  [b]  =  2"""
        end
        @test d isa SA{Int, 1, Tuple{Symbol}}
        d2 = @inferred SA(Dict((:a,) => 2, (:b,) => 4))
        d3 = @inferred SA(Dict((:a,) => 3, (:b,) => 6))
        dsqr = @inferred SA(Dict((:a,) => 1, (:b,) => 4))
        da = @inferred SA(Dict((:b,) => 2))
        dc = @inferred SA(Dict((:a,) => 1, (:b,) => 2, (:c,) => 3))
        sparse_test(d, 3, d2, d3, dsqr, [da, dc])
        @testset "Broadcasting with type unstability" begin
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
        end
    end
    @testset "2-dimensional" begin
        SA = SparseAxisArray
        d = @inferred SA(Dict((:a, 'u') => 2.0, (:b, 'v') => 0.5))
        @test d isa SA{Float64, 2, Tuple{Symbol, Char}}
        @testset "Printing" begin
            @test sprint(show, "text/plain", d) == """
SparseAxisArray{Float64,2,Tuple{Symbol,Char}} with 2 entries:
  [b, v]  =  0.5
  [a, u]  =  2.0"""
        end
        d2 = @inferred SA(Dict((:b, 'v') => 1.0, (:a, 'u') => 4.0))
        d3 = @inferred SA(Dict((:a, 'u') => 6.0, (:b, 'v') => 1.5))
        dsqr = @inferred SA(Dict((:a, 'u') => 4.0, (:b, 'v') => 0.25))
        da = @inferred SA(Dict((:b, 'v') => 2.0))
        db = @inferred SA(Dict((:a, 'u') => 1.0, (:b, 'u') => 2.0))
        dc = @inferred SA(Dict((:a, 'u') => 1.0, (:b, 'v') => 2.0,
                               (:c, 'w') => 3.0))
        sparse_test(d, 2.5, d2, d3, dsqr, [da, db, dc])
    end
end
