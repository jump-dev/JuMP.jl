@testset "DenseAxisArray" begin
    @testset "undef constructor" begin
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
    end

    @testset "undef constructor (ii)" begin
        A = @inferred DenseAxisArray{String}(undef, 1:2)
        @test !isassigned(A, 1)
        @test !isassigned(A, 2)
        @test !isassigned(A, 3)
        A[1] = "abc"
        @test isassigned(A, 1)
        @test !isassigned(A, 2)
        @test !isassigned(A, 3)
        @test !isassigned(A, 2, 2)
    end

    @testset "Range index set" begin
        A = @inferred DenseAxisArray([1.0,2.0], 2:3)
        @test size(A) == (2,)
        @test size(A, 1) == 2
        @test @inferred A[2] == 1.0
        @test A[3] == 2.0
        @test A[2,1] == 1.0
        @test A[3,1,1,1,1] == 2.0
        @test isassigned(A, 2)
        @test !isassigned(A, 1)
        @test length.(axes(A)) == (2,)
        correct_answer = DenseAxisArray([2.0, 3.0], 2:3)
        @test sprint(show, correct_answer) == """
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, 2:3
And data, a 2-element Array{Float64,1}:
 2.0
 3.0"""
        plus1(x) = x + 1
        @test plus1.(A) == correct_answer
        @test A .+ 1 == correct_answer
        @test 1 .+ A == correct_answer
    end

    @testset "Symbol index set" begin
        A = @inferred DenseAxisArray([1.0,2.0], [:a, :b])
        @test size(A) == (2,)
        @test size(A, 1) == 2
        @test @inferred A[:a] == 1.0
        @test A[:b] == 2.0
        @test length.(axes(A)) == (2,)
        correct_answer = DenseAxisArray([2.0, 3.0], [:a, :b])
        @test sprint(show, correct_answer) == """
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
 2.0
 3.0"""
        plus1(x) = x + 1
        @test plus1.(A) == correct_answer
        @test A .+ 1 == correct_answer
        @test 1 .+ A == correct_answer
    end

    @testset "Mixed range/symbol index sets" begin
        A = @inferred DenseAxisArray([1 2; 3 4], 2:3, [:a, :b])
        @test size(A) == (2, 2)
        @test size(A, 1) == 2
        @test size(A, 2) == 2
        @test length.(axes(A)) == (2,2)
        @test @inferred A[2,:a] == 1
        @test A[3,:a] == 3
        @test A[2,:b] == 2
        @test A[3,:b] == 4
        @test A[2,:a,1] == 1
        @test A[2,:a,1,1] == 1
        @test A[3,:a,1,1,1] == 3
        @test @inferred A[:,:a] == DenseAxisArray([1,3], 2:3)
        @test A[2, :] == DenseAxisArray([1,2], [:a, :b])
        @test sprint(show, A) == """
2-dimensional DenseAxisArray{$Int,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Symbol[:a, :b]
And data, a 2×2 Array{$Int,2}:
 1  2
 3  4"""
    end

    @testset "4-dimensional DenseAxisArray" begin
        # TODO: This inference tests fails on 0.7. Investigate and fix.
        A = DenseAxisArray(zeros(2,2,2,2), 2:3, [:a, :b], -1:0, ["a","b"])
        @test size(A) == (2, 2, 2, 2)
        @test size(A, 1) == 2
        @test size(A, 2) == 2
        @test size(A, 3) == 2
        @test size(A, 4) == 2
        A[2,:a,-1,"a"] = 1.0
        f = 0.0
        for I in eachindex(A)
            f += A[I]
        end
        @test f == 1.0
        @test isassigned(A, 2, :a, -1, "a")
        @test A[:,:,-1,"a"] == DenseAxisArray([1.0 0.0; 0.0 0.0], 2:3, [:a,:b])
        @test_throws KeyError A[2,:a,-1,:a]
        @test sprint(show, A) == """
4-dimensional DenseAxisArray{Float64,4,...} with index sets:
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
    end

    @testset "0-dimensional DenseAxisArray" begin
        a = Array{Int,0}(undef)
        a[] = 10
        A = DenseAxisArray(a)
        @test size(A) == tuple()
        @test A[] == 10
        A[] = 1
        @test sprint(show, A) == """
0-dimensional DenseAxisArray{$Int,0,...} with index sets:
And data, a 0-dimensional Array{$Int,0}:
1"""
    end

    @testset "DenseAxisArray keys" begin
        A = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, [:a,:b])
        A_keys = collect(keys(A))
        @test A[A_keys[3]] == 6.0
        @test A[A_keys[4]] == 8.0
        @test A_keys[3][1] == 2
        @test A_keys[3][2] == :b
        @test A_keys[4][1] == 3
        @test A_keys[4][2] == :b

        B = DenseAxisArray([5.0 6.0; 7.0 8.0], 2:3, Set([:a,:b]))
        B_keys = keys(B)
        @test Containers.DenseAxisArrayKey((2, :a)) in B_keys
        @test Containers.DenseAxisArrayKey((2, :b)) in B_keys
        @test Containers.DenseAxisArrayKey((3, :a)) in B_keys
        @test Containers.DenseAxisArrayKey((3, :b)) in B_keys
    end
end
