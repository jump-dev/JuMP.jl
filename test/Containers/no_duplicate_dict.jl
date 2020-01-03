@testset "Iterator with constant eltype" begin
    f(ij) = ij => sum(ij)
    g = Base.Generator(f, Iterators.product(1:2, 1:2))
    dict = Containers.NoDuplicateDict(g)
    @test length(dict) == 4
    for i in 1:2
        for j in 1:2
            @test dict[(i, j)] == i + j
        end
    end
    err = ErrorException("Repeated index (1, 2). Index sets must have unique elements.")
    g = Base.Generator(f, [(1, 1), (1, 2), (1, 2)])
    @test_throws err Containers.NoDuplicateDict(g)
    g = Base.Generator(f, [(1, 2), (1, 2)])
    @test_throws err Containers.NoDuplicateDict(g)
end
@testset "Iterator with varying eltype" begin
    f(ij) = ij => ==(ij...) ? nothing : 1.0
    g = Base.Generator(f, Iterators.product(1:2, 1:2))
    dict = Containers.NoDuplicateDict(g)
    @test length(dict) == 4
    for i in 1:2
        for j in 1:2
            if i == j
                @test dict[(i, j)] === nothing
            else
                @test dict[(i, j)] == 1.0
            end
        end
    end
    g = Base.Generator(f, [(1, 1), (1, 2), (1, 1)])
    err = ErrorException("Repeated index (1, 1). Index sets must have unique elements.")
    @test_throws err Containers.NoDuplicateDict(g)
    g = Base.Generator(f, [(1, 1), (1, 1)])
    @test_throws err Containers.NoDuplicateDict(g)
end
