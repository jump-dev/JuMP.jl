@testset "Nested Iterator" begin
    iterators = (() -> 1:3, i -> 1:i)
    condition(i, j) = j > i
    @test isempty(Containers.nested(iterators..., condition=condition))
    @test isempty(collect(Containers.nested(iterators..., condition=condition)))
    condition(i, j) = isodd(i) || isodd(j)
    @test collect(Containers.nested(iterators..., condition=condition)) == [
        (1, 1)
        (2, 1)
        (3, 1)
        (3, 2)
        (3, 3)
    ]
end
