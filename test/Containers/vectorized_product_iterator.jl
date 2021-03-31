using JuMP.Containers
using Test

@testset "Vectorized Product Iterator" begin
    I = [
        1 2
        3 4
    ]
    @test isempty(Containers.vectorized_product(2, I, 1:0))
    @test isempty(collect(Containers.vectorized_product(2, I, 1:0)))
    @test collect(Containers.vectorized_product(2, I)) ==
          [(2, 1) (2, 3) (2, 2) (2, 4)]
end

@testset "Unknown size" begin
    f = Iterators.filter(k -> isodd(k), 1:10)
    v = Containers.vectorized_product(f)
    @test axes(v) == (Base.OneTo(5),)
end
