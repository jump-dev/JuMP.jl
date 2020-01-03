@testset "Vectorized Product Iterator" begin
    I = [1 2
         3 4]
    @test isempty(Containers.vectorized_product(2, I, 1:0))
    @test isempty(collect(Containers.vectorized_product(2, I, 1:0)))
    @test collect(Containers.vectorized_product(2, I)) == [
        (2, 1) (2, 3) (2, 2) (2, 4)
    ]
end
