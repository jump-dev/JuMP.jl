# TODO: Copy over tests that are still relevant from old/macros.jl.

@testset "Macros" begin
    @testset "Nested tuple destructuring" begin
        m = Model()
        d = Dict((1,2) => 3)
        ex = @expression(m, sum(i+j+k for ((i,j),k) in d))
        @test ex == 6
    end
end
