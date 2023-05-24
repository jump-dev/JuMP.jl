module TestContainersDimensionalData

if isdefined(Base,:get_extension)
    using JuMP, DimensionalData
    @testset "DimensionalData.jl integration" begin 
        model = JuMP.Model();
        @test @variable(model, x[i=2:4, j=["a", "b"]], container = DimArray) isa DimArray
        @test_throws ArgumentError @variable(model, z[i=1:3, j=i:2], container = DimArray)
        @test_throws ArgumentError @variable(model, w[1:3, j=1:2], container = DimArray)
    end   
end

end