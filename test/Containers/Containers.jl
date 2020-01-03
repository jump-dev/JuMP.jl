using Test
using JuMP
using JuMP.Containers

@testset "Containers" begin
    include("DenseAxisArray.jl")
    include("SparseAxisArray.jl")
    include("generate_container.jl")
    include("vectorized_product_iterator.jl")
    include("nested_iterator.jl")
    include("no_duplicate_dict.jl")
    include("macro.jl")
end
