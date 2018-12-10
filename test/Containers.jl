using Compat
using Compat.Test
using JuMP
using JuMP.Containers

@testset "Containers" begin
    include("DenseAxisArray.jl")
    include("SparseAxisArray.jl")
end
