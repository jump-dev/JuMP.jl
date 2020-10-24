using Test

@testset "Containers" begin
    @testset "$(file)" for file in filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
        if file in [
            "Containers.jl",
        ]
            continue
        end
        include(joinpath(@__DIR__, file))
    end
end
