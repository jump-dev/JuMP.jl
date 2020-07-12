using Test

@testset "Containers" begin
    @testset "$(file)" for file in filter(f -> endswith(f, ".jl"), readdir(@__DIR__))
        if file in [
            "Containers.jl",
        ]
            continue
        end
        filename = joinpath(@__DIR__, file)
        t = time()
        include(filename)
        println("$(filename) took $(round(time() - t; digits = 1)) seconds.")
    end
end
