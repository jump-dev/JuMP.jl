using Documenter
using Literate
using JuMP
using Test

const _EXAMPLE_DIR = joinpath(@__DIR__, "src", "examples")

for file in readdir(_EXAMPLE_DIR)
    if !endswith(file, ".jl")
        continue
    end
    filename = joinpath(_EXAMPLE_DIR, file)
    # `include` the file to test it before `#src` lines are removed. It is in a
    # testset to isolate local variables between files.
    @testset "$(file)" begin
        include(filename)
    end
    Literate.markdown(filename, _EXAMPLE_DIR; documenter = true)
end

makedocs(
    sitename = "JuMP",
    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "UA-44252521-1",
    ),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    authors = "Miles Lubin, Iain Dunning, and Joey Huchette",
    pages = [
        "Introduction" => "index.md",
        "Installation Guide" => "installation.md",
        "Quick Start Guide" => "quickstart.md",
        "Variables" => "variables.md",
        "Expressions" => "expressions.md",
        "Objective" => "objective.md",
        "Constraints" => "constraints.md",
        "Containers" => "containers.md",
        "Solvers" => "solvers.md",
        "Query Solutions" => "solutions.md",
        "Nonlinear Modeling" => "nlp.md",
        "Callbacks" => "callbacks.md",
        "Style Guide" => "style.md",
        "Extensions" => "extensions.md",
        "Development Roadmap" => "roadmap.md",
        "Examples" => map(
            file -> joinpath("examples", file),
            filter(
                file -> endswith(file, ".md"),
                sort(readdir(_EXAMPLE_DIR)),
            )
        ),
    ],
)

deploydocs(
    repo   = "github.com/jump-dev/JuMP.jl.git",
)
