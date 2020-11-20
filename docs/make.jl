using Documenter
using Literate
using JuMP

const _EXAMPLE_INPUT_DIR = joinpath(dirname(@__DIR__), "examples")
const _EXAMPLE_OUTPUT_DIR = joinpath(@__DIR__, "src", "examples")

# Delete all files in the output directory to make way for fresh versions.
for file in readdir(_EXAMPLE_OUTPUT_DIR)
    if file == ".gitkeep"
        continue
    end
    rm(joinpath(_EXAMPLE_OUTPUT_DIR, file))
end

# Replace with this for-loop once all examples are migrated.
# for file in sort(readdir(_EXAMPLE_INPUT_DIR))
for file in [
    "basic.jl",
    "callbacks.jl",
    "cannery.jl",
    "clnlbeam.jl",
    "mle.jl"
]
    if !endswith(file, ".jl")
        continue
    end
    Literate.markdown(
        joinpath(_EXAMPLE_INPUT_DIR, file),
        _EXAMPLE_OUTPUT_DIR;
        documenter = true,
        execute = true,
    )
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
           sort(readdir(_EXAMPLE_OUTPUT_DIR)),
        ),
    ],
)

deploydocs(
    repo   = "github.com/jump-dev/JuMP.jl.git",
)
