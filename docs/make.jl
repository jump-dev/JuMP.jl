using Documenter
using Literate
using JuMP
using Test

const _EXAMPLE_DIR = joinpath(@__DIR__, "src", "examples")

function link_example(content)
    edit_url = match(r"EditURL = \"(.+?)\"", content)[1]
    footer = match(r"^(---\n\n\*This page was generated using)"m, content)[1]
    content = replace(
        content, footer => "[View this file on Github]($(edit_url)).\n\n" * footer
    )
    return content
end

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
    Literate.markdown(
        filename, _EXAMPLE_DIR; documenter = true, postprocess = link_example
    )
end

makedocs(
    sitename = "JuMP",
    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "UA-44252521-1",
        collapselevel = 1,
    ),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    authors = "Miles Lubin, Iain Dunning, and Joey Huchette",
    pages = [
        "Introduction" => "index.md",
        "Manual" => [
            "installation.md",
            "quickstart.md",
            "variables.md",
            "expressions.md",
            "objective.md",
            "constraints.md",
            "containers.md",
            "solvers.md",
            "solutions.md",
            "nlp.md",
            "callbacks.md",
            "Extensions" => "extensions.md",
        ],
        "API Reference" => map(
            file -> joinpath("reference", file),
            sort(readdir(joinpath(@__DIR__, "src", "reference"))),
        ),
        "Examples" => map(
            file -> joinpath("examples", file),
            filter(
                file -> endswith(file, ".md"),
                sort(readdir(_EXAMPLE_DIR)),
            )
        ),
        "Style Guide" => "style.md",
        "Development Roadmap" => "roadmap.md",
    ],
)

deploydocs(
    repo   = "github.com/jump-dev/JuMP.jl.git",
)
