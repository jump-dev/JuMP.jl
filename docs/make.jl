using Documenter, JuMP

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
        "Development Roadmap" => "roadmap.md"
    ],
)

deploydocs(
    repo   = "github.com/jump-dev/JuMP.jl.git",
)
