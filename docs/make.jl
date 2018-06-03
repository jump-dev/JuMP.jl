using Documenter, JuMP

makedocs(
    format = :html,
    sitename = "JuMP",
    authors = "Miles Lubin, Iain Dunning, and Joey Huchette",
    analytics = "UA-44252521-1",
    pages = [
        "Introduction" => "index.md",
        "Installation Guide" => "installation.md",
        "Quick Start Guide" => "quickstart.md",
        "Concepts" => "concepts.md",
        "Variables" => "variables.md",
        "Expressions" => "expressions.md",
        "Constraints" => "constraints.md",
        "Solvers" => "solvers.md",
        "Style Guide" => "style.md",
        "Extensions" => "extensions.md",
        "Nonlinear Modeling" => "nlp.md",
        "Updating Guide" => "updating.md",
        "How do I ...? (FAQ)" => "howdoi.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/JuMP.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
