if VERSION < v"1.0.1"
    # Workaround for JuliaLang/julia/pull/28625
    if Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end
end

using Documenter, JuMP

makedocs(
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
        "Containers" => "containers.md",
        "Names" => "names.md",
        "Solvers" => "solvers.md",
        "Nonlinear Modeling" => "nlp.md",
        "Style Guide" => "style.md",
        "Extensions" => "extensions.md",
        "Updating Guide" => "updating.md",
        "How do I ...? (FAQ)" => "howdoi.md"
    ],
    assets = [
        "assets/jump-logo-with-text.svg",
        "assets/numfocus-logo.png"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/JuMP.jl.git",
)
