if VERSION < v"1.0.1"
    # Workaround for JuliaLang/julia/pull/28625
    if Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end
end

using Documenter, JuMP

makedocs(
    sitename = "JuMP",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    authors = "Miles Lubin, Iain Dunning, and Joey Huchette",
    analytics = "UA-44252521-1",
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
        "Style Guide" => "style.md",
        "Extensions" => "extensions.md",
        "Development Roadmap" => "roadmap.md"
    ],
    assets = [
        "assets/jump-logo-with-text.svg",
        "assets/numfocus-logo.png"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/JuMP.jl.git",
)
