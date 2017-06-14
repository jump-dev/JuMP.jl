using Documenter, JuMP

makedocs(
    format = :html,
    sitename = "JuMP",
    pages = [
        "Constraints" => "constraints.md"
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
