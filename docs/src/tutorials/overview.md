# Tutorial overview

To make it easier to find relevant tutorials, this page lists all the tutorials
by the packages they import.

For example, to find tutorials that solve nonlinear programs, look at the
[Ipopt.jl](@ref tutorial_Ipopt) section.

```@eval
package_to_headers = Dict{String,Vector{String}}()
for (root, dirs, files) in walkdir(@__DIR__)
    for file in filter(f -> endswith(f, ".jl") | endswith(f, ".md"), files)
        contents = read(joinpath(root, file), String)
        if startswith(contents, "```@meta")
            continue  # Already literated
        end
        m = match(r"\n# # (.+)\n", contents)
        if m === nothing
            m = match(r"^# (.+)\n", contents)
        end
        title = m[1]
        anchor = "[$title](@ref)"
        if occursin("@id", title)
            m = match(r"\[(.+)\]\(@id (.+)\)", title)
            anchor = "[$(m[1])](@ref $(m[2]))"
        end
        for m in eachmatch(r"\nimport ([a-z][A-Z]+)\n"i, contents; overlap = true)
            push!(get!(package_to_headers, m[1], String[]), anchor)
        end
    end
end
io = IOBuffer()
for (pkg, tutorials) in sort!(collect(package_to_headers))
    println(io, "## [$(pkg).jl](@id tutorial_$pkg)\n", join(sort!(tutorials), ", "))
end
seekstart(io)
import Markdown
Markdown.parse(read(io, String))
```
