#  Copyright 2023, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module DocumenterReference

import Documenter
import Markdown
import MarkdownAST

@enum(
    DocType,
    DOCTYPE_ABSTRACT_TYPE,
    DOCTYPE_CONSTANT,
    DOCTYPE_FUNCTION,
    DOCTYPE_MACRO,
    DOCTYPE_MODULE,
    DOCTYPE_STRUCT,
)

struct _Config
    current_module::Module
    root::String
    subdirectory::String
    modules::Dict{Module,<:Vector}
    sort_by::Function
end

const CONFIG = _Config[]

abstract type APIBuilder <: Documenter.Builder.DocumentPipeline end

Documenter.Selectors.order(::Type{APIBuilder}) = 0.0

"""
    automatic_reference_documentation(;
        root::String,
        subdirectory::String,
        modules::Dict{Module,Vector{Pair{String,DocType}}},
        sort_by::Function = identity,
    )

Automatically creates the API reference documentation for `current_module` and
returns a `Vector` which can be used in the `pages` argument of
`Documenter.makedocs`.

## Arguments

 * `current_module`: the module from which to create an API reference.
 * `subdirectory`: the directory relative to the documentation root in which to
   write the API files.
 * `modules`: a dictionary mapping modules to a vector of non-exported
   docstrings to include in the API reference. Each element is a pair which maps
   the docstring signature to a [`DocumenterReference.DocType`](@ref) enum.

## Multiple instances

Each time you call this function, a new object is added to the global variable
`DocumenterReference.CONFIG`.
"""
function automatic_reference_documentation(;
    root::String,
    subdirectory::String,
    modules::Vector,
    sort_by::Function = identity,
)
    _to_extras(m::Module) = m => Any[]
    _to_extras(m::Pair) = m
    _modules = Dict(_to_extras(m) for m in modules)
    list_of_pages = Any[]
    for m in modules
        current_module = first(_to_extras(m))
        pages = _automatic_reference_documentation(
            current_module;
            root,
            subdirectory,
            modules = _modules,
            sort_by,
        )
        push!(list_of_pages, "$current_module" => pages)
    end
    return "API Reference" => list_of_pages
end

function _automatic_reference_documentation(
    current_module::Module;
    root::String,
    subdirectory::String,
    modules::Dict{Module,<:Vector},
    sort_by::Function,
)
    push!(CONFIG, _Config(current_module, root, subdirectory, modules, sort_by))
    return "$subdirectory/$current_module.md"
end

function _exported_symbols(mod)
    contents = Pair{Symbol,DocType}[]
    for n in names(mod)
        f = getfield(mod, n)
        f_str = string(f)
        if startswith(f_str, "@")
            push!(contents, n => DOCTYPE_MACRO)
        elseif startswith(f_str, "Abstract")
            push!(contents, n => DOCTYPE_ABSTRACT_TYPE)
        elseif f isa Type
            push!(contents, n => DOCTYPE_STRUCT)
        elseif f isa Function
            if islowercase(f_str[1])
                push!(contents, n => DOCTYPE_FUNCTION)
            else
                push!(contents, n => DOCTYPE_STRUCT)
            end
        elseif f isa Module
            push!(contents, n => DOCTYPE_MODULE)
        else
            push!(contents, n => DOCTYPE_CONSTANT)
        end
    end
    order = Dict(
        DOCTYPE_MODULE => 0,
        DOCTYPE_MACRO => 1,
        DOCTYPE_FUNCTION => 2,
        DOCTYPE_ABSTRACT_TYPE => 3,
        DOCTYPE_STRUCT => 4,
        DOCTYPE_CONSTANT => 5,
    )
    return sort(contents; by = x -> (order[x[2]], "$(x[1])"))
end

function _iterate_over_symbols(f, config)
    current_module = config.current_module
    modules = get(config.modules, config.current_module, Any[])
    key_types = vcat(_exported_symbols(current_module), modules)
    for (key, type) in sort!(key_types; by = config.sort_by)
        if key isa Symbol
            doc = Base.Docs.doc(Base.Docs.Binding(current_module, key))
            if occursin("No documentation found.", string(doc))
                if type == DOCTYPE_MODULE
                    mod = getfield(current_module, key)
                    if mod == current_module || !haskey(config.modules, mod)
                        continue
                    end
                else
                    error("Documentation missing for $key")
                end
            end
        end
        f(key, type)
    end
    return
end

function _to_string(x::DocType)
    if x == DOCTYPE_ABSTRACT_TYPE
        return "abstract type"
    elseif x == DOCTYPE_CONSTANT
        return "constant"
    elseif x == DOCTYPE_FUNCTION
        return "function"
    elseif x == DOCTYPE_MACRO
        return "macro"
    elseif x == DOCTYPE_MODULE
        return "module"
    elseif x == DOCTYPE_STRUCT
        return "struct"
    end
end

function _build_api_page(document::Documenter.Document, config::_Config)
    subdir = config.subdirectory
    overview_md = """
    ```@meta
    EditURL = nothing
    ```

    # [$(config.current_module)](@id DocumenterReference_$(config.current_module))

    This page lists the public API of `$(config.current_module)`.

    !!! info
        This page is an unstructured list of the $(config.current_module) API. For a
        more structured overview, read the Manual or Tutorial parts of this
        documentation.

    Load all of the public the API into the current scope with:
    ```julia
    using $(config.current_module)
    ```
    Alternatively, load only the module with:
    ```julia
    import $(config.current_module)
    ```
    and then prefix all calls with `$(config.current_module).` to create
    `$(config.current_module).<NAME>`.
    """
    list_of_docstrings = String[]
    _iterate_over_symbols(config) do key, type
        if type == DOCTYPE_MODULE
            return
        end
        push!(
            list_of_docstrings,
            "## `$key`\n\n```@docs\n$(config.current_module).$key\n```\n\n",
        )
        return
    end
    md_page = Markdown.parse(overview_md * join(list_of_docstrings, "\n"))
    filename = "$subdir/$(config.current_module).md"
    document.blueprint.pages[filename] = Documenter.Page(
        joinpath(document.user.source, filename),
        joinpath(document.user.build, filename),
        document.user.build,
        md_page.content,
        Documenter.Globals(),
        convert(MarkdownAST.Node, md_page),
    )
    return
end

function Documenter.Selectors.runner(
    ::Type{APIBuilder},
    document::Documenter.Document,
)
    @info "APIBuilder: creating API reference"
    for config in CONFIG
        _build_api_page(document, config)
    end
    return
end

end  # module
