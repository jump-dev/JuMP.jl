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
    DOCTYPE_STRUCT,
)

struct _Config
    current_module::Module
    root::String
    subdirectory::String
    extras::Vector{Pair{String,DocType}}
end

const CONFIG = _Config[]

abstract type APIBuilder <: Documenter.Builder.DocumentPipeline end

Documenter.Selectors.order(::Type{APIBuilder}) = 0.0

"""
    automatic_reference_documentation(
        current_module::Module;
        root::String,
        subdirectory::String,
        extras::Vector{Pair{String,DocType}} = Pair{String,DocType}[],
    )

Automatically creates the API reference documentation for `current_module` and
returns a `Vector` which can be used in the `pages` argument of
`Documenter.makedocs`.

## Arguments

 * `current_module`: the module from which to create an API reference.
 * `subdirectory`: the directory relative to the documentation root in which to
   write the API files.
 * `extras`: a vector of non-exported docstrings to include in the API
   reference. Each element is a pair which maps the docstring signature to a
   [`DocumenterReference.DocType`](@ref) enum.

## Multiple instances

Each time you call this function, a new object is added to the global variable
`DocumenterReference.CONFIG`.
"""
function automatic_reference_documentation(
    current_module::Module;
    root::String,
    subdirectory::String,
    extras::Vector{Pair{String,DocType}} = Pair{String,DocType}[],
)
    config = _Config(current_module, root, subdirectory, extras)
    pages = Any["Overview"=>"$subdirectory/$current_module.md"]
    _iterate_over_symbols(config) do key, type
        push!(pages, Documenter.hide("`$key`" => "$subdirectory/$key.md"))
        return
    end
    push!(CONFIG, config)
    return pages
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
            # Skip
        else
            push!(contents, n => DOCTYPE_CONSTANT)
        end
    end
    order = Dict(
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
    for (key, type) in vcat(_exported_symbols(current_module), config.extras)
        if key isa Symbol
            doc = Base.Docs.doc(Base.Docs.Binding(current_module, key))
            if occursin("No documentation found.", string(doc))
                if type == DOCTYPE_CONSTANT
                    continue  # It's okay not to document every constant.
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
    elseif x == DOCTYPE_STRUCT
        return "struct"
    end
end

function _add_page(document::Documenter.Document, filename, contents)
    mdpage = Markdown.parse(contents)
    document.blueprint.pages[filename] = Documenter.Page(
        filename, # source, gets ignored because of `EditURL = nothing`.
        "$(document.user.build)/$filename", # build,
        "$(document.user.build)/", # workdir,
        mdpage.content,
        Documenter.Globals(),
        convert(MarkdownAST.Node, mdpage),
    )
    return
end

function _build_api_page(document::Documenter.Document, config::_Config)
    subdir = config.subdirectory
    overview_md = """
    ```@meta
    EditURL = nothing
    ```

    # [API](@id DocumenterReference_$(config.current_module))

    | NAME | KIND |
    | :--- | :--- |
    """
    _iterate_over_symbols(config) do key, type
        _add_page(
            document,
            "$subdir/$key.md",
            """
            ```@meta
            EditURL = nothing
            ```

            ```@docs
            $key
            ```
            """,
        )
        overview_md *= "| [`$key`](@ref) | `$(_to_string(type))` |\n"
        return
    end
    _add_page(document, "$subdir/$(config.current_module).md", overview_md)
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
