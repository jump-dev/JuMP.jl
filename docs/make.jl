#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

import Documenter
import DocumenterCitations
import Downloads
import Literate
import MathOptInterface
import Pkg
import Test
import TOML

using JuMP
using JuMP.Containers

# Pass --fast as an argument to skip rebuilding the examples and running
# doctests. Only use this argument to rapidly test small changes to the
# Markdown. _Never_ set it in production.
const _FAST = findfirst(isequal("--fast"), ARGS) !== nothing

# Pass --fix to run with `doctests=:fix`.
const _FIX = findfirst(isequal("--fix"), ARGS) !== nothing

# A flag to check if we are running in a GitHub action.
const _IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Pass --pdf to build the PDF. On GitHub actions, we always build the PDF.
const _PDF = findfirst(isequal("--pdf"), ARGS) !== nothing || _IS_GITHUB_ACTIONS

# ==============================================================================
#  Run literate.jl
# ==============================================================================

function _link_example(content)
    edit_url = match(r"EditURL = \"(.+?)\"", content)[1]
    if !_IS_GITHUB_ACTIONS
        # The link won't work locally. So hard-code in a URL.
        edit_url = replace(
            edit_url,
            "<unknown>" => "https://github.com/jump-dev/JuMP.jl/tree/master",
        )
    end
    title_line = findfirst(r"\n# .+?\n", content)
    line = content[title_line]
    new_title = string(
        line,
        "\n",
        "_This tutorial was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl)._\n",
        "[_Download the source as a `.jl` file_]($edit_url).\n",
    )
    return replace(content, line => new_title)
end

function _file_list(full_dir, relative_dir, extension)
    return map(
        file -> joinpath(relative_dir, file),
        filter(file -> endswith(file, extension), sort(readdir(full_dir))),
    )
end

"""
    _include_sandbox(filename)

Include the `filename` in a temporary module that acts as a sandbox. (Ensuring
no constants or functions leak into other files.)
"""
function _include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

function _literate_directory(dir)
    for filename in _file_list(dir, dir, ".md")
        if endswith(filename, "introduction.md")
            continue
        elseif endswith(filename, "parallelism.md")
            continue
        else
            rm(filename)
        end
    end
    for filename in _file_list(dir, dir, ".jl")
        # `include` the file to test it before `#src` lines are removed. It is
        # in a testset to isolate local variables between files.
        Test.@testset "$(filename)" begin
            _include_sandbox(filename)
        end
        Literate.markdown(
            filename,
            dir;
            documenter = true,
            postprocess = _link_example,
            # Turn off the footer. We manually add a modified one.
            credit = false,
        )
    end
    return nothing
end

if !_FAST
    for (root, dir, files) in walkdir(joinpath(@__DIR__, "src", "tutorials"))
        _literate_directory.(joinpath.(root, dir))
    end
    # Convert `@example` blocks into `@repl` blocks in the following files:
    for file in [
        joinpath("getting_started", "getting_started_with_julia.md"),
        joinpath("getting_started", "getting_started_with_JuMP.md"),
        joinpath("getting_started", "debugging.md"),
        joinpath("getting_started", "performance_tips.md"),
        joinpath("linear", "tips_and_tricks.md"),
    ]
        filename = joinpath(@__DIR__, "src", "tutorials", file)
        content = read(filename, String)
        content = replace(content, "@example" => "@repl")
        write(filename, content)
    end
end

function _add_edit_url(filename, url)
    contents = read(filename, String)
    open(filename, "w") do io
        write(io, "```@meta\nEditURL = \"$url\"\n```\n\n")
        write(io, contents)
        return
    end
    return
end

# ==============================================================================
#  Add solver README
# ==============================================================================

const _LIST_OF_SOLVERS = Pair{String,String}[]
const _LIST_OF_EXTENSIONS = Pair{String,String}[]
for (solver, data) in TOML.parsefile(joinpath(@__DIR__, "packages.toml"))
    user = get(data, "user", "jump-dev")
    tag = data["rev"]
    filename = get(data, "filename", "README.md")
    out_filename = joinpath(@__DIR__, "src", "packages", "$solver.md")
    Downloads.download(
        "https://raw.githubusercontent.com/$user/$solver.jl/$tag/$filename",
        out_filename,
    )
    _add_edit_url(
        out_filename,
        "https://github.com/$user/$solver.jl/blob/$tag/$filename",
    )
    if get(data, "has_html", false) == true
        # Very simple detector of HTML to wrap in ```@raw html
        lines = readlines(out_filename)
        open(out_filename, "w") do io
            closing_tag = nothing
            math_open = false
            for line in lines
                line = replace(
                    line,
                    "#gh-light-mode-only\"" => "\" class=\"display-light-only\"",
                )
                line = replace(
                    line,
                    "#gh-dark-mode-only\"" => "\" class=\"display-dark-only\"",
                )
                for ext in ("svg", "png", "gif")
                    line = replace(line, ".$ext\"" => ".$ext?raw=true\"")
                end
                if line == "\$\$"
                    line = math_open ? "```" : "```math"
                    math_open = !math_open
                end
                tag = if startswith(line, "<img")
                    "/>"
                else
                    m = match(r"\<([a-z0-9]+)", line)
                    m === nothing ? nothing : "</$(m[1])>"
                end
                if closing_tag === nothing && tag !== nothing
                    println(io, "```@raw html")
                    closing_tag = tag
                end
                println(io, line)
                if closing_tag !== nothing && endswith(line, closing_tag)
                    println(io, "```")
                    closing_tag = nothing
                end
            end
        end
    end
    if get(data, "extension", false)
        push!(_LIST_OF_EXTENSIONS, "$user/$solver.jl" => "packages/$solver.md")
    else
        push!(_LIST_OF_SOLVERS, "$user/$solver.jl" => "packages/$solver.md")
    end
end
push!(
    _LIST_OF_EXTENSIONS,
    "rafaqz/DimensionalData.jl" => "extensions/DimensionalData.md",
)

# Sort, with jump-dev repos at the start.
sort!(_LIST_OF_SOLVERS; by = x -> (!startswith(x[1], "jump-dev/"), x[1]))
sort!(_LIST_OF_EXTENSIONS; by = x -> (!startswith(x[1], "jump-dev/"), x[1]))
pushfirst!(_LIST_OF_SOLVERS, "Introduction" => "packages/solvers.md")
pushfirst!(_LIST_OF_EXTENSIONS, "Introduction" => "extensions/introduction.md")

# ==============================================================================
#  JuMP API
# ==============================================================================

include(joinpath(@__DIR__, "DocumenterReference.jl"))

function sort_by_api_fn((key, type))
    deprecated_methods = [
        :add_nonlinear_constraint,
        :add_nonlinear_expression,
        :add_nonlinear_parameter,
        :all_nonlinear_constraints,
        :get_optimizer_attribute,
        :nonlinear_constraint_string,
        :nonlinear_dual_start_value,
        :nonlinear_expr_string,
        :nonlinear_model,
        :num_nonlinear_constraints,
        :register,
        :set_nonlinear_dual_start_value,
        :set_nonlinear_objective,
        :set_normalized_coefficients,
        :set_optimizer_attribute,
        :set_optimizer_attributes,
        :set_value,
        :NonlinearConstraintIndex,
        :NonlinearConstraintRef,
        :NonlinearExpression,
        :NonlinearParameter,
    ]
    return startswith("$key", "@NL") || key in deprecated_methods
end

jump_api_reference = DocumenterReference.automatic_reference_documentation(;
    root = joinpath(@__DIR__, "src"),
    subdirectory = "api",
    modules = [
        JuMP => [
            "Base.empty!(::GenericModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Base.isempty(::GenericModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Base.copy(::AbstractModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Base.write(::IO, ::GenericModel; ::MOI.FileFormats.FileFormat)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            # Documenter tries to add `read(::IO, ::Type)` instead of this
            # specific method
            # "Base.read(::IO, ::Type{<:GenericModel}; ::MOI.FileFormats.FileFormat)" =>
            #     DocumenterReference.DOCTYPE_FUNCTION,
            "MOI.Utilities.reset_optimizer(::GenericModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "MOI.Utilities.drop_optimizer(::GenericModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "MOI.Utilities.attach_optimizer(::GenericModel)" =>
                DocumenterReference.DOCTYPE_FUNCTION,
        ],
        JuMP.Containers => [
            # TODO(odow): consider exporting these from JuMP.Containers.
            "Containers.@container" => DocumenterReference.DOCTYPE_MACRO,
            "Containers.container" => DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.rowtable" => DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.default_container" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.nested" => DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.vectorized_product" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.build_error_fn" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.parse_macro_arguments" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.parse_ref_sets" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.build_name_expr" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.add_additional_args" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.container_code" =>
                DocumenterReference.DOCTYPE_FUNCTION,
            "Containers.AutoContainerType" =>
                DocumenterReference.DOCTYPE_STRUCT,
            "Containers.NestedIterator" =>
                DocumenterReference.DOCTYPE_STRUCT,
            "Containers.VectorizedProductIterator" =>
                DocumenterReference.DOCTYPE_STRUCT,
        ],
    ],
    sort_by = sort_by_api_fn,
)

# ==============================================================================
#  JuMP documentation structure
# ==============================================================================

# This constant dictates the layout of the documentation. It is manually
# constructed so that we can have control over the order in which pages are
# shown. If you add a new page to the documentation, make sure to add it here!
#
# !!! warning
#     If you move any of the top-level chapters around, make sure to update the
#     index of the "release_notes.md" in the section which builds the PDF.
const _PAGES = [
    "Introduction" => ["index.md", "should_i_use.md", "installation.md"],
    "Tutorials" => [
        "Getting started" => [
            "tutorials/getting_started/introduction.md",
            "tutorials/getting_started/getting_started_with_julia.md",
            "tutorials/getting_started/getting_started_with_JuMP.md",
            "tutorials/getting_started/getting_started_with_sets_and_indexing.md",
            "tutorials/getting_started/getting_started_with_data_and_plotting.md",
            "tutorials/getting_started/debugging.md",
            "tutorials/getting_started/tolerances.md",
            "tutorials/getting_started/design_patterns_for_larger_models.md",
            "tutorials/getting_started/performance_tips.md",
            "tutorials/getting_started/sum_if.md",
        ],
        "Transitioning" =>
            ["tutorials/transitioning/transitioning_from_matlab.md"],
        "Linear programs" => [
            "tutorials/linear/introduction.md",
            "tutorials/linear/knapsack.md",
            "tutorials/linear/diet.md",
            "tutorials/linear/cannery.md",
            "tutorials/linear/factory_schedule.md",
            "tutorials/linear/multi.md",
            "tutorials/linear/multi_commodity_network.md",
            "tutorials/linear/tips_and_tricks.md",
            "tutorials/linear/piecewise_linear.md",
            "tutorials/linear/facility_location.md",
            "tutorials/linear/finance.md",
            "tutorials/linear/geographic_clustering.md",
            "tutorials/linear/network_flows.md",
            "tutorials/linear/transp.md",
            "tutorials/linear/multi_objective_knapsack.md",
            "tutorials/linear/multi_objective_examples.md",
            "tutorials/linear/sudoku.md",
            "tutorials/linear/n-queens.md",
            "tutorials/linear/constraint_programming.md",
            "tutorials/linear/callbacks.md",
            "tutorials/linear/lp_sensitivity.md",
            "tutorials/linear/basis.md",
            "tutorials/linear/mip_duality.md",
        ],
        "Nonlinear programs" => [
            "tutorials/nonlinear/introduction.md",
            "tutorials/nonlinear/simple_examples.md",
            "tutorials/nonlinear/tips_and_tricks.md",
            "tutorials/nonlinear/operator_ad.md",
            "tutorials/nonlinear/user_defined_hessians.md",
            "tutorials/nonlinear/nested_problems.md",
            "tutorials/nonlinear/querying_hessians.md",
            "tutorials/nonlinear/complementarity.md",
            "tutorials/nonlinear/classifiers.md",
            "tutorials/nonlinear/portfolio.md",
            "tutorials/nonlinear/rocket_control.md",
            "tutorials/nonlinear/space_shuttle_reentry_trajectory.md",
        ],
        "Conic programs" => [
            "tutorials/conic/introduction.md",
            "tutorials/conic/tips_and_tricks.md",
            "tutorials/conic/dualization.md",
            "tutorials/conic/arbitrary_precision.md",
            "tutorials/conic/start_values.md",
            "tutorials/conic/simple_examples.md",
            "tutorials/conic/logistic_regression.md",
            "tutorials/conic/experiment_design.md",
            "tutorials/conic/min_ellipse.md",
            "tutorials/conic/ellipse_approx.md",
            "tutorials/conic/quantum_discrimination.md",
        ],
        "Algorithms" => [
            "tutorials/algorithms/benders_decomposition.md",
            "tutorials/algorithms/cutting_stock_column_generation.md",
            "tutorials/algorithms/tsp_lazy_constraints.md",
            "tutorials/algorithms/rolling_horizon.md",
            "tutorials/algorithms/parallelism.md",
        ],
        "Applications" => [
            "tutorials/applications/power_systems.md",
            "tutorials/applications/optimal_power_flow.md",
            "tutorials/applications/web_app.md",
            "tutorials/applications/two_stage_stochastic.md",
        ],
    ],
    "Manual" => [
        "manual/models.md",
        "manual/variables.md",
        "manual/constraints.md",
        "manual/expressions.md",
        "manual/objective.md",
        "manual/containers.md",
        "manual/solutions.md",
        "manual/callbacks.md",
        "manual/complex.md",
        "manual/nonlinear.md",
        "manual/nlp.md",
    ],
    jump_api_reference,
    "Background Information" => [
        "background/algebraic_modeling_languages.md",
        "background/bibliography.md",
    ],
    "Developer Docs" => [
        "Contributing" => "developers/contributing.md",
        "Extensions" => "developers/extensions.md",
        "Custom binaries" => "developers/custom_solver_binaries.md",
        "Style Guide" => "developers/style.md",
        "Roadmap" => "developers/roadmap.md",
        "Checklists" => "developers/checklists.md",
    ],
    "Solvers" => _LIST_OF_SOLVERS,
    "Extensions" => _LIST_OF_EXTENSIONS,
    # "release_notes.md",  # To be added later
]

# ==============================================================================
#  Modify the release notes
# ==============================================================================

function fix_release_line(
    line::String,
    url::String = "https://github.com/jump-dev/JuMP.jl",
)
    # (#XXXX) -> ([#XXXX](url/issue/XXXX))
    while (m = match(r"\(\#([0-9]+)\)", line)) !== nothing
        id = m.captures[1]
        line = replace(line, m.match => "([#$id]($url/issues/$id))")
    end
    # ## Version X.Y.Z -> [Version X.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# Version ([0-9]+.[0-9]+.[0-9]+)", line)) !== nothing
        tag = m.captures[1]
        line = replace(
            line,
            m.match => "## [Version $tag]($url/releases/tag/v$tag)",
        )
    end
    # ## vX.Y.Z -> [vX.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# (v[0-9]+.[0-9]+.[0-9]+)", line)) !== nothing
        tag = m.captures[1]
        line = replace(line, m.match => "## [$tag]($url/releases/tag/$tag)")
    end
    return line
end

function _fix_release_lines(changelog, release_notes, args...)
    open(release_notes, "w") do io
        for line in readlines(changelog; keep = true)
            write(io, fix_release_line(line, args...))
        end
    end
    return
end

_fix_release_lines(
    joinpath(@__DIR__, "src", "changelog.md"),
    joinpath(@__DIR__, "src", "release_notes.md"),
)

_add_edit_url(joinpath(@__DIR__, "src", "release_notes.md"), "changelog.md")

# ==============================================================================
#  Embed MathOptInterface.jl documentation
# ==============================================================================

function _add_moi_pages()
    moi_dir = joinpath(@__DIR__, "src", "moi")
    try
        rm(moi_dir; recursive = true)
    catch
    end
    moi_docs = joinpath(dirname(dirname(pathof(MOI))), "docs")
    cp(joinpath(moi_docs, "src"), moi_dir; force = true)
    # Files in `moi_docs` are probably in read-only mode (`0o444`). Let's give
    # ourselves write permission.
    chmod(moi_dir, 0o777; recursive = true)
    make = read(joinpath(moi_docs, "make.jl"), String)
    # Match from `_PAGES = [` until the start of in `# =====`
    s = strip(match(r"_PAGES = (\[.+?)\#"s, make)[1])
    # Rename every file to the `moi/` directory.
    for m in eachmatch(r"\"([a-zA-Z\_\/]+?\.md)\"", s)
        s = replace(s, m[1] => "moi/" * m[1])
    end
    push!(_PAGES, "MathOptInterface" => eval(Meta.parse(s)))
    # Update the intro of the MOI docs.
    src = """# Introduction

    Welcome to the documentation for MathOptInterface.

    !!! note
        This documentation is also available in PDF format:
        [MathOptInterface.pdf](MathOptInterface.pdf)."""
    status = sprint(io -> Pkg.status("MathOptInterface"; io = io))
    version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
    dest = """# [Introduction](@id moi_documentation)

    !!! warning
        This documentation in this section is a copy of the official
        MathOptInterface documentation available at
        [https://jump.dev/MathOptInterface.jl/$version](https://jump.dev/MathOptInterface.jl/$version).
        It is included here to make it easier to link concepts between JuMP and
        MathOptInterface.
    """
    index_filename = joinpath(moi_dir, "index.md")
    content = replace(read(index_filename, String), src => dest)
    write(index_filename, content)
    _fix_release_lines(
        joinpath(moi_dir, "changelog.md"),
        joinpath(moi_dir, "release_notes.md"),
        "https://github.com/jump-dev/MathOptInterface.jl",
    )
    for (root, dirs, files) in walkdir(moi_dir)
        for file in files
            if endswith(file, ".md")
                filename = joinpath(root, file)
                moi_filename = replace(filename, moi_dir => "")
                _add_edit_url(
                    filename,
                    "https://github.com/jump-dev/MathOptInterface.jl/blob/$version/docs/src$moi_filename",
                )
            end
        end
    end
    return
end

_add_moi_pages()

# ==============================================================================
#  Check that we have included all the markdown files in _PAGES!
# ==============================================================================

# For Documenter.hide arguments
function _add_to_set(
    set::Set{String},
    arg::Tuple{Bool,String,String,Vector{Any}},
)
    _add_to_set(set, arg[3])
    return
end
_add_to_set(set, filename::String) = push!(set, filename)
_add_to_set(set, filename::Pair) = _add_to_set(set, filename[2])
_add_to_set(set, filename::Vector) = _add_to_set.(Ref(set), filename)

function _validate_pages()
    set = Set{String}()
    for page in _PAGES
        _add_to_set(set, page)
    end
    missing_files = String[]
    doc_src = joinpath(@__DIR__, "src", "")
    for (root, dir, files) in walkdir(doc_src)
        for file in files
            if file == "changelog.md" || file == "release_notes.md"
                continue
            end
            filename = replace(joinpath(root, file), doc_src => "")
            if endswith(filename, ".md") && !(filename in set)
                push!(missing_files, filename)
            end
        end
    end
    if !isempty(missing_files)
        error("Some files missing from documentation: $(missing_files)")
    end
    return
end

_validate_pages()

# ==============================================================================
#  Build the HTML docs
# ==============================================================================

Documenter.DocMeta.setdocmeta!(
    JuMP,
    :DocTestSetup,
    :(using JuMP, JuMP.Containers);
    recursive = true,
)

# Needed to make Documenter think that there is a PDF in the right place when
# link checking. Inn production we replace this by running the LaTeX build.
write(joinpath(@__DIR__, "src", "JuMP.pdf"), "")

# Create remotes for Documenter
status = sprint(io -> Pkg.status("MathOptInterface"; io = io))
version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
gh_moi = Documenter.Remotes.GitHub("jump-dev", "MathOptInterface.jl")
remotes = Dict(pkgdir(MOI) => (gh_moi, version))

@time Documenter.makedocs(
    sitename = "JuMP",
    authors = "The JuMP core developers and contributors",
    format = Documenter.HTML(;
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "G-0RZ8X3D3D0",
        mathengine = Documenter.MathJax2(),
        collapselevel = 1,
        assets = ["assets/extra_styles.css", "assets/citations.css"],
        sidebar_sitename = false,
        # Do no check for large pages.
        size_threshold_ignore = [
            "changelog.md",
            "release_notes.md",
            "api/JuMP.md",
            "tutorials/getting_started/getting_started_with_data_and_plotting.md",
            "moi/changelog.md",
            "moi/release_notes.md",
            "moi/reference/models.md",
            "moi/reference/standard_form.md",
            "moi/submodules/Bridges/list_of_bridges.md",
            "moi/submodules/Bridges/reference.md",
            "moi/submodules/Utilities/reference.md",
        ],
    ),
    modules = [JuMP, MOI],
    checkdocs = :none,
    # Skip doctests if --fast provided.
    doctest = _FIX ? :fix : !_FAST,
    pages = vcat(_PAGES, "release_notes.md"),
    remotes = remotes,
    plugins = [
        DocumenterCitations.CitationBibliography(
            joinpath(@__DIR__, "src", "references.bib");
            style = :authoryear,
        ),
    ],
)

# ==============================================================================
#  Build the LaTeX docs (if needed)
# ==============================================================================

function _remove_literate_footer(dir)
    for filename in _file_list(dir, dir, ".md")
        file = read(filename, String)
        index = findfirst(
            "---\n\n!!! tip\n    This tutorial was generated using [Literate",
            file,
        )
        if index !== nothing
            write(filename, file[1:(first(index)-1)])
        end
    end
    return
end

if _PDF
    for (root, dir, files) in walkdir(joinpath(@__DIR__, "src", "tutorials"))
        _remove_literate_footer.(joinpath.(root, dir))
    end
    moi = pop!(_PAGES)   # remove /MathOptInterface
    pop!(moi[2])        # remove /MathOptInterface/release_notes.md
    pop!(_PAGES)        # remove /Extensions
    pop!(_PAGES)        # remove /Solvers
    push!(_PAGES, moi)  # Re-add /MathOptInterface
    section_title, contents = _PAGES[4]
    @assert section_title == "API Reference"
    # `contents` is a big list of docstrings. By default, they'll
    # show up at the `\chapter` level. That's too high.
    _PAGES[4] = section_title => ["Docstrings" => contents]
    latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"
    @time Documenter.makedocs(
        sitename = "JuMP",
        authors = "The JuMP core developers and contributors",
        format = Documenter.LaTeX(; platform = latex_platform),
        build = "latex_build",
        pages = _PAGES,
        debug = true,
        plugins = [
            DocumenterCitations.CitationBibliography(
                joinpath(@__DIR__, "src", "references.bib");
                style = :authoryear,
            ),
        ],
    )
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
    cp(
        joinpath(@__DIR__, "latex_build", "JuMP.pdf"),
        joinpath(@__DIR__, "build", "JuMP.pdf");
        force = true,
    )
end

# ==============================================================================
#  Deploy everything in `build`
# ==============================================================================

Documenter.deploydocs(;
    repo = "github.com/jump-dev/JuMP.jl.git",
    push_preview = true,
)
