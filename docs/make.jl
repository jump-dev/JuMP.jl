import Documenter
import Literate
import Test

using JuMP
const MathOptInterface = MOI

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
    return content *
           "---\n\n!!! tip\n    This tutorial was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl). [View the source `.jl` file on GitHub]($(edit_url)).\n"
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
        if !endswith(filename, "introduction.md")
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
end

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
            "tutorials/getting_started/performance_tips.md",
            "tutorials/getting_started/design_patterns_for_larger_models.md",
        ],
        "Linear programs" => [
            "tutorials/linear/introduction.md",
            "tutorials/linear/tips_and_tricks.md",
            "tutorials/linear/diet.md",
            "tutorials/linear/cannery.md",
            "tutorials/linear/facility_location.md",
            "tutorials/linear/factory_schedule.md",
            "tutorials/linear/finance.md",
            "tutorials/linear/geographic_clustering.md",
            "tutorials/linear/knapsack.md",
            "tutorials/linear/multi.md",
            "tutorials/linear/n-queens.md",
            "tutorials/linear/network_flows.md",
            "tutorials/linear/prod.md",
            "tutorials/linear/steelT3.md",
            "tutorials/linear/sudoku.md",
            "tutorials/linear/transp.md",
            "tutorials/linear/urban_plan.md",
            "tutorials/linear/callbacks.md",
        ],
        "Nonlinear programs" => [
            "tutorials/nonlinear/introduction.md",
            "tutorials/nonlinear/tips_and_tricks.md",
            "tutorials/nonlinear/portfolio.md",
            "tutorials/nonlinear/qcp.md",
            "tutorials/nonlinear/space_shuttle_reentry_trajectory.md",
            "tutorials/nonlinear/rocket_control.md",
            "tutorials/nonlinear/rosenbrock.md",
            "tutorials/nonlinear/mle.md",
            "tutorials/nonlinear/clnlbeam.md",
            "tutorials/nonlinear/querying_hessians.md",
        ],
        "Conic programs" => [
            "tutorials/conic/introduction.md",
            "tutorials/conic/start_values.md",
            "tutorials/conic/tips_and_tricks.md",
            "tutorials/conic/logistic_regression.md",
            "tutorials/conic/cluster.md",
            "tutorials/conic/corr_sdp.md",
            "tutorials/conic/experiment_design.md",
            "tutorials/conic/max_cut_sdp.md",
            "tutorials/conic/min_distortion.md",
            "tutorials/conic/min_ellipse.md",
            "tutorials/conic/robust_uncertainty.md",
        ],
        "Algorithms" => [
            "tutorials/algorithms/benders_decomposition.md",
            "tutorials/algorithms/cutting_stock_column_generation.md",
            "tutorials/algorithms/tsp_lazy_constraints.md",
        ],
        "Applications" => ["tutorials/applications/power_systems.md"],
    ],
    "Manual" => [
        "manual/models.md",
        "manual/variables.md",
        "manual/constraints.md",
        "manual/expressions.md",
        "manual/objective.md",
        "manual/containers.md",
        "manual/solutions.md",
        "manual/nlp.md",
        "manual/callbacks.md",
    ],
    "API Reference" => [
        "reference/models.md",
        "reference/variables.md",
        "reference/expressions.md",
        "reference/objectives.md",
        "reference/constraints.md",
        "reference/containers.md",
        "reference/solutions.md",
        "reference/nlp.md",
        "reference/callbacks.md",
        "reference/extensions.md",
        "reference/nonlinear.md",
    ],
    "Background Information" =>
        ["background/algebraic_modeling_languages.md"],
    "Developer Docs" => [
        "Contributing" => "developers/contributing.md",
        "Extensions" => "developers/extensions.md",
        "Nonlinear" => "developers/nonlinear.md",
        "Custom binaries" => "developers/custom_solver_binaries.md",
        "Style Guide" => "developers/style.md",
        "Roadmap" => "developers/roadmap.md",
    ],
    "release_notes.md",
]

# ==============================================================================
#  Embed MathOptInterface.jl documentation
# ==============================================================================

function _add_moi_pages()
    moi_docs = joinpath(dirname(dirname(pathof(MOI))), "docs")
    cp(
        joinpath(moi_docs, "src"),
        joinpath(@__DIR__, "src", "moi");
        force = true,
    )
    # Files in `moi_docs` are probably in read-only mode (`0o444`). Let's give
    # ourselves write permission.
    chmod(joinpath(@__DIR__, "src", "moi"), 0o777; recursive = true)
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
    dest = """# [Introduction](@id moi_documentation)

    !!! warning
        This documentation in this section is a copy of the official
        MathOptInterface documentation available at
        [https://jump.dev/MathOptInterface.jl/v1.1.1](https://jump.dev/MathOptInterface.jl/v1.1.1).
        It is included here to make it easier to link concepts between JuMP and
        MathOptInterface.
    """
    index_filename = joinpath(@__DIR__, "src", "moi", "index.md")
    content = replace(read(index_filename, String), src => dest)
    write(index_filename, content)
    return
end

try
    rm(joinpath(@__DIR__, "src", "moi"); recursive = true)
catch
end
_add_moi_pages()

# ==============================================================================
#  Check that we have included all the markdown files in _PAGES!
# ==============================================================================

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

@time Documenter.makedocs(
    sitename = "JuMP",
    authors = "The JuMP core developers and contributors",
    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        analytics = "UA-44252521-1",
        collapselevel = 1,
        assets = ["assets/extra_styles.css"],
        sidebar_sitename = false,
    ),
    strict = true,
    # ==========================================================================
    # `modules = [JuMP]`, along with `checkdocs = :exports` causes Documenter to
    # throw an error if exported functions with docstrings are not contained in
    # the Documentation. However, problematically, we include some MOI docs,
    # which forces us to include MOI in `modules`, despite the fact that we
    # don't necessarily want to document every MOI method.
    #
    # This is should be fine for now, because MOI doesn't export anything.
    # However, also problematically, some doctests in MOI are not checked and
    # are failing. Until they are fixed, we can't enable these options.
    #
    # TODO(odow): uncomment when possible.
    # modules = [JuMP, MOI],
    # checkdocs = :exports,
    # ==========================================================================
    # Skip doctests if --fast provided.
    doctest = _FIX ? :fix : !_FAST,
    pages = _PAGES,
)

# ==============================================================================
#  Build the LaTeX docs (if needed)
# ==============================================================================

if _PDF
    # Remove release notes from PDF
    splice!(_PAGES, 7)   # JuMP release notes
    pop!(_PAGES[end][2]) # MOI release notes
    latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"
    @time Documenter.makedocs(
        sitename = "JuMP",
        authors = "The JuMP core developers and contributors",
        format = Documenter.LaTeX(platform = latex_platform),
        build = "latex_build",
        pages = _PAGES,
    )
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
    cp(
        joinpath(@__DIR__, "latex_build", "JuMP.pdf"),
        joinpath(@__DIR__, "build", "JuMP.pdf"),
    )
end

# ==============================================================================
#  Deploy everything in `build`
# ==============================================================================

Documenter.deploydocs(
    repo = "github.com/jump-dev/JuMP.jl.git",
    push_preview = true,
)
