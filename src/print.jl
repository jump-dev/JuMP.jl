#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################
# print.jl
# All "pretty printers" for JuMP types.
# - Delegates to appropriate handler methods for REPL or IJulia.
# - These handler methods then pass the correct symbols to use into a
#   generic string builder. The IJulia handlers will also wrap in MathJax
#   start/close tags.
# - To find printing code for a type in this file, search for `## TypeName`
# - Code here does not need to be fast, in fact simplicity trumps speed
#   within reason as this code is thorny enough as it is.
# - Corresponding tests are in test/print.jl, although test/operator.jl
#   is also testing the constraint/expression code extensively as well.
# - Base.print and Base.string both delegate to Base.show, if they are not
#   separately defined.
#############################################################################

using Printf

# Used for dispatching
"""
    PrintMode

An abstract type used for dispatching printing.
"""
abstract type PrintMode end

"""
    REPLMode

A type used for dispatching printing. Produces plain-text representations.
"""
abstract type REPLMode <: PrintMode end

"""
    IJuliaMode

A type used for dispatching printing. Produces LaTeX representations.
"""
abstract type IJuliaMode <: PrintMode end

Base.show(io::IO, model::AbstractModel) = _print_summary(io, model)

struct _LatexModel{T<:AbstractModel}
    model::T
end

"""
    latex_formulation(model::AbstractModel)

Wrap `model` in a type so that it can be pretty-printed as `text/latex` in a
notebook like IJulia, or in Documenter.

To render the model, end the cell with `latex_formulation(model)`, or call
`display(latex_formulation(model))` in to force the display of the model from
inside a function.
"""
latex_formulation(model::AbstractModel) = _LatexModel(model)

Base.show(io::IO, model::_LatexModel) = _print_latex(io, model.model)
Base.show(io::IO, ::MIME"text/latex", model::_LatexModel) = show(io, model)

function Base.print(model::AbstractModel)
    for d in Base.Multimedia.displays
        if Base.Multimedia.displayable(d, "text/latex") &&
           startswith("$(typeof(d))", "IJulia.")
            return display(d, "text/latex", latex_formulation(model))
        end
    end
    return _print_model(stdout, model)
end
Base.print(io::IO, model::AbstractModel) = _print_model(io, model)

# Whether something is zero or not for the purposes of printing it
# oneunit is useful e.g. if coef is a Unitful quantity.
_is_zero_for_printing(coef) = abs(coef) < 1e-10 * oneunit(coef)
# Whether something is one or not for the purposes of printing it.
_is_one_for_printing(coef) = _is_zero_for_printing(abs(coef) - oneunit(coef))
_sign_string(coef) = coef < zero(coef) ? " - " : " + "

# Helper function that rounds carefully for the purposes of printing
# e.g.   5.3  =>  5.3
#        1.0  =>  1
function _string_round(f::Float64)
    iszero(f) && return "0" # strip sign off zero
    str = string(f)
    return length(str) >= 2 && str[end-1:end] == ".0" ? str[1:end-2] : str
end
_string_round(f) = string(f)

# REPL-specific symbols
# Anything here: https://en.wikipedia.org/wiki/Windows-1252
# should probably work fine on Windows
function _math_symbol(::Type{REPLMode}, name::Symbol)
    if name == :leq
        return Sys.iswindows() ? "<=" : "≤"
    elseif name == :geq
        return Sys.iswindows() ? ">=" : "≥"
    elseif name == :eq
        return Sys.iswindows() ? "==" : "="
    elseif name == :times
        return "*"
    elseif name == :sq
        return "²"
    elseif name == :ind_open
        return "["
    elseif name == :ind_close
        return "]"
    elseif name == :for_all
        return Sys.iswindows() ? "for all" : "∀"
    elseif name == :in
        return Sys.iswindows() ? "in" : "∈"
    elseif name == :open_set
        return "{"
    elseif name == :dots
        return Sys.iswindows() ? ".." : "…"
    elseif name == :close_set
        return "}"
    elseif name == :union
        return Sys.iswindows() ? "or" : "∪"
    elseif name == :infty
        return Sys.iswindows() ? "Inf" : "∞"
    elseif name == :open_rng
        return "["
    elseif name == :close_rng
        return "]"
    elseif name == :integer
        return "integer"
    elseif name == :succeq0
        return " is semidefinite"
    elseif name == :Vert
        return Sys.iswindows() ? "||" : "‖"
    elseif name == :sub2
        return Sys.iswindows() ? "_2" : "₂"
    else
        error("Internal error: Unrecognized symbol $name.")
    end
end

# IJulia-specific symbols.
function _math_symbol(::Type{IJuliaMode}, name::Symbol)
    if name == :leq
        return "\\leq"
    elseif name == :geq
        return "\\geq"
    elseif name == :eq
        return "="
    elseif name == :times
        return "\\times "
    elseif name == :sq
        return "^2"
    elseif name == :ind_open
        return "_{"
    elseif name == :ind_close
        return "}"
    elseif name == :for_all
        return "\\quad\\forall"
    elseif name == :in
        return "\\in"
    elseif name == :open_set
        return "\\{"
    elseif name == :dots
        return "\\dots"
    elseif name == :close_set
        return "\\}"
    elseif name == :union
        return "\\cup"
    elseif name == :infty
        return "\\infty"
    elseif name == :open_rng
        return "\\["
    elseif name == :close_rng
        return "\\]"
    elseif name == :integer
        return "\\in \\mathbb{Z}"
    elseif name == :succeq0
        return "\\succeq 0"
    elseif name == :Vert
        return "\\Vert"
    elseif name == :sub2
        return "_2"
    else
        error("Internal error: Unrecognized symbol $name.")
    end
end

_wrap_in_math_mode(str) = "\$\$ $str \$\$"
_wrap_in_inline_math_mode(str) = "\$ $str \$"

_plural(n) = (isone(n) ? "" : "s")

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------

"""
    _print_summary(io::IO, model::AbstractModel)

Print a plain-text summary of `model` to `io`.

An `AbstractModel` subtype should implement `show_objective_function_summary`,
`show_constraints_summary` and `show_backend_summary` for this method to work.
"""
function _print_summary(io::IO, model::AbstractModel)
    println(io, "A JuMP Model")
    sense = objective_sense(model)
    if sense == MOI.MAX_SENSE
        print(io, "Maximization")
    elseif sense == MOI.MIN_SENSE
        print(io, "Minimization")
    else
        print(io, "Feasibility")
    end
    println(io, " problem with:")
    # TODO: Use MOI.Name for the name of a JuMP model.
    println(
        io,
        "Variable",
        _plural(num_variables(model)),
        ": ",
        num_variables(model),
    )
    if sense != MOI.FEASIBILITY_SENSE
        show_objective_function_summary(io, model)
    end
    show_constraints_summary(io, model)
    show_backend_summary(io, model)
    names_in_scope = sort(collect(keys(object_dictionary(model))))
    if !isempty(names_in_scope)
        println(io)
        print(
            io,
            "Names registered in the model: ",
            join(string.(names_in_scope), ", "),
        )
    end
end

"""
    show_backend_summary(io::IO, model::Model)

Print a summary of the optimizer backing `model`.

`AbstractModel`s should implement this method.
"""
function show_backend_summary(io::IO, model::Model)
    model_mode = mode(model)
    println(io, "Model mode: ", model_mode)
    if model_mode == MANUAL || model_mode == AUTOMATIC
        println(io, "CachingOptimizer state: ", MOIU.state(backend(model)))
    end
    # The last print shouldn't have a new line
    print(io, "Solver name: ", solver_name(model))
    return
end

"""
    _print_model(io::IO, model::AbstractModel)

Print a plain-text formulation of `model` to `io`.

An `AbstractModel` subtype should implement `objective_function_string`,
`constraints_string` and `_nl_subexpression_string` for this method to work.
"""
function _print_model(io::IO, model::AbstractModel)
    sense = objective_sense(model)
    if sense == MOI.MAX_SENSE
        print(io, "Max ")
        println(io, objective_function_string(REPLMode, model))
    elseif sense == MOI.MIN_SENSE
        print(io, "Min ")
        println(io, objective_function_string(REPLMode, model))
    else
        println(io, "Feasibility")
    end
    println(io, "Subject to")
    for constraint in constraints_string(REPLMode, model)
        println(io, " ", replace(constraint, '\n' => "\n "))
    end
    # TODO: Generalize this when similar functionality is needed for
    # `AbstractModel`.
    nl_subexpressions = _nl_subexpression_string(REPLMode, model)
    if !isempty(nl_subexpressions)
        println(io, "With NL expressions")
    end
    for expr in nl_subexpressions
        println(io, " ", expr)
    end
    return
end

"""
    _print_latex(io::IO, model::AbstractModel)

Print a LaTeX formulation of `model` to `io`.

An `AbstractModel` subtype should implement `objective_function_string`,
`constraints_string` and `_nl_subexpression_string` for this method to work.
"""
function _print_latex(io::IO, model::AbstractModel)
    println(io, "\$\$ \\begin{aligned}")
    sense = objective_sense(model)
    if sense == MOI.MAX_SENSE
        print(io, "\\max\\quad & ")
        println(io, objective_function_string(IJuliaMode, model), "\\\\")
    elseif sense == MOI.MIN_SENSE
        print(io, "\\min\\quad & ")
        println(io, objective_function_string(IJuliaMode, model), "\\\\")
    else
        println(io, "\\text{feasibility}\\\\")
    end
    constraints = constraints_string(IJuliaMode, model)
    if !isempty(constraints)
        print(io, "\\text{Subject to} \\quad")
    end
    for constraint in constraints
        println(io, " & ", constraint, "\\\\")
    end
    # TODO: Generalize this when similar functionality is needed for
    # `AbstractModel`.
    nl_subexpressions = _nl_subexpression_string(IJuliaMode, model)
    if !isempty(nl_subexpressions)
        print(io, "\\text{With NL expressions} \\quad")
    end
    for expr in nl_subexpressions
        println(io, " & ", expr, "\\\\")
    end
    return print(io, "\\end{aligned} \$\$")
end

# TODO(odow): deprecate this?
function model_string(print_mode, model::AbstractModel)
    if print_mode == IJuliaMode
        return sprint(_print_latex, model)
    else
        return sprint(_print_model, model)
    end
end

_nl_subexpression_string(print_mode, ::AbstractModel) = String[]

function _nl_subexpression_string(print_mode, model::Model)
    strings = String[]
    if model.nlp_data !== nothing
        num_subexpressions = length(model.nlp_data.nlexpr)::Int
        for k in 1:num_subexpressions
            ex = model.nlp_data.nlexpr[k]
            expr_string = _tape_to_expr(
                model,
                1, # start index in the expression
                ex.nd,
                adjmat(ex.nd),
                ex.const_values,
                [], # parameter_values (not used)
                [], # subexpressions (not needed because !splat_subexpressions)
                model.nlp_data.user_operators,
                false, # generic_variable_names
                false, # splat_subexpressions
                print_mode,
            )
            if print_mode == IJuliaMode
                expr_name = "subexpression_{$k}"
            else
                expr_name = "subexpression[$k]"
            end
            push!(strings, "$expr_name: $expr_string")
        end
    end
    return strings
end

"""
    show_objective_function_summary(io::IO, model::AbstractModel)

Write to `io` a summary of the objective function type.
"""
function show_objective_function_summary(io::IO, model::Model)
    nlobj = _nlp_objective_function(model)
    print(io, "Objective function type: ")
    if nlobj === nothing
        println(io, objective_function_type(model))
    else
        println(io, "Nonlinear")
    end
end

"""
    objective_function_string(print_mode, model::AbstractModel)::String

Return a `String` describing the objective function of the model.
"""
function objective_function_string(print_mode, model::Model)
    nlobj = _nlp_objective_function(model)
    if nlobj === nothing
        return function_string(print_mode, objective_function(model))
    else
        return nl_expr_string(model, print_mode, nlobj)
    end
end

struct _SolutionSummary
    verbose::Bool
    solver::String
    # Status
    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    raw_status::String
    result_count::Int
    has_values::Bool
    has_duals::Bool
    # Candidate solution
    objective_value::Float64
    objective_bound::Union{Missing,Float64}
    dual_objective_value::Union{Missing,Float64}
    primal_solution::Union{Missing,Dict{String,Float64}}
    dual_solution::Union{Missing,Dict{String,Float64}}
    # Work counters
    solve_time::Float64
    barrier_iterations::Union{Missing,Int}
    simplex_iterations::Union{Missing,Int}
    node_count::Union{Missing,Int}
end

"""
    solution_summary(model::Model; verbose::Bool = false)

Return a struct that can be used print a summary of the solution.

If `verbose=true`, write out the primal solution for every variable and the
dual solution for every constraint, excluding those with empty names.

## Examples

When called at the REPL, the summary is automatically printed:
```julia
julia> solution_summary(model)
[...]
```

Use `print` to force the printing of the summary from inside a function:
```julia
function foo(model)
    print(solution_summary(model))
    return
end
```
"""
function solution_summary(model::Model; verbose::Bool = false)
    return _SolutionSummary(
        verbose,
        solver_name(model),
        termination_status(model),
        primal_status(model),
        dual_status(model),
        raw_status(model),
        result_count(model),
        has_values(model),
        has_duals(model),
        objective_value(model),
        _try_get(objective_bound, model),
        _try_get(dual_objective_value, model),
        verbose ? _get_solution_dict(model) : missing,
        verbose ? _get_constraint_dict(model) : missing,
        solve_time(model),
        _try_get(simplex_iterations, model),
        _try_get(barrier_iterations, model),
        _try_get(node_count, model),
    )
end

"""
    Base.show([io::IO], summary::SolutionSummary; verbose::Bool = false)

Write a summary of the solution results to `io` (or to `stdout` if `io` is not
given).
"""
function Base.show(io::IO, summary::_SolutionSummary)
    println(io, "* Solver : ", summary.solver)
    println(io)
    _show_status_summary(io, summary)
    _show_candidate_solution_summary(io, summary)
    _show_work_counters_summary(io, summary)
    return
end

function _show_status_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Status")
    println(io, "  Termination status : ", summary.termination_status)
    println(io, "  Primal status      : ", summary.primal_status)
    println(io, "  Dual status        : ", summary.dual_status)
    if summary.verbose
        println(io, "  Result count       : ", summary.result_count)
        println(io, "  Has duals          : ", summary.has_duals)
    end
    println(io, "  Message from the solver:")
    println(io, "  \"", summary.raw_status, "\"")
    println(io)
    return
end

function _show_candidate_solution_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Candidate solution")
    println(io, "  Objective value      : ", summary.objective_value)
    _print_if_not_missing(
        io,
        "  Objective bound      : ",
        summary.objective_bound,
    )
    _print_if_not_missing(
        io,
        "  Dual objective value : ",
        summary.dual_objective_value,
    )
    if summary.verbose && summary.has_values
        println(io, "  Primal solution :")
        for variable_name in sort(collect(keys(summary.primal_solution)))
            println(
                io,
                "    ",
                variable_name,
                " : ",
                summary.primal_solution[variable_name],
            )
        end
    end
    if summary.verbose && summary.has_duals
        println(io, "  Dual solution :")
        for constraint_name in sort(collect(keys(summary.dual_solution)))
            println(
                io,
                "    ",
                constraint_name,
                " : ",
                summary.dual_solution[constraint_name],
            )
        end
    end
    println(io)
    return
end

function _show_work_counters_summary(io::IO, summary::_SolutionSummary)
    println(io, "* Work counters")
    println(io, "  Solve time (sec)   : ", @sprintf("%.5f", summary.solve_time))
    _print_if_not_missing(
        io,
        "  Simplex iterations : ",
        summary.simplex_iterations,
    )
    _print_if_not_missing(
        io,
        "  Barrier iterations : ",
        summary.barrier_iterations,
    )
    _print_if_not_missing(io, "  Node count         : ", summary.node_count)
    return
end

function _get_solution_dict(model)
    dict = Dict{String,Float64}()
    if has_values(model)
        for x in all_variables(model)
            variable_name = name(x)
            if !isempty(variable_name)
                dict[variable_name] = value(x)
            end
        end
    end
    return dict
end

function _get_constraint_dict(model)
    dict = Dict{String,Float64}()
    if has_duals(model)
        for (F, S) in list_of_constraint_types(model)
            for constraint in all_constraints(model, F, S)
                constraint_name = name(constraint)
                if !isempty(constraint_name)
                    dict[constraint_name] = dual(constraint)
                end
            end
        end
    end
    return dict
end

function _try_get(f, model)
    try
        return f(model)
    catch
        return missing
    end
end

_print_if_not_missing(io, header, ::Missing) = nothing
_print_if_not_missing(io, header, value) = println(io, header, value)

#------------------------------------------------------------------------
## VariableRef
#------------------------------------------------------------------------

function function_string(::Type{REPLMode}, v::AbstractVariableRef)
    var_name = name(v)
    if !isempty(var_name)
        return var_name
    else
        return "noname"
    end
end
function function_string(::Type{IJuliaMode}, v::AbstractVariableRef)
    var_name = name(v)
    if isempty(var_name)
        return "noname"
    end
    # We need to escape latex math characters that appear in the name.
    # However, it's probably impractical to catch everything, so let's just
    # escape the common ones:
    # Escape underscores to prevent them being treated as subscript markers.
    var_name = replace(var_name, "_" => "\\_")
    # Escape carets to prevent them being treated as superscript markers.
    var_name = replace(var_name, "^" => "\\^")
    # Convert any x[args] to x_{args} so that indices on x print as subscripts.
    m = match(r"^(.*)\[(.+)\]$", var_name)
    if m !== nothing
        var_name = m[1] * "_{" * m[2] * "}"
    end

    return var_name
end

#------------------------------------------------------------------------
## GenericAffExpr
#------------------------------------------------------------------------

function function_string(mode, a::GenericAffExpr, show_constant = true)
    # If the expression is empty, return the constant (or 0)
    if length(linear_terms(a)) == 0
        return show_constant ? _string_round(a.constant) : "0"
    end

    term_str = Array{String}(undef, 2 * length(linear_terms(a)))
    elm = 1

    for (coef, var) in linear_terms(a)
        pre = _is_one_for_printing(coef) ? "" : _string_round(abs(coef)) * " "

        term_str[2*elm-1] = _sign_string(coef)
        term_str[2*elm] = string(pre, function_string(mode, var))
        elm += 1
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? _string_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if !_is_zero_for_printing(a.constant) && show_constant
            ret = string(
                ret,
                _sign_string(a.constant),
                _string_round(abs(a.constant)),
            )
        end
        return ret
    end
end

#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------

function function_string(mode, q::GenericQuadExpr)
    length(quad_terms(q)) == 0 && return function_string(mode, q.aff)

    # Odd terms are +/i, even terms are the variables/coeffs
    term_str = Array{String}(undef, 2 * length(quad_terms(q)))
    elm = 1
    if length(term_str) > 0
        for (coef, var1, var2) in quad_terms(q)
            pre =
                _is_one_for_printing(coef) ? "" : _string_round(abs(coef)) * " "

            x = function_string(mode, var1)
            y = function_string(mode, var2)

            term_str[2*elm-1] = _sign_string(coef)
            term_str[2*elm] = "$pre$x"
            if x == y
                term_str[2*elm] *= _math_symbol(mode, :sq)
            else
                term_str[2*elm] *= string(_math_symbol(mode, :times), y)
            end
            if elm == 1
                # Correction for first term as there is no space
                # between - and variable coefficient/name
                term_str[1] = coef < zero(coef) ? "-" : ""
            end
            elm += 1
        end
    end
    ret = join(term_str[1:2*(elm-1)])

    aff_str = function_string(mode, q.aff)
    if aff_str == "0"
        return ret
    else
        if aff_str[1] == '-'
            return string(ret, " - ", aff_str[2:end])
        else
            return string(ret, " + ", aff_str)
        end
    end
end

#------------------------------------------------------------------------
## Constraints
#------------------------------------------------------------------------

"""
    show_constraints_summary(io::IO, model::AbstractModel)

Write to `io` a summary of the number of constraints.
"""
function show_constraints_summary(io::IO, model::Model)
    for (F, S) in list_of_constraint_types(model)
        n_constraints = num_constraints(model, F, S)
        println(
            io,
            "`$F`-in-`$S`: $n_constraints constraint",
            _plural(n_constraints),
        )
    end
    if !iszero(num_nl_constraints(model))
        println(
            io,
            "Nonlinear: ",
            num_nl_constraints(model),
            " constraint",
            _plural(num_nl_constraints(model)),
        )
    end
end

"""
    constraints_string(print_mode, model::AbstractModel)::Vector{String}

Return a list of `String`s describing each constraint of the model.
"""
function constraints_string(print_mode, model::Model)
    strings = String[]
    for (F, S) in list_of_constraint_types(model)
        for cref in all_constraints(model, F, S)
            push!(
                strings,
                constraint_string(print_mode, cref, in_math_mode = true),
            )
        end
    end
    if model.nlp_data !== nothing
        for nl_constraint in model.nlp_data.nlconstr
            push!(
                strings,
                nl_constraint_string(model, print_mode, nl_constraint),
            )
        end
    end
    return strings
end

## Notes for extensions
# For a `ConstraintRef{ModelType, IndexType}` where `ModelType` is not
# `JuMP.Model` or `IndexType` is not `MathOptInterface.ConstraintIndex`, the
# methods `JuMP.name` and `JuMP.constraint_object` should be implemented for
# printing to work. If the `AbstractConstraint` returned by `constraint_object`
# is not `JuMP.ScalarConstraint` nor `JuMP.VectorConstraint`, then either
# `JuMP.jump_function` or `JuMP.function_string` and either `JuMP.moi_set` or
# `JuMP.in_set_string` should be implemented.
function Base.show(io::IO, ref::ConstraintRef)
    return print(io, constraint_string(REPLMode, ref))
end
function Base.show(io::IO, ::MIME"text/latex", ref::ConstraintRef)
    return print(io, constraint_string(IJuliaMode, ref))
end

"""
    function_string(print_mode::Type{<:JuMP.PrintMode},
                    func::Union{JuMP.AbstractJuMPScalar,
                                Vector{<:JuMP.AbstractJuMPScalar}})

Return a `String` representing the function `func` using print mode
`print_mode`.
"""
function function_string end

function Base.show(io::IO, f::AbstractJuMPScalar)
    return print(io, function_string(REPLMode, f))
end
function Base.show(io::IO, ::MIME"text/latex", f::AbstractJuMPScalar)
    return print(io, _wrap_in_math_mode(function_string(IJuliaMode, f)))
end

function function_string(print_mode, vector::Vector{<:AbstractJuMPScalar})
    return "[" * join(function_string.(print_mode, vector), ", ") * "]"
end

function function_string(
    ::Type{REPLMode},
    A::AbstractMatrix{<:AbstractJuMPScalar},
)
    str = sprint(show, MIME"text/plain"(), A)
    lines = split(str, '\n')
    # We drop the first line with the signature "m×n Array{...}:"
    lines = lines[2:end]
    # We replace the first space by an opening `[`
    lines[1] = '[' * lines[1][2:end]
    for i in 1:length(lines)
        lines[i] = lines[i] * (i == length(lines) ? ']' : ';')
    end
    return join(lines, '\n')
end

function function_string(
    print_mode::Type{IJuliaMode},
    A::AbstractMatrix{<:AbstractJuMPScalar},
)
    str = "\\begin{bmatrix}\n"
    for i in 1:size(A, 1)
        line = ""
        for j in 1:size(A, 2)
            if j != 1
                line *= " & "
            end
            if A isa Symmetric && i > j
                line *= "\\cdot"
            else
                line *= function_string(print_mode, A[i, j])
            end
        end
        str *= line * "\\\\\n"
    end
    return str * "\\end{bmatrix}"
end

"""
    function_string(print_mode::{<:JuMP.PrintMode},
                    constraint::JuMP.AbstractConstraint)

Return a `String` representing the function of the constraint `constraint`
using print mode `print_mode`.
"""
function function_string(print_mode, constraint::AbstractConstraint)
    f = reshape_vector(jump_function(constraint), shape(constraint))
    return function_string(print_mode, f)
end

"""
    in_set_string(print_mode::Type{<:PrintMode}, set)

Return a `String` representing the membership to the set `set` using print mode
`print_mode`.
"""
function in_set_string end

function in_set_string(print_mode, set::MOI.LessThan)
    return string(_math_symbol(print_mode, :leq), " ", set.upper)
end

function in_set_string(print_mode, set::MOI.GreaterThan)
    return string(_math_symbol(print_mode, :geq), " ", set.lower)
end

function in_set_string(print_mode, set::MOI.EqualTo)
    return string(_math_symbol(print_mode, :eq), " ", set.value)
end

function in_set_string(print_mode, set::MOI.Interval)
    return string(
        _math_symbol(print_mode, :in),
        " ",
        _math_symbol(print_mode, :open_rng),
        set.lower,
        ", ",
        set.upper,
        _math_symbol(print_mode, :close_rng),
    )
end

in_set_string(print_mode, ::MOI.ZeroOne) = "binary"
in_set_string(print_mode, ::MOI.Integer) = "integer"

in_set_string(::Type{IJuliaMode}, ::MOI.ZeroOne) = "\\in \\{0, 1\\}"
in_set_string(::Type{IJuliaMode}, ::MOI.Integer) = "\\in \\mathbb{Z}"

function in_set_string(print_mode, set::Union{PSDCone,MOI.AbstractSet})
    # Use an `if` here instead of multiple dispatch to avoid ambiguity errors.
    if print_mode == REPLMode
        return _math_symbol(print_mode, :in) * " $(set)"
    else
        set_str = replace(replace(string(set), "{" => "\\{"), "}" => "\\}")
        return "\\in \\text{$(set_str)}"
    end
end

"""
    in_set_string(print_mode::Type{<:PrintMode}, constraint::AbstractConstraint)

Return a `String` representing the membership to the set of the constraint
`constraint` using print mode `print_mode`.
"""
function in_set_string(print_mode, constraint::AbstractConstraint)
    set = reshape_set(moi_set(constraint), shape(constraint))
    return in_set_string(print_mode, set)
end

function constraint_string(print_mode, constraint_object::AbstractConstraint)
    func_str = function_string(print_mode, constraint_object)
    in_set_str = in_set_string(print_mode, constraint_object)
    if print_mode == REPLMode
        lines = split(func_str, '\n')
        lines[1+div(length(lines), 2)] *= " " * in_set_str
        return join(lines, '\n')
    else
        return func_str * " " * in_set_str
    end
end

function constraint_string(
    print_mode,
    constraint_name::String,
    constraint_object::AbstractConstraint;
    in_math_mode::Bool = false,
)
    prefix = isempty(constraint_name) ? "" : constraint_name * " : "
    constraint_str = constraint_string(print_mode, constraint_object)
    if print_mode == IJuliaMode
        if in_math_mode
            return constraint_str
        elseif isempty(prefix)
            return _wrap_in_math_mode(constraint_str)
        else
            return prefix * _wrap_in_inline_math_mode(constraint_str)
        end
    else
        return prefix * constraint_str
    end
end

"""
    constraint_string(
        print_mode,
        ref::ConstraintRef;
        in_math_mode::Bool = false)

Return a string representation of the constraint `ref`, given the `print_mode`.

`print_mode` should be `IJuliaMode` or `REPLMode`.
"""
function constraint_string(print_mode, ref::ConstraintRef; in_math_mode = false)
    return constraint_string(
        print_mode,
        name(ref),
        constraint_object(ref),
        in_math_mode = in_math_mode,
    )
end

#------------------------------------------------------------------------
## _NonlinearExprData
#------------------------------------------------------------------------

"""
    nl_expr_string(model::Model, print_mode, c::_NonlinearExprData)

Return a string representation of the nonlinear expression `c` belonging to
`model`, given the `print_mode`. `print_mode` should be `IJuliaMode` or
`REPLMode`.
"""
function nl_expr_string(model::Model, print_mode, c::_NonlinearExprData)
    ex = _tape_to_expr(
        model,
        1,
        c.nd,
        adjmat(c.nd),
        c.const_values,
        [],
        [],
        model.nlp_data.user_operators,
        false,
        false,
        print_mode,
    )
    if print_mode == IJuliaMode
        ex = _latexify_exponentials(ex)
    end
    return string(ex)
end

# Change x ^ -2.0 to x ^ {-2.0}
# x ^ (x ^ 2.0) to x ^ {x ^ {2.0}}
# and so on
_latexify_exponentials(ex) = ex
function _latexify_exponentials(ex::Expr)
    for i in 1:length(ex.args)
        ex.args[i] = _latexify_exponentials(ex.args[i])
    end
    if length(ex.args) == 3 && ex.args[1] == :^
        ex.args[3] = Expr(:braces, ex.args[3])
    end
    return ex
end
#------------------------------------------------------------------------
## _NonlinearConstraint
#------------------------------------------------------------------------
const NonlinearConstraintRef = ConstraintRef{Model,NonlinearConstraintIndex}

function Base.show(io::IO, c::NonlinearConstraintRef)
    return print(
        io,
        nl_constraint_string(
            c.model,
            REPLMode,
            c.model.nlp_data.nlconstr[c.index.value],
        ),
    )
end

function Base.show(io::IO, ::MIME"text/latex", c::NonlinearConstraintRef)
    constraint = c.model.nlp_data.nlconstr[c.index.value]
    return print(
        io,
        _wrap_in_math_mode(
            nl_constraint_string(c.model, IJuliaMode, constraint),
        ),
    )
end

# TODO: Printing is inconsistent between regular constraints and nonlinear
# constraints because nonlinear constraints don't have names.
"""
    nl_constraint_string(model::Model, mode, c::_NonlinearConstraint)

Return a string representation of the nonlinear constraint `c` belonging to
`model`, given the `print_mode`. `print_mode` should be `IJuliaMode` or
`REPLMode`.
"""
function nl_constraint_string(model::Model, mode, c::_NonlinearConstraint)
    s = _sense(c)
    nl = nl_expr_string(model, mode, c.terms)
    if s == :range
        out_str =
            "$(_string_round(c.lb)) " *
            _math_symbol(mode, :leq) *
            " $nl " *
            _math_symbol(mode, :leq) *
            " " *
            _string_round(c.ub)
    else
        if s == :<=
            rel = _math_symbol(mode, :leq)
        elseif s == :>=
            rel = _math_symbol(mode, :geq)
        else
            rel = _math_symbol(mode, :eq)
        end
        out_str = string(nl, " ", rel, " ", _string_round(_rhs(c)))
    end
    return out_str
end

#------------------------------------------------------------------------
## Opaque nonlinear objects
#------------------------------------------------------------------------
# TODO: This could be pretty printed.

function function_string(::Type{<:PrintMode}, p::NonlinearExpression)
    return "Reference to nonlinear expression #$(p.index)"
end

function function_string(::Type{<:PrintMode}, p::NonlinearParameter)
    relevant_parameters = filter(
        i -> i[2] isa NonlinearParameter && i[2].index == p.index,
        p.m.obj_dict,
    )
    if length(relevant_parameters) == 1
        par_name = first(relevant_parameters)[1]
        return "Reference to nonlinear parameter $(par_name)"
    else
        return "Reference to nonlinear parameter #$(p.index)"
    end
end

function Base.show(io::IO, ex::Union{NonlinearExpression,NonlinearParameter})
    return Base.show(io, function_string(REPLMode, ex))
end

function Base.show(
    io::IO,
    ::MIME"text/latex",
    ex::Union{NonlinearExpression,NonlinearParameter},
)
    return print(io, function_string(IJuliaMode, ex))
end

# TODO: Print the status of the NLPEvaluator, features available, etc.
function Base.show(io::IO, evaluator::NLPEvaluator)
    return Base.show(io, "A JuMP.NLPEvaluator")
end
