#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

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
# oneunit is useful, for example, if coef is a Unitful quantity.
_is_zero_for_printing(coef) = abs(coef) < 1e-10 * oneunit(coef)

# Whether something is one or not for the purposes of printing it.
function _is_one_for_printing(coef)
    return _is_zero_for_printing(abs(coef) - oneunit(coef))
end

function _is_one_for_printing(coef::Complex{T}) where {T}
    r, i = reim(coef)
    return _is_one_for_printing(r) && _is_zero_for_printing(i)
end

function _is_zero_for_printing(coef::Complex)
    return _is_zero_for_printing(real(coef)) &&
           _is_zero_for_printing(imag(coef))
end

_is_im_for_printing(coef) = false

function _is_im_for_printing(coef::Complex)
    r, i = reim(coef)
    return _is_zero_for_printing(r) && _is_one_for_printing(i)
end

_escape_if_scientific(::MIME, x::String) = x

function _escape_if_scientific(::MIME"text/latex", x::String)
    m = match(r"([0-9]+.[0-9]+)e(-?[0-9]+)", x)
    if m === nothing
        return x
    end
    return "$(m[1]) \\times 10^{$(m[2])}"
end

# Helper function that rounds carefully for the purposes of printing Reals
# for example, 5.3  =>  5.3, and 1.0  =>  1
function _string_round(mode, x::Union{Float32,Float64})
    if isinteger(x) && typemin(Int64) <= x <= typemax(Int64)
        return string(round(Int64, x))
    end
    return _escape_if_scientific(mode, string(x))
end

_string_round(mode, ::typeof(abs), x::Real) = _string_round(mode, abs(x))

_sign_string(x::Real) = x < zero(x) ? " - " : " + "

function _string_round(mode, ::typeof(abs), x::Complex)
    r, i = reim(x)
    if _is_zero_for_printing(r)
        return _string_round(mode, Complex(r, abs(i)))
    elseif _is_zero_for_printing(i)
        return _string_round(mode, Complex(abs(r), i))
    else
        return _string_round(mode, x)
    end
end

function _sign_string(x::Complex)
    r, i = reim(x)
    if _is_zero_for_printing(r)
        return _sign_string(i)
    elseif _is_zero_for_printing(i)
        return _sign_string(r)
    else
        return " + "
    end
end

# Fallbacks for other number types

_string_round(mode, x::Any) = string(x)

_string_round(mode, ::typeof(abs), x::Any) = _string_round(mode, x)

_sign_string(::Any) = " + "

function _string_round(mode, x::Complex)
    r, i = reim(x)
    r_str = _string_round(mode, r)
    if _is_zero_for_printing(i)
        return r_str
    elseif _is_zero_for_printing(r)
        if _is_one_for_printing(i)
            if i < 0
                return "-im"
            else
                return "im"
            end
        else
            return string(_string_round(mode, i), "im")
        end
    elseif _is_one_for_printing(i)
        return string("(", r_str, _sign_string(i), "im)")
    else
        abs_i = _string_round(mode, abs, i)
        return string("(", r_str, _sign_string(i), abs_i, "im)")
    end
end

# REPL-specific symbols
# Anything here: https://en.wikipedia.org/wiki/Windows-1252
# should probably work fine on Windows
function _math_symbol(::MIME"text/plain", name::Symbol)
    if name == :leq
        return Sys.iswindows() ? "<=" : "≤"
    elseif name == :geq
        return Sys.iswindows() ? ">=" : "≥"
    elseif name == :eq
        return Sys.iswindows() ? "==" : "="
    elseif name == :sq
        return "²"
    else
        @assert name == :in
        return Sys.iswindows() ? "in" : "∈"
    end
end

# IJulia-specific symbols.
function _math_symbol(::MIME"text/latex", name::Symbol)
    if name == :leq
        return "\\leq"
    elseif name == :geq
        return "\\geq"
    elseif name == :eq
        return "="
    else
        @assert name == :sq
        return "^2"
    end
end

_wrap_in_math_mode(str) = "\$\$ $str \$\$"
_wrap_in_inline_math_mode(str) = "\$ $str \$"

_plural(n) = isone(n) ? "" : "s"

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------

"""
    name(model::AbstractModel)

Return the [`MOI.Name`](@ref) attribute of `model`'s [`backend`](@ref), or a
default if empty.

## Example

```jldoctest
julia> model = Model();

julia> name(model)
"A JuMP Model"
```
"""
name(model::AbstractModel) = "An Abstract JuMP Model"

function name(model::GenericModel)
    if MOI.supports(backend(model), MOI.Name())
        ret = MOI.get(model, MOI.Name())
        return isempty(ret) ? "A JuMP Model" : ret
    end
    return "A JuMP Model"
end

"""
    _print_summary(io::IO, model::AbstractModel)

Print a plain-text summary of `model` to `io`.

For this method to work, an `AbstractModel` subtype should implement:
 * `name(::AbstractModel)`
 * `show_objective_function_summary`
 * `show_constraints_summary`
 * `show_backend_summary`
"""
function _print_summary(io::IO, model::AbstractModel)
    println(io, name(model))
    sense = objective_sense(model)
    if sense == MAX_SENSE
        println(io, "Maximization problem with:")
    elseif sense == MIN_SENSE
        println(io, "Minimization problem with:")
    else
        println(io, "Feasibility problem with:")
    end
    n = num_variables(model)
    println(io, "Variable", _plural(n), ": ", n)
    if sense != FEASIBILITY_SENSE
        show_objective_function_summary(io, model)
    end
    show_constraints_summary(io, model)
    show_backend_summary(io, model)
    if !isempty(object_dictionary(model))
        print(io, "\nNames registered in the model: ")
        print(io, join(sort!(collect(keys(object_dictionary(model)))), ", "))
    end
    return
end

function _print_summary(io::IO, model::GenericModel{T}) where {T}
    println(io, name(model))
    if T != Float64
        println(io, "├ value_type: ", T)
    end
    if mode(model) != AUTOMATIC
        println(io, "├ mode: ", mode(model))
    end
    println(io, "├ solver: ", _try_solver_name(model))
    sense = objective_sense(model)
    println(io, "├ objective_sense: ", sense)
    if sense != FEASIBILITY_SENSE
        F = objective_function_type(model)
        println(io, "│ └ objective_function_type: ", F)
    end
    println(io, "├ num_variables: ", num_variables(model))
    n_constraints = 0
    constraint_lines = String[]
    for (F, S) in list_of_constraint_types(model)
        n = num_constraints(model, F, S)
        n_constraints += n
        line = sprint(MOI.Utilities.print_with_acronym, "$F in $S: $n")
        push!(constraint_lines, line)
    end
    n = num_nonlinear_constraints(model)
    if n > 0
        n_constraints += n
        push!(constraint_lines, "Nonlinear: $n")
    end
    println(io, "├ num_constraints: $n_constraints")
    for (i, line) in enumerate(constraint_lines)
        tag = i == length(constraint_lines) ? "└" : "├"
        println(io, "│ $tag $line")
    end
    print(io, "└ Names registered in the model")
    if isempty(object_dictionary(model))
        print(io, ": none")
    else
        names = sort!(collect(keys(object_dictionary(model))))
        print(io, "\n  └ :", join(names, ", :"))
    end
    return
end

function _try_solver_name(model)
    if mode(model) != DIRECT &&
       MOI.Utilities.state(backend(model)) == MOI.Utilities.NO_OPTIMIZER
        return "none"
    end
    try
        return MOI.get(backend(model), MOI.SolverName())
    catch
        return "unknown"
    end
end

"""
    show_objective_function_summary(io::IO, model::AbstractModel)

Write to `io` a summary of the objective function type.

## Extensions

`AbstractModel`s should implement this method.

## Example

```jldoctest
julia> model = Model();

julia> show_objective_function_summary(stdout, model)
Objective function type: AffExpr
```
"""
function show_objective_function_summary(io::IO, model::GenericModel)
    nlobj = _nlp_objective_function(model)
    print(io, "Objective function type: ")
    if nlobj === nothing
        println(io, objective_function_type(model))
    else
        println(io, "Nonlinear")
    end
    return
end

"""
    show_constraints_summary(io::IO, model::AbstractModel)

Write to `io` a summary of the number of constraints.

## Extensions

`AbstractModel`s should implement this method.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> show_constraints_summary(stdout, model)
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 1 constraint
```
"""
function show_constraints_summary(io::IO, model::GenericModel)
    for (F, S) in list_of_constraint_types(model)
        n = num_constraints(model, F, S)
        println(io, "`$F`-in-`$S`: $n constraint", _plural(n))
    end
    n = num_nonlinear_constraints(model)
    if n > 0
        println(io, "Nonlinear: ", n, " constraint", _plural(n))
    end
    return
end

"""
    show_backend_summary(io::IO, model::GenericModel)

Print a summary of the optimizer backing `model`.

## Extensions

`AbstractModel`s should implement this method.

## Example

```jldoctest
julia> model = Model();

julia> show_backend_summary(stdout, model)
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
```
"""
function show_backend_summary(io::IO, model::GenericModel)
    model_mode = mode(model)
    println(io, "Model mode: ", model_mode)
    if model_mode == MANUAL || model_mode == AUTOMATIC
        println(io, "CachingOptimizer state: ", MOIU.state(backend(model)))
    end
    # The last print shouldn't have a new line
    name = try
        solver_name(model)
    catch
        "unknown"
    end
    if name == "SolverName() attribute not implemented by the optimizer."
        name = "unknown"
    end
    print(io, "Solver name: ", name)
    return
end

"""
    const _CONSTRAINT_LIMIT_FOR_PRINTING = Ref{Int}(100)

A global constant used to control when constraints are omitted when printing
the model.

Get and set this value using `_CONSTRAINT_LIMIT_FOR_PRINTING[]`.

```julia
julia> _CONSTRAINT_LIMIT_FOR_PRINTING[]
100

julia> _CONSTRAINT_LIMIT_FOR_PRINTING[] = 10
10
```
"""
const _CONSTRAINT_LIMIT_FOR_PRINTING = Ref{Int}(100)

"""
    _print_model(io::IO, model::AbstractModel)

Print a plain-text formulation of `model` to `io`.

For this method to work, an `AbstractModel` subtype must implement:
 * `objective_function_string`
 * `constraints_string`
 * `_nl_subexpression_string`
"""
function _print_model(io::IO, model::AbstractModel)
    mode = MIME("text/plain")
    sense = objective_sense(model)
    if sense == MAX_SENSE
        println(io, "Max ", objective_function_string(mode, model))
    elseif sense == MIN_SENSE
        println(io, "Min ", objective_function_string(mode, model))
    else
        println(io, "Feasibility")
    end
    println(io, "Subject to")
    constraints = constraints_string(mode, model)
    m = div(_CONSTRAINT_LIMIT_FOR_PRINTING[], 2)
    skip_start = _CONSTRAINT_LIMIT_FOR_PRINTING[] - m + 1
    skip_stop = length(constraints) - m
    n = skip_stop - skip_start + 1
    for (i, constraint) in enumerate(constraints)
        if i < skip_start || skip_stop < i
            println(io, " ", replace(constraint, '\n' => "\n "))
        end
        if n > 0 && i == skip_start
            println(io, "[[...$n constraints skipped...]]")
        end
    end
    nl_subexpressions = _nl_subexpression_string(mode, model)
    if !isempty(nl_subexpressions)
        println(io, "With NL expressions")
        for expr in nl_subexpressions
            println(io, " ", expr)
        end
    end
    return
end

"""
    _print_latex(io::IO, model::AbstractModel)

Print a LaTeX formulation of `model` to `io`.

For this method to work, an `AbstractModel` subtype must implement:
 * `objective_function_string`
 * `constraints_string`
 * `_nl_subexpression_string`
"""
function _print_latex(io::IO, model::AbstractModel)
    mode = MIME("text/latex")
    println(io, "\$\$ \\begin{aligned}")
    sense = objective_sense(model)
    if sense == MAX_SENSE
        print(io, "\\max\\quad & ")
        println(io, objective_function_string(mode, model), "\\\\")
    elseif sense == MIN_SENSE
        print(io, "\\min\\quad & ")
        println(io, objective_function_string(mode, model), "\\\\")
    else
        println(io, "\\text{feasibility}\\\\")
    end
    constraints = constraints_string(mode, model)
    if !isempty(constraints)
        print(io, "\\text{Subject to} \\quad")
        m = div(_CONSTRAINT_LIMIT_FOR_PRINTING[], 2)
        skip_start = _CONSTRAINT_LIMIT_FOR_PRINTING[] - m + 1
        skip_stop = length(constraints) - m
        n = skip_stop - skip_start + 1
        for (i, constraint) in enumerate(constraints)
            if i < skip_start || skip_stop < i
                println(io, " & ", constraint, "\\\\")
            end
            if n > 0 && i == skip_start
                println(
                    io,
                    " & [[\\ldots\\text{$n constraints skipped}\\ldots]] \\\\",
                )
            end
        end
    end
    nl_subexpressions = _nl_subexpression_string(mode, model)
    if !isempty(nl_subexpressions)
        print(io, "\\text{With NL expressions} \\quad")
        for expr in nl_subexpressions
            println(io, " & ", expr, "\\\\")
        end
    end
    print(io, "\\end{aligned} \$\$")
    return
end

"""
    model_string(mode::MIME, model::AbstractModel)

Return a `String` representation of `model` given the `mode`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> print(model_string(MIME("text/plain"), model))
Feasibility
Subject to
 x ≥ 0
```
"""
function model_string(mode::MIME, model::AbstractModel)
    if mode == MIME("text/latex")
        return sprint(_print_latex, model)
    end
    return sprint(_print_model, model)
end

"""
    objective_function_string(mode, model::AbstractModel)::String

Return a `String` describing the objective function of the model.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @objective(model, Min, 2 * x);

julia> objective_function_string(MIME("text/plain"), model)
"2 x"
```
"""
function objective_function_string(mode, model::GenericModel)
    nlobj = _nlp_objective_function(model)
    if nlobj === nothing
        return function_string(mode, objective_function(model))
    end
    return nonlinear_expr_string(model, mode, nlobj)
end

_set_rhs(s::MOI.LessThan) = :leq, s.upper
_set_rhs(s::MOI.GreaterThan) = :geq, s.lower
_set_rhs(s::MOI.EqualTo) = :eq, s.value
_set_rhs(s::MOI.Interval) = :leq, s.upper
_set_lhs(s::MOI.Interval) = :leq, s.lower
_set_lhs(::Any) = nothing

"""
    nonlinear_constraint_string(
        model::GenericModel,
        mode::MIME,
        c::_NonlinearConstraint,
    )

Return a string representation of the nonlinear constraint `c` belonging to
`model`, given the `mode`.

!!! compat
    This function is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref).
"""
function nonlinear_constraint_string(
    model::GenericModel,
    mode::MIME,
    c::MOI.Nonlinear.ConstraintIndex,
)
    nlp = nonlinear_model(model)::MOI.Nonlinear.Model
    constraint = nlp[c]
    body = nonlinear_expr_string(model, mode, constraint.expression)
    lhs = _set_lhs(constraint.set)
    rhs = _set_rhs(constraint.set)
    output = "$body $(_math_symbol(mode, rhs[1])) $(_string_round(mode, rhs[2]))"
    if lhs === nothing
        return output
    end
    return "$(_string_round(mode, lhs[2])) $(_math_symbol(mode, lhs[1])) $output"
end

"""
    constraints_string(mode, model::AbstractModel)::Vector{String}

Return a list of `String`s describing each constraint of the model.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0);

julia> @constraint(model, c, 2 * x <= 1);

julia> constraints_string(MIME("text/plain"), model)
2-element Vector{String}:
 "c : 2 x ≤ 1"
 "x ≥ 0"
```
"""
function constraints_string(mode, model::GenericModel)
    strings = String[
        constraint_string(mode, cref; in_math_mode = true) for
        (F, S) in list_of_constraint_types(model) for
        cref in all_constraints(model, F, S)
    ]
    for c in all_nonlinear_constraints(model)
        push!(strings, nonlinear_constraint_string(model, mode, index(c)))
    end
    return strings
end

"""
    nonlinear_expr_string(
        model::GenericModel,
        mode::MIME,
        c::MOI.Nonlinear.Expression,
    )

Return a string representation of the nonlinear expression `c` belonging to
`model`, given the `mode`.

!!! compat
    This function is part of the legacy nonlinear interface. Consider using the
    new nonlinear interface documented in [Nonlinear Modeling](@ref).
"""
function nonlinear_expr_string(
    model::GenericModel,
    mode::MIME,
    c::MOI.Nonlinear.Expression,
)
    nlp = nonlinear_model(model)::MOI.Nonlinear.Model
    expr = MOI.Nonlinear.convert_to_expr(nlp, c)
    # Walk terms, and replace
    #    MOI.VariableIndex => VariableRef
    #    MOI.Nonlinear.ExpressionIndex => _NonlinearExpressionIO
    #    MOI.Nonlinear.ParameterIndex => _NonlinearParameterIO
    expr = _replace_expr_terms(model, mode, expr)
    return string(_latexify_exponentials(mode, expr))
end

function _replace_expr_terms(model, mode, expr::Expr)
    for i in 1:length(expr.args)
        expr.args[i] = _replace_expr_terms(model, mode, expr.args[i])
    end
    return expr
end

_replace_expr_terms(::Any, ::Any, expr::Any) = expr

function _replace_expr_terms(model, ::Any, x::MOI.VariableIndex)
    return GenericVariableRef(model, x)
end

# By default, JuMP will print NonlinearExpression objects with some preamble and
# their fill expression. But when printed nested expressions, we only want to
# print `subexpression[i]`. To create this behavior, we create a new object and
# overload `Base.show`, and we replace any ExpressionIndex with this new type.
struct _NonlinearExpressionIO
    model::GenericModel
    mode::MIME
    value::Int
end

function Base.show(io::IO, x::_NonlinearExpressionIO)
    if x.mode == MIME("text/latex")
        return print(io, "subexpression_{$(x.value)}")
    end
    return print(io, "subexpression[$(x.value)]")
end

function _replace_expr_terms(model, mode, x::MOI.Nonlinear.ExpressionIndex)
    return _NonlinearExpressionIO(model, mode, x.value)
end

# We do a similar thing for nonlinear parameters.
struct _NonlinearParameterIO
    model::GenericModel
    mode::MIME
    value::Int
end

function Base.show(io::IO, x::_NonlinearParameterIO)
    for (k, v) in object_dictionary(x.model)
        if v == NonlinearParameter(x.model, x.value)
            return print(io, "$k")
        end
    end
    if x.mode == MIME("text/latex")
        return print(io, "parameter_{$(x.value)}")
    end
    return print(io, "parameter[$(x.value)]")
end

function _replace_expr_terms(model, mode, x::MOI.Nonlinear.ParameterIndex)
    return _NonlinearParameterIO(model, mode, x.value)
end

# Change x ^ -2.0 to x ^ {-2.0}
# x ^ (x ^ 2.0) to x ^ {x ^ {2.0}}
# and so on
_latexify_exponentials(::MIME, ex) = ex

function _latexify_exponentials(mode::MIME"text/latex", ex::Expr)
    for i in 1:length(ex.args)
        ex.args[i] = _latexify_exponentials(mode, ex.args[i])
    end
    if length(ex.args) == 3 && ex.args[1] == :^
        ex.args[3] = Expr(:braces, ex.args[3])
    end
    return ex
end

_nl_subexpression_string(::Any, ::AbstractModel) = String[]

function _nl_subexpression_string(mode::MIME, model::GenericModel)
    nlp_model = nonlinear_model(model)
    strings = String[]
    if nlp_model === nothing
        return strings
    end
    for (k, ex) in enumerate(nlp_model.expressions)
        expr = nonlinear_expr_string(model, mode, ex)
        if mode == MIME("text/latex")
            push!(strings, "subexpression_{$k}: $expr")
        else
            push!(strings, "subexpression[$k]: $expr")
        end
    end
    return strings
end

"""
    anonymous_name(::MIME, x::AbstractVariableRef)

The name to use for an anonymous variable `x` when printing.

## Example

```jldoctest
julia> model = Model();

julia> x = @variable(model);

julia> anonymous_name(MIME("text/plain"), x)
"_[1]"
```
"""
anonymous_name(::Any, x::AbstractVariableRef) = "anon"

function anonymous_name(::MIME"text/plain", x::GenericVariableRef)
    return "_[$(index(x).value)]"
end

function anonymous_name(::MIME"text/latex", x::GenericVariableRef)
    return "{\\_}_{$(index(x).value)}"
end

"""
    function_string(
        mode::MIME,
        func::Union{JuMP.AbstractJuMPScalar,Vector{<:JuMP.AbstractJuMPScalar}},
    )

Return a `String` representing the function `func` using print mode `mode`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> function_string(MIME("text/plain"), 2 * x + 1)
"2 x + 1"
```
"""
function function_string(mode::MIME"text/plain", v::AbstractVariableRef)
    if !is_valid(owner_model(v), v)
        return "InvalidVariableRef"
    end
    var_name = name(v)
    if isempty(var_name)
        return anonymous_name(mode, v)
    end
    return var_name
end

function function_string(mode::MIME"text/latex", v::AbstractVariableRef)
    if !is_valid(owner_model(v), v)
        return "InvalidVariableRef"
    end
    var_name = name(v)
    if isempty(var_name)
        return anonymous_name(mode, v)
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
        return string(m[1]::AbstractString, "_{", m[2]::AbstractString, "}")
    end
    return var_name
end

function _term_string(mode, coef, factor)
    if _is_one_for_printing(coef)
        return factor
    elseif _is_im_for_printing(coef)
        return string(factor, " ", _string_round(mode, abs, coef))
    else
        return string(_string_round(mode, abs, coef), " ", factor)
    end
end

"""
    const _TERM_LIMIT_FOR_PRINTING = Ref{Int}(60)

A global constant used to control when terms are omitted when printing
expressions.

Get and set this value using `_TERM_LIMIT_FOR_PRINTING[]`.

```julia
julia> _TERM_LIMIT_FOR_PRINTING[]
60

julia> _TERM_LIMIT_FOR_PRINTING[] = 10
10
```
"""
const _TERM_LIMIT_FOR_PRINTING = Ref{Int}(60)

_terms_omitted(::MIME, n::Int) = "[[...$n terms omitted...]]"

function _terms_omitted(::MIME"text/latex", n::Int)
    return "[[\\ldots\\text{$n terms omitted}\\ldots]]"
end

function _terms_to_truncated_string(mode, terms)
    m = _TERM_LIMIT_FOR_PRINTING[]
    if length(terms) <= 2 * m
        return join(terms)
    end
    k_l = iseven(m) ? m + 1 : m + 2
    k_r = iseven(m) ? m - 1 : m - 2
    block = _terms_omitted(mode, div(length(terms), 2) - m)
    return string(join(terms[1:k_l]), block, join(terms[(end-k_r):end]))
end

# TODO(odow): remove show_constant in JuMP 1.0
function function_string(mode, a::GenericAffExpr, show_constant = true)
    if length(linear_terms(a)) == 0
        return show_constant ? _string_round(mode, a.constant) : "0"
    end
    terms = fill("", 2 * length(linear_terms(a)))
    for (elm, (coef, var)) in enumerate(linear_terms(a))
        terms[2*elm-1] = _sign_string(coef)
        terms[2*elm] = _term_string(mode, coef, function_string(mode, var))
    end
    terms[1] = terms[1] == " - " ? "-" : ""
    ret = _terms_to_truncated_string(mode, terms)
    if show_constant && !_is_zero_for_printing(a.constant)
        ret = string(
            ret,
            _sign_string(a.constant),
            _string_round(mode, abs, a.constant),
        )
    end
    return ret
end

function function_string(mode, q::GenericQuadExpr)
    if length(quad_terms(q)) == 0
        return function_string(mode, q.aff)
    end
    terms = fill("", 2 * length(quad_terms(q)))
    for (elm, (coef, var1, var2)) in enumerate(quad_terms(q))
        x = function_string(mode, var1)
        y = function_string(mode, var2)
        terms[2*elm-1] = _sign_string(coef)
        if x == y
            factor = x * _math_symbol(mode, :sq)
        else
            times = mode == MIME("text/latex") ? "\\times " : "*"
            factor = string(x, times, y)
        end
        terms[2*elm] = _term_string(mode, coef, factor)
    end
    terms[1] = terms[1] == " - " ? "-" : ""
    ret = _terms_to_truncated_string(mode, terms)
    aff_str = function_string(mode, q.aff)
    if aff_str == "0"
        return ret
    elseif aff_str[1] == '-'
        return string(ret, " - ", aff_str[2:end])
    else
        return string(ret, " + ", aff_str)
    end
end

function function_string(mode, vector::Vector{<:AbstractJuMPScalar})
    strings = function_string.(Ref(mode), vector)
    n = max(div(_TERM_LIMIT_FOR_PRINTING[], 2), 2)
    if length(vector) <= n
        return string("[", join(strings, ", "), "]")
    end
    n_lhs = div(n, 2)
    i_rhs = length(vector) - (n - n_lhs) + 1
    block = _terms_omitted(mode, length(vector) - n)
    return string(
        "[",
        join(vcat(strings[1:n_lhs], block, strings[i_rhs:end]), ", "),
        "]",
    )
end

function function_string(
    ::MIME"text/plain",
    A::AbstractMatrix{<:AbstractJuMPScalar},
)
    str = sprint() do io
        return show(IOContext(io, :limit => true), MIME("text/plain"), A)
    end
    lines = split(str, '\n')
    # We drop the first line with the signature "m×n Array{...}:"
    popfirst!(lines)
    # We replace the first space by an opening `[`
    lines[1] = '[' * lines[1][2:end]
    return join(lines, "\n") * "]"
end

function function_string(
    mode::MIME"text/latex",
    A::AbstractMatrix{<:AbstractJuMPScalar},
)
    str = "\\begin{bmatrix}\n"
    for i in axes(A, 1)
        line = ""
        for j in axes(A, 2)
            if j != 1
                line *= " & "
            end
            if A isa LinearAlgebra.Symmetric && i > j
                line *= "\\cdots"
            elseif A isa LinearAlgebra.UpperTriangular && i > j
                line *= "\\cdots"
            else
                line *= function_string(mode, A[i, j])
            end
        end
        str *= line * "\\\\\n"
    end
    return str * "\\end{bmatrix}"
end

function function_string(mode, constraint::AbstractConstraint)
    f = reshape_vector(jump_function(constraint), shape(constraint))
    return function_string(mode, f)
end

# A special case for symmetric matrix constraints. Since the shape is
# SymmetricMatrixShape, we know that MOI has been passed the upper triangle of
# the matrix. We can make this clearer to users by printing the
# LinearAlgebra.UpperTriangular. There shouldn't be any cases in which the
# constraint function becomes ambiguous.
function function_string(
    mode,
    constraint::VectorConstraint{F,S,SymmetricMatrixShape},
) where {F,S}
    f = reshape_vector(jump_function(constraint), shape(constraint))
    str = function_string(mode, LinearAlgebra.UpperTriangular(f))
    return replace(str, "⋅" => "⋯")
end

function function_string(mode::MIME, p::NonlinearExpression)
    nlp = nonlinear_model(p.model)::MOI.Nonlinear.Model
    expr = nlp[index(p)]
    s = nonlinear_expr_string(p.model, mode, expr)
    return "subexpression[$(p.index)]: " * s
end

function function_string(::MIME, p::NonlinearParameter)
    for (k, v) in object_dictionary(p.model)
        if v == p
            return "$k == $(value(p))"
        end
    end
    return "parameter[$(p.index)] == $(value(p))"
end

"""
    in_set_string(mode::MIME, set)

Return a `String` representing the membership to the set `set` using print mode
`mode`.

## Extensions

JuMP extensions may extend this method for new `set` types to improve the
legibility of their printing.

## Example

```jldoctest
julia> in_set_string(MIME("text/plain"), MOI.Interval(1.0, 2.0))
"∈ [1, 2]"
```
"""
function in_set_string end

function in_set_string(mode::MIME, set::MOI.LessThan)
    return string(_math_symbol(mode, :leq), " ", _string_round(mode, set.upper))
end

function in_set_string(mode::MIME, set::MOI.GreaterThan)
    return string(_math_symbol(mode, :geq), " ", _string_round(mode, set.lower))
end

function in_set_string(mode::MIME, set::MOI.EqualTo)
    return string(_math_symbol(mode, :eq), " ", _string_round(mode, set.value))
end

function in_set_string(::MIME"text/latex", set::MOI.Interval)
    lower = _string_round(mode, set.lower)
    upper = _string_round(mode, set.upper)
    return string("\\in [", lower, ", ", upper, "]")
end

function in_set_string(mode::MIME"text/plain", set::MOI.Interval)
    in = _math_symbol(mode, :in)
    lower = _string_round(mode, set.lower)
    upper = _string_round(mode, set.upper)
    return string("$in [", lower, ", ", upper, "]")
end

in_set_string(::MIME"text/plain", ::MOI.ZeroOne) = "binary"
in_set_string(::MIME"text/latex", ::MOI.ZeroOne) = "\\in \\{0, 1\\}"

in_set_string(::MIME"text/plain", ::MOI.Integer) = "integer"
in_set_string(::MIME"text/latex", ::MOI.Integer) = "\\in \\mathbb{Z}"

function in_set_string(
    mode,
    set::Union{
        PSDCone,
        HermitianPSDCone,
        Zeros,
        Nonnegatives,
        Nonpositives,
        MOI.AbstractSet,
    },
)
    # Use an `if` here instead of multiple dispatch to avoid ambiguity errors.
    if mode == MIME("text/plain")
        return _math_symbol(mode, :in) * " $(set)"
    else
        @assert mode == MIME("text/latex")
        set_str = replace(replace(string(set), "{" => "\\{"), "}" => "\\}")
        return "\\in \\text{$(set_str)}"
    end
end

"""
    in_set_string(mode::MIME, constraint::AbstractConstraint)

Return a `String` representing the membership to the set of the constraint
`constraint` using print mode `mode`.
"""
function in_set_string(mode, constraint::AbstractConstraint)
    # Leave `mode` untyped to avoid ambiguities!
    set = reshape_set(moi_set(constraint), shape(constraint))
    return in_set_string(mode, set)
end

"""
    constraint_string(
        mode::MIME,
        ref::ConstraintRef;
        in_math_mode::Bool = false,
    )

Return a string representation of the constraint `ref`, given the `mode`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, c, 2 * x <= 1);

julia> constraint_string(MIME("text/plain"), c)
"c : 2 x ≤ 1"
```
"""
function constraint_string(mode::MIME, ref::ConstraintRef; in_math_mode = false)
    if !is_valid(owner_model(ref), ref)
        return "InvalidConstraintRef"
    end
    return constraint_string(
        mode,
        name(ref),
        constraint_object(ref);
        in_math_mode = in_math_mode,
    )
end

function constraint_string(mode, constraint_object::AbstractConstraint)
    # Leave `mode` untyped to avoid ambiguities!
    func_str = function_string(mode, constraint_object)
    in_set_str = in_set_string(mode, constraint_object)
    return func_str * " " * in_set_str
end

function constraint_string(
    mode,  # Leave mode untyped to avoid ambiguities
    constraint_name::String,
    constraint_object::AbstractConstraint;
    in_math_mode::Bool = false,
)
    constraint_str = constraint_string(mode, constraint_object)
    if mode == MIME("text/latex")
        # Do not print names in text/latex mode.
        if in_math_mode
            return constraint_str
        else
            return _wrap_in_math_mode(constraint_str)
        end
    end
    prefix = isempty(constraint_name) ? "" : constraint_name * " : "
    return prefix * constraint_str
end

function Base.show(io::IO, ref::ConstraintRef)
    return print(io, constraint_string(MIME("text/plain"), ref))
end

function Base.show(io::IO, ::MIME"text/latex", ref::ConstraintRef)
    return print(io, constraint_string(MIME("text/latex"), ref))
end

function Base.show(io::IO, f::AbstractJuMPScalar)
    return print(io, function_string(MIME("text/plain"), f))
end

function Base.show(io::IO, ::MIME"text/latex", f::AbstractJuMPScalar)
    str = function_string(MIME("text/latex"), f)
    return print(io, _wrap_in_inline_math_mode(str))
end

function Base.show(io::IO, ex::Union{NonlinearExpression,NonlinearParameter})
    return print(io, function_string(MIME("text/plain"), ex))
end

function Base.show(
    io::IO,
    ::MIME"text/latex",
    ex::Union{NonlinearExpression,NonlinearParameter},
)
    return print(io, function_string(MIME("text/latex"), ex))
end

function Base.show(io::IO, c::NonlinearConstraintRef)
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    str = nonlinear_constraint_string(c.model, MIME("text/plain"), index)
    return print(io, str)
end

function Base.show(io::IO, ::MIME"text/latex", c::NonlinearConstraintRef)
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    mode = MIME("text/latex")
    s = _wrap_in_math_mode(nonlinear_constraint_string(c.model, mode, index))
    return print(io, s)
end
