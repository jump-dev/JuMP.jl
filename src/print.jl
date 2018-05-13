#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
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

# Used for dispatching
abstract type PrintMode end
abstract type REPLMode <: PrintMode end
abstract type IJuliaMode <: PrintMode end

# Whether something is zero or not for the purposes of printing it
const PRINT_ZERO_TOL = 1e-10

# List of indices available for variable printing
const DIMS = ["i","j","k","l","m","n"]

# Helper function that rounds carefully for the purposes of printing
# e.g.   5.3  =>  5.3
#        1.0  =>  1
function str_round(f::Float64)
    abs(f) == 0.0 && return "0" # strip sign off zero
    str = string(f)
    length(str) >= 2 && str[end-1:end] == ".0" ? str[1:end-2] : str
end

# TODO: get rid of this! This is only a helper, and should be Base.values
# (and probably live there, as well)
_values(x::Array) = x
_values(x) = Base.values(x)

# REPL-specific symbols
# Anything here: https://en.wikipedia.org/wiki/Windows-1252
# should probably work fine on Windows
const repl = Dict{Symbol,String}(
    :leq        => (Compat.Sys.iswindows() ? "<=" : "≤"),
    :geq        => (Compat.Sys.iswindows() ? ">=" : "≥"),
    :eq         => (Compat.Sys.iswindows() ? "==" : "="),
    :times      => "*",
    :sq         => "²",
    :ind_open   => "[",
    :ind_close  => "]",
    :for_all    => (Compat.Sys.iswindows() ? "for all" : "∀"),
    :in         => (Compat.Sys.iswindows() ? "in" : "∈"),
    :open_set   => "{",
    :dots       => (Compat.Sys.iswindows() ? ".." : "…"),
    :close_set  => "}",
    :union      => (Compat.Sys.iswindows() ? "or" : "∪"),
    :infty      => (Compat.Sys.iswindows() ? "Inf" : "∞"),
    :open_rng   => "[",
    :close_rng  => "]",
    :integer    => "integer",
    :succeq0    => " is semidefinite",
    :Vert       => (Compat.Sys.iswindows() ? "||" : "‖"),
    :sub2       => (Compat.Sys.iswindows() ? "_2" : "₂"))

# IJulia-specific symbols
const ijulia = Dict{Symbol,String}(
    :leq        => "\\leq",
    :geq        => "\\geq",
    :eq         => "=",
    :times      => "\\times ",
    :sq         => "^2",
    :ind_open   => "_{",
    :ind_close  => "}",
    :for_all    => "\\quad\\forall",
    :in         => "\\in",
    :open_set   => "\\{",
    :dots       => "\\dots",
    :close_set  => "\\}",
    :union      => "\\cup",
    :infty      => "\\infty",
    :open_rng   => "\\[",
    :close_rng  => "\\]",
    :integer    => "\\in \\mathbb{Z}",
    :succeq0    => "\\succeq 0",
    :Vert       => "\\Vert",
    :sub2       => "_2")

const PrintSymbols = Dict{Symbol,String}

# If not already mathmode, then wrap in MathJax start/close tags
math(s,mathmode) = mathmode ? s : "\$\$ $s \$\$"

#------------------------------------------------------------------------
## Model
#------------------------------------------------------------------------
function Base.show(io::IO, m::Model) # TODO temporary
    print(io, "A JuMP Model")
end

#------------------------------------------------------------------------
## VariableRef
#------------------------------------------------------------------------
Base.show(io::IO, v::VariableRef) = print(io, var_str(REPLMode,v))
Base.show(io::IO, ::MIME"text/latex", v::VariableRef) =
    print(io, var_str(IJuliaMode,v,mathmode=false))
function var_str(::Type{REPLMode}, v::VariableRef; mathmode=true)
    name = MOI.get(v.m, MOI.VariableName(), v)
    if name != ""
        return name
    else
        return "noname"
    end
end
function var_str(::Type{IJuliaMode}, v::VariableRef; mathmode=true)
    name = MOI.get(v.m, MOI.VariableName(), v)
    if name != ""
        # TODO: This is wrong if variable name constains extra "]"
        return math(replace(replace(name,"[","_{",1),"]","}"), mathmode)
    else
        return math("noname", mathmode)
    end
end


#------------------------------------------------------------------------
## AffExpr  (not GenericAffExpr)
#------------------------------------------------------------------------
Base.show(io::IO, a::AffExpr) = print(io, aff_str(REPLMode,a))
Base.show(io::IO, ::MIME"text/latex", a::AffExpr) =
    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::AffExpr, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(linearterms(a)) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    term_str = Array{String}(undef,2*length(linearterms(a)))
    elm = 1
    # For each non-zero for this model
    for (coef, var) in linearterms(a)
        abs(coef) < PRINT_ZERO_TOL && continue  # e.g. x - x

        pre = abs(abs(coef)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(coef)) * " "

        term_str[2*elm-1] = coef < 0 ? " - " : " + "
        term_str[2*elm  ] = string(pre, var_str(mode, var))
        elm += 1
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? str_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if abs(a.constant) >= PRINT_ZERO_TOL && show_constant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end
# Precompile for faster boot times
Base.precompile(aff_str, (Type{JuMP.REPLMode}, AffExpr, Bool))
Base.precompile(aff_str, (Type{JuMP.IJuliaMode}, AffExpr, Bool))
Base.precompile(aff_str, (Type{JuMP.REPLMode}, AffExpr))
Base.precompile(aff_str, (Type{JuMP.IJuliaMode}, AffExpr))


#------------------------------------------------------------------------
## GenericQuadExpr
#------------------------------------------------------------------------
Base.show(io::IO, q::GenericQuadExpr) = print(io, quad_str(REPLMode,q))
Base.show(io::IO, ::MIME"text/latex", q::GenericQuadExpr) =
    print(io, quad_str(IJuliaMode,q,mathmode=false))
# Generic string converter, called by mode-specific handlers
function quad_str(mode, q::GenericQuadExpr, sym)
    length(quadterms(q)) == 0 && return aff_str(mode,q.aff)

    # Odd terms are +/i, even terms are the variables/coeffs
    term_str = Array{String}(undef,2*length(quadterms(q)))
    elm = 1
    if length(term_str) > 0
        for (coef, var1, var2) in quadterms(q)
            abs(coef) < PRINT_ZERO_TOL && continue  # e.g. x - x

            pre = abs(abs(coef)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(coef)) * " "

            x = var_str(mode,var1)
            y = var_str(mode,var2)

            term_str[2*elm-1] = coef < 0 ? " - " : " + "
            term_str[2*elm  ] = "$pre$x" * (x == y ? sym[:sq] : "$(sym[:times])$y")
            if elm == 1
                # Correction for first term as there is no space
                # between - and variable coefficient/name
                term_str[1] = coef < 0 ? "-" : ""
            end
            elm += 1
        end
    end
    ret = join(term_str[1:2*(elm-1)])

    aff_string = aff_str(mode, q.aff)
    if aff_string == "0"
        return ret
    else
        if aff_string[1] == '-'
            return string(ret, " - ", aff_string[2:end])
        else
            return string(ret, " + ", aff_string)
        end
    end
end

# Handlers to use correct symbols
quad_str(::Type{REPLMode}, q::GenericQuadExpr) =
    quad_str(REPLMode, q, repl)
quad_str(::Type{IJuliaMode}, q::GenericQuadExpr; mathmode=true) =
    math(quad_str(IJuliaMode, q, ijulia), mathmode)

#------------------------------------------------------------------------
## NonlinearExprData
#------------------------------------------------------------------------
#Base.show(io::IO, c::NonlinearExprData) = print(io, expr_str(REPLMode, c))
#Base.show(io::IO, ::MIME"text/latex", c::NonlinearExprData) =
#    print(io, expr_str(IJuliaMode, c))
function expr_str(m::Model, mode, c::NonlinearExprData)
    return string(tapeToExpr(m, 1, c.nd, adjmat(c.nd), c.const_values, [], [], m.nlpdata.user_operators, false, false, mode))
end

# TODO: Print SingleVariableConstraint, VectorOfVariablesConstraint, AffExprConstraint, VectorAffExprConstraint, QuadExprConstraint


#------------------------------------------------------------------------
## NonlinearConstraint
#------------------------------------------------------------------------
Base.show(io::IO, c::NonlinearConstraint) = print(io, con_str(REPLMode,c))
Base.show(io::IO, ::MIME"text/latex", c::NonlinearConstraint) =
    print(io, con_str(IJuliaMode,c,mathmode=false))
# Generic string converter, called by mode-specific handlers
function con_str(m::Model, mode, c::NonlinearConstraint, sym)
    s = sense(c)
    nl = expr_str(m, mode, c.terms)
    if s == :range
        out_str = "$(str_round(c.lb)) $(sym[:leq]) $nl $(sym[:leq]) $(str_round(c.ub))"
    else
        rel = s == :<= ? sym[:leq] : (s == :>= ? sym[:geq] : sym[:eq])
        out_str = string(nl," ",rel," ",str_round(rhs(c)))
    end
    out_str
end
# Handlers to use correct symbols
# con_str(m::Model, ::Type{REPLMode}, c::GenericRangeConstraint; args...) =
#     con_str(m, REPLMode, c, repl)
# con_str(m::Model, ::Type{IJuliaMode}, c::GenericRangeConstraint; mathmode=true) =
#     math(con_str(m, IJuliaMode, c, ijulia), mathmode)

# TODO: Print ConstraintRef

#------------------------------------------------------------------------
## Nonlinear expression/parameter reference
#------------------------------------------------------------------------
Base.show(io::IO, ex::NonlinearExpression) = Base.show(io, "Reference to nonlinear expression #$(ex.index)")
Base.show(io::IO, p::NonlinearParameter) = Base.show(io, "Reference to nonlinear parameter #$(p.index)")
