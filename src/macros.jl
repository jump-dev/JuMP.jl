#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    op_ifelse(a, x, y)

A function that falls back to `ifelse(a, x, y)`, but when called with a JuMP
variables or expression in the first argument, returns a
[`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_ifelse(true, 1.0, 2.0)
1.0

julia> op_ifelse(x, 1.0, 2.0)
ifelse(x, 1.0, 2.0)

julia> op_ifelse(true, x, 2.0)
x
```
"""
op_ifelse(a, x, y) = ifelse(a, x, y)

# We can't make this a generic `NonlinearOperator` because we only want to
# intercept `ifelse` if the first argument is an `AbstractJuMPScalar` (if it's a
# `Bool`, we want to return the correct branch).
op_ifelse(a::AbstractJuMPScalar, x, y) = NonlinearExpr(:ifelse, Any[a, x, y])

"""
    op_and(x, y)

A function that falls back to `x & y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_and(true, false)
false

julia> op_and(true, x)
true && x
```
"""
const op_and = NonlinearOperator(&, :&&)
# Note that the function is `&` instead of `&&` because `&&` is special lowering
# syntax and is not a regular Julia function, but the MOI operator is `:&&`.

"""
    op_or(x, y)

A function that falls back to `x | y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_or(true, false)
true

julia> op_or(true, x)
true || x
```
"""
const op_or = NonlinearOperator(|, :||)
# Note that the function is `|` instead of `||` because `||` is special lowering
# syntax and is not a regular Julia function, but the MOI operator is `:||`.

"""
    op_strictly_less_than(x, y)

A function that falls back to `x < y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_strictly_less_than(1, 2)
true

julia> op_strictly_less_than(x, 2)
x < 2
```
"""
const op_strictly_less_than = NonlinearOperator(<, :<)

"""
    op_strictly_greater_than(x, y)

A function that falls back to `x > y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_strictly_greater_than(1, 2)
false

julia> op_strictly_greater_than(x, 2)
x > 2
```
"""
const op_strictly_greater_than = NonlinearOperator(>, :>)

"""
    op_less_than_or_equal_to(x, y)

A function that falls back to `x <= y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_less_than_or_equal_to(2, 2)
true

julia> op_less_than_or_equal_to(x, 2)
x <= 2
```
"""
const op_less_than_or_equal_to = NonlinearOperator(<=, :<=)

"""
    op_greater_than_or_equal_to(x, y)

A function that falls back to `x >= y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_greater_than_or_equal_to(2, 2)
true

julia> op_greater_than_or_equal_to(x, 2)
x >= 2
```
"""
const op_greater_than_or_equal_to = NonlinearOperator(>=, :>=)

"""
    op_equal_to(x, y)

A function that falls back to `x == y`, but when called with JuMP variables or
expressions, returns a [`GenericNonlinearExpr`](@ref).

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> op_equal_to(2, 2)
true

julia> op_equal_to(x, 2)
x == 2
```
"""
const op_equal_to = NonlinearOperator(==, :(==))

function _rewrite_to_jump_logic(x)
    if Meta.isexpr(x, :call)
        op = if x.args[1] == :ifelse
            return Expr(:call, op_ifelse, x.args[2:end]...)
        elseif x.args[1] == :<
            return Expr(:call, op_strictly_less_than, x.args[2:end]...)
        elseif x.args[1] == :>
            return Expr(:call, op_strictly_greater_than, x.args[2:end]...)
        elseif x.args[1] == :<=
            return Expr(:call, op_less_than_or_equal_to, x.args[2:end]...)
        elseif x.args[1] == :>=
            return Expr(:call, op_greater_than_or_equal_to, x.args[2:end]...)
        elseif x.args[1] == :(==)
            return Expr(:call, op_equal_to, x.args[2:end]...)
        end
    elseif Meta.isexpr(x, :||)
        # Take special care to ensure short-circuiting behavior of operator is
        # retained. We don't want to evaluate the second argument if the first
        # is `true`.
        @assert length(x.args) == 2
        return Expr(
            :if,
            Expr(:call, ===, x.args[1], true),
            true,
            Expr(:call, op_or, x.args[1], x.args[2]),
        )
    elseif Meta.isexpr(x, :&&)
        # Take special care to ensure short-circuiting behavior of operator is
        # retained. We don't want to evaluate the second argument if the first
        # is `false`.
        @assert length(x.args) == 2
        return Expr(
            :if,
            Expr(:call, ===, x.args[1], false),
            false,
            Expr(:call, op_and, x.args[1], x.args[2]),
        )
    elseif Meta.isexpr(x, :comparison)
        lhs = Expr(:call, x.args[2], x.args[1], x.args[3])
        rhs = Expr(:call, x.args[4], x.args[3], x.args[5])
        return Expr(
            :call,
            op_and,
            _rewrite_to_jump_logic(lhs),
            _rewrite_to_jump_logic(rhs),
        )
    end
    return x
end

"""
    _rewrite_expression(expr)

A helper function so that we can change how we rewrite expressions in a single
place and have it cascade to all locations in the JuMP macros that rewrite
expressions.
"""
function _rewrite_expression(expr::Expr)
    new_expr = MacroTools.postwalk(_rewrite_to_jump_logic, expr)
    new_aff, parse_aff = _MA.rewrite(new_expr; move_factors_into_sums = false)
    ret = gensym()
    code = quote
        $parse_aff
        $ret = $flatten!($new_aff)
    end
    return ret, code
end

"""
    _rewrite_expression(expr)

If `expr` is not an `Expr`, then rewriting it won't do anything. We just need to
copy if it is mutable so that future operations do not modify the user's data.
"""
function _rewrite_expression(expr)
    ret = gensym()
    return ret, :($ret = $_MA.copy_if_mutable($(esc(expr))))
end

"""
    model_convert(
        model::AbstractModel,
        rhs::Union{
            AbstractConstraint,
            Number,
            AbstractJuMPScalar,
            MOI.AbstractSet,
        },
    )

Convert the coefficients and constants of functions and sets in the `rhs` to the
coefficient type `value_type(typeof(model))`.

## Purpose

Creating and adding a constraint is a two-step process. The first step calls
[`build_constraint`](@ref), and the result of that is passed to
[`add_constraint`](@ref).

However, because [`build_constraint`](@ref) does not take the `model` as an
argument, the coefficients and constants of the function or set might be
different than `value_type(typeof(model))`.

Therefore, the result of [`build_constraint`](@ref) is converted in a call to
`model_convert` before the result is passed to [`add_constraint`](@ref).
"""
model_convert(::AbstractModel, rhs::Any) = rhs

function model_convert(model::AbstractModel, set::MOI.AbstractScalarSet)
    if MOI.Utilities.supports_shift_constant(typeof(set))
        T = value_type(typeof(model))
        return MOI.Utilities.shift_constant(set, zero(T))
    end
    return set
end

function model_convert(model::AbstractModel, α::Number)
    T = value_type(typeof(model))
    V = variable_ref_type(model)
    C = _complex_convert_type(T, typeof(α))
    return convert(GenericAffExpr{C,V}, α)
end

function model_convert(model::AbstractModel, con::BridgeableConstraint)
    return BridgeableConstraint(
        model_convert(model, con.constraint),
        con.bridge_type;
        con.coefficient_type,
    )
end

function model_convert(model::AbstractModel, con::ScalarConstraint)
    return ScalarConstraint(
        model_convert(model, con.func),
        model_convert(model, con.set),
    )
end

function model_convert(model::AbstractModel, con::VectorConstraint)
    return VectorConstraint(
        model_convert.(model, con.func),
        model_convert(model, con.set),
        con.shape,
    )
end

function model_convert(model::AbstractModel, x::VariableConstrainedOnCreation)
    return VariableConstrainedOnCreation(
        x.scalar_variable,
        model_convert(model, x.set),
    )
end

function model_convert(
    model::AbstractModel,
    x::AbstractArray{<:VariableConstrainedOnCreation},
)
    return model_convert.(model, x)
end

_valid_model(::AbstractModel, ::Any) = nothing

function _valid_model(m::M, name) where {M}
    return error("Expected $name to be a JuMP model, but it has type $M")
end

"""
    _finalize_macro(
        model,
        code,
        source::LineNumberNode;
        register_name::Union{Nothing,Symbol} = nothing,
        wrap_let::Bool = false,
    )

Wraps the `code` generated by a macro in a code block with the first argument as
`source`, the `LineNumberNode` of where the macro was called from in the user's
code. This results in better stacktraces in error messages.

In addition, this function adds a check that `model` is a valid `AbstractModel`.

If `register_name` is a `Symbol`, register the result of `code` in `model` under
the name `register_name`.

If `wrap_let`, wraps `code` in a `let model = model` block to enforce the model
as a local variable.
"""
function _finalize_macro(
    model::Expr,
    code::Any,
    source::LineNumberNode;
    register_name::Union{Nothing,Symbol} = nothing,
    wrap_let::Bool = false,
)
    @assert Meta.isexpr(model, :escape)
    if wrap_let && model.args[1] isa Symbol
        code = quote
            let $model = $model
                $code
            end
        end
    end
    if register_name !== nothing
        sym_name = Meta.quot(register_name)
        code = quote
            _error_if_cannot_register($model, $sym_name)
            $(esc(register_name)) = $model[$sym_name] = $code
        end
    end
    is_valid_code = :(_valid_model($model, $(Meta.quot(model.args[1]))))
    return Expr(:block, source, is_valid_code, code)
end

function _error_if_cannot_register(model::AbstractModel, name::Symbol)
    obj_dict = object_dictionary(model)
    if haskey(obj_dict, name)
        error(
            """An object of name $name is already attached to this model. If this
          is intended, consider using the anonymous construction syntax, for example,
          `x = @variable(model, [1:N], ...)` where the name of the object does
          not appear inside the macro.

          Alternatively, use `unregister(model, :$(name))` to first unregister
          the existing name from the model. Note that this will not delete the
          object; it will just remove the reference at `model[:$(name)]`.
      """,
        )
    end
    return
end

"""
    _replace_zero(model::M, x) where {M<:AbstractModel}

Replaces `_MA.Zero` with a floating point `zero(value_type(M))`.
"""
_replace_zero(::M, ::_MA.Zero) where {M<:AbstractModel} = zero(value_type(M))

_replace_zero(::AbstractModel, x::Any) = x

function _plural_macro_code(model, block, macro_sym)
    if !Meta.isexpr(block, :block)
        error(
            "Invalid syntax for $(macro_sym)s. The second argument must be a " *
            "`begin end` block. For example:\n" *
            "```julia\n$(macro_sym)s(model, begin\n    # ... lines here ...\nend)\n```.",
        )
    end
    @assert block.args[1] isa LineNumberNode
    last_line = block.args[1]
    code = Expr(:tuple)
    jump_macro = Expr(:., JuMP, QuoteNode(macro_sym))
    for arg in block.args
        if arg isa LineNumberNode
            last_line = arg
        elseif Meta.isexpr(arg, :tuple)  # Line with commas.
            macro_call = Expr(:macrocall, jump_macro, last_line, model)
            # Because of the precedence of "=", Keyword arguments have to appear
            # like: `x, (start = 10, lower_bound = 5)`
            for ex in arg.args
                if Meta.isexpr(ex, :tuple) # embedded tuple
                    append!(macro_call.args, ex.args)
                else
                    push!(macro_call.args, ex)
                end
            end
            push!(code.args, esc(macro_call))
        else  # Stand-alone symbol or expression.
            macro_call = Expr(:macrocall, jump_macro, last_line, model, arg)
            push!(code.args, esc(macro_call))
        end
    end
    return code
end

for file in readdir(joinpath(@__DIR__, "macros"))
    # The check for .jl is necessary because some users may have other files
    # like .cov from running code coverage. See JuMP.jl#3746.
    if endswith(file, ".jl")
        include(joinpath(@__DIR__, "macros", file))
    end
end
