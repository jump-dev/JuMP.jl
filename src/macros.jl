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
        return Expr(:call, op_or, x.args...)
    elseif Meta.isexpr(x, :&&)
        return Expr(:call, op_and, x.args...)
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
function _rewrite_expression(expr)
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

"""
    _add_kw_args(call, kw_args)

Add the keyword arguments `kw_args` to the function call expression `call`,
escaping the expressions. The elements of `kw_args` should be expressions of the
form `:(key = value)`. The `kw_args` vector can be extracted from the arguments
of a macro with [`Containers._extract_kw_args`](@ref).

## Example

```jldoctest
julia> call = :(f(1, a=2))
:(f(1, a = 2))

julia> JuMP._add_kw_args(call, [:(b=3), :(c=4)])

julia> call
:(f(1, a = 2, \$(Expr(:escape, :(\$(Expr(:kw, :b, 3))))), \$(Expr(:escape, :(\$(Expr(:kw, :c, 4)))))))
```
"""
function _add_kw_args(call, kw_args; exclude = Symbol[])
    for kw in kw_args
        @assert Meta.isexpr(kw, :(=))
        if kw.args[1] in exclude
            continue
        end
        push!(call.args, esc(Expr(:kw, kw.args...)))
    end
    return
end

"""
    _add_positional_args(call, args)::Nothing

Add the positional arguments `args` to the function call expression `call`,
escaping each argument expression. The elements of `args` should be ones that
were extracted via [`Containers._extract_kw_args`](@ref) and had appropriate
arguments filtered out (e.g., the model argument). This is able to incorporate
additional positional arguments to `call`s that already have keyword arguments.

## Example

```jldoctest
julia> call = :(f(1, a=2))
:(f(1, a = 2))

julia> JuMP._add_positional_args(call, [:(x)])

julia> call
:(f(1, $(Expr(:escape, :x)), a = 2))
```
"""
function _add_positional_args(call, args)
    call_args = call.args
    if Meta.isexpr(call, :.)
        # call is broadcasted
        call_args = call.args[2].args
    end
    # Cache all keyword arguments
    kw_args = filter(Base.Fix2(Meta.isexpr, :kw), call_args)
    # Remove keyowrd arguments from the end
    filter!(!Base.Fix2(Meta.isexpr, :kw), call_args)
    # Add the new positional arguments
    append!(call_args, esc.(args))
    # Re-add the cached keyword arguments back to the end
    append!(call_args, kw_args)
    return
end

function _reorder_parameters(args)
    if !Meta.isexpr(args[1], :parameters)
        return args
    end
    args = collect(args)
    p = popfirst!(args)
    for arg in p.args
        @assert arg.head == :kw
        push!(args, Expr(:(=), arg.args[1], arg.args[2]))
    end
    return args
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
    if wrap_let
        code = _wrap_let(model, code)
    end
    if register_name !== nothing
        sym_name = Meta.quot(register_name)
        code = quote
            _error_if_cannot_register($model, $sym_name)
            $(esc(register_name)) = $model[$sym_name] = $code
        end
    end
    return Expr(
        :block,
        source,
        :(_valid_model($model, $(Meta.quot(model.args[1])))),
        code,
    )
end

function _error_if_cannot_register(model::AbstractModel, name::Symbol)
    obj_dict = object_dictionary(model)
    if haskey(obj_dict, name)
        error(
            """An object of name $name is already attached to this model. If this
          is intended, consider using the anonymous construction syntax, e.g.,
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

function _check_vectorized(sense::Symbol)
    sense_str = string(sense)
    if startswith(sense_str, '.')
        return Symbol(sense_str[2:end]), true
    end
    return sense, false
end

"""
    _desparsify(x)

If `x` is an `AbstractSparseArray`, return the dense equivalent, otherwise just
return `x`.

This function is used in `_build_constraint`.

## Why is this needed?

When broadcasting `f.(x)` over an `AbstractSparseArray` `x`, Julia first calls
the equivalent of `f(zero(eltype(x))`. Here's an example:

```jldoctest
julia> import SparseArrays

julia> foo(x) = (println("Calling \$(x)"); x)
foo (generic function with 1 method)

julia> foo.(SparseArrays.sparsevec([1, 2], [1, 2]))
Calling 1
Calling 2
2-element SparseArrays.SparseVector{Int64, Int64} with 2 stored entries:
  [1]  =  1
  [2]  =  2
```

However, if `f` is mutating, this can have serious consequences! In our case,
broadcasting `build_constraint` will add a new `0 = 0` constraint.

Sparse arrays most-often arise when some input data to the constraint is sparse
(e.g., a constant vector or matrix). Due to promotion and arithmetic, this
results in a constraint function that is represented by an `AbstractSparseArray`,
but is actually dense. Thus, we can safely `collect` the matrix into a dense
array.

If the function is sparse, it's not obvious what to do. What is the "zero"
element of the result? What does it mean to broadcast `build_constraint` over a
sparse array adding scalar constraints? This likely means that the user is using
the wrong data structure. For simplicity, let's also call `collect` into a dense
array, and wait for complaints.
"""
_desparsify(x::SparseArrays.AbstractSparseArray) = collect(x)

_desparsify(x) = x

function _functionize(v::V) where {V<:AbstractVariableRef}
    return convert(GenericAffExpr{value_type(V),V}, v)
end

_functionize(v::AbstractArray{<:AbstractVariableRef}) = _functionize.(v)

function _functionize(
    v::LinearAlgebra.Symmetric{V},
) where {V<:AbstractVariableRef}
    return LinearAlgebra.Symmetric(_functionize(v.data))
end

_functionize(x) = x

_functionize(::_MA.Zero) = false

"""
    reverse_sense(::Val{T}) where {T}

Given an (in)equality symbol `T`, return a new `Val` object with the opposite
(in)equality symbol.
"""
function reverse_sense end
reverse_sense(::Val{:<=}) = Val(:>=)
reverse_sense(::Val{:≤}) = Val(:≥)
reverse_sense(::Val{:>=}) = Val(:<=)
reverse_sense(::Val{:≥}) = Val(:≤)
reverse_sense(::Val{:(==)}) = Val(:(==))

# This method is needed because Julia v1.10 prints LineNumberNode in the string
# representation of an expression.
function _strip_LineNumberNode(x::Expr)
    if Meta.isexpr(x, :block)
        return Expr(:block, filter(!Base.Fix2(isa, LineNumberNode), x.args)...)
    end
    return x
end

_strip_LineNumberNode(x) = x

function _macro_error(macroname, args, source, str...)
    str_args = join(_strip_LineNumberNode.(args), ", ")
    return error(
        "At $(source.file):$(source.line): `@$macroname($str_args)`: ",
        str...,
    )
end

# Given a base_name and idxvars, returns an expression that constructs the name
# of the object. For use within macros only.
function _name_call(base_name, idxvars)
    if isempty(idxvars) || base_name == ""
        return base_name
    end
    ex = Expr(:call, :string, base_name, "[")
    for i in 1:length(idxvars)
        # Converting the arguments to strings before concatenating is faster:
        # https://github.com/JuliaLang/julia/issues/29550.
        esc_idxvar = esc(idxvars[i])
        push!(ex.args, :(string($esc_idxvar)))
        i < length(idxvars) && push!(ex.args, ",")
    end
    push!(ex.args, "]")
    return ex
end

_esc_non_constant(x::Number) = x
_esc_non_constant(x::Expr) = Meta.isexpr(x, :quote) ? x : esc(x)
_esc_non_constant(x) = esc(x)

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

function _wrap_let(model, code)
    if Meta.isexpr(model, :escape) && model.args[1] isa Symbol
        return quote
            let $model = $model
                $code
            end
        end
    end
    return code
end

function _get_kwarg_value(kwargs, key::Symbol; default = nothing)
    for kwarg in kwargs
        if kwarg.args[1] == key
            return esc(kwarg.args[2])
        end
    end
    return default
end

include("macros/@objective.jl")
include("macros/@expression.jl")
include("macros/@constraint.jl")
include("macros/@variable.jl")
include("macros/@NL.jl")
