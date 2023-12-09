#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @constraint(model::GenericModel, expr, kwargs...)

Add a constraint described by the expression `expr`.

    @constraint(model::GenericModel, ref[i=..., j=..., ...], expr, kwargs...)

Add a group of constraints described by the expression `expr` parametrized by
`i`, `j`, ...

The expression `expr` can either be

* of the form `func in set` constraining the function `func` to belong to the
  set `set` which is either a [`MOI.AbstractSet`](@ref)
  or one of the JuMP shortcuts [`SecondOrderCone`](@ref),
  [`RotatedSecondOrderCone`](@ref) and [`PSDCone`](@ref), e.g.
  `@constraint(model, [1, x-1, y-2] in SecondOrderCone())` constrains the norm
  of `[x-1, y-2]` be less than 1;

* of the form `a sign b`, where `sign` is one of `==`, `≥`, `>=`, `≤` and
  `<=` building the single constraint enforcing the comparison to hold for the
  expression `a` and `b`, e.g. `@constraint(model, x^2 + y^2 == 1)` constrains
  `x` and `y` to lie on the unit circle;

* of the form `a ≤ b ≤ c` or `a ≥ b ≥ c` (where `≤` and `<=` (resp. `≥` and
  `>=`) can be used interchangeably) constraining the paired the expression
  `b` to lie between `a` and `c`;

* of the forms `@constraint(m, a .sign b)` or
  `@constraint(m, a .sign b .sign c)` which broadcast the constraint creation to
  each element of the vectors.

The recognized keyword arguments in `kwargs` are the following:

* `base_name`: Sets the name prefix used to generate constraint names. It
  corresponds to the constraint name for scalar constraints, otherwise, the
  constraint names are set to `base_name[...]` for each index `...` of the axes
  `axes`.

* `container`: Specify the container type.

* `set_string_name::Bool = true`: control whether to set the
  [`MOI.ConstraintName`](@ref) attribute. Passing `set_string_name = false` can
  improve performance.

## Note for extending the constraint macro

Each constraint will be created using
`add_constraint(m, build_constraint(error_fn, func, set))` where

* `error_fn` is an error function showing the constraint call in addition to the
  error message given as argument,

* `func` is the expression that is constrained
* and `set` is the set in which it is constrained to belong.

For `expr` of the first type (i.e. `@constraint(m, func in set)`), `func` and
`set` are passed unchanged to `build_constraint` but for the other types, they
are determined from the expressions and signs. For instance,
`@constraint(m, x^2 + y^2 == 1)` is transformed into
`add_constraint(m, build_constraint(error_fn, x^2 + y^2, MOI.EqualTo(1.0)))`.

To extend JuMP to accept new constraints of this form, it is necessary to add
the corresponding methods to `build_constraint`. Note that this will likely mean
that either `func` or `set` will be some custom type, rather than e.g. a
`Symbol`, since we will likely want to dispatch on the type of the function or
set appearing in the constraint.

For extensions that need to create constraints with more information than just
`func` and `set`, an additional positional argument can be specified to
`@constraint` that will then be passed on `build_constraint`. Hence, we can
enable this syntax by defining extensions of
`build_constraint(error_fn, func, set, my_arg; kwargs...)`. This produces the
user syntax: `@constraint(model, ref[...], expr, my_arg, kwargs...)`.
"""
macro constraint(args...)
    return _constraint_macro(args, :constraint, parse_constraint, __source__)
end

"""
    @constraints(model, args...)

Adds groups of constraints at once, in the same fashion as the
[`@constraint`](@ref) macro.

The model must be the first argument, and multiple constraints can be added on
multiple lines wrapped in a `begin ... end` block.

The macro returns a tuple containing the constraints that were defined.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, w);

julia> @variable(model, x);

julia> @variable(model, y);

julia> @variable(model, z[1:3]);

julia> @constraints(model, begin
           x >= 1
           y - w <= 2
           sum_to_one[i=1:3], z[i] + y == 1
       end);

julia> print(model)
Feasibility
Subject to
 sum_to_one[1] : y + z[1] = 1
 sum_to_one[2] : y + z[2] = 1
 sum_to_one[3] : y + z[3] = 1
 x ≥ 1
 -w + y ≤ 2
```
"""
macro constraints(model, block)
    return _plural_macro_code(model, block, Symbol("@constraint"))
end

"""
    @build_constraint(constraint_expr)

Constructs a `ScalarConstraint` or `VectorConstraint` using the same
machinery as [`@constraint`](@ref) but without adding the constraint to a model.

Constraints using broadcast operators like `x .<= 1` are also supported and will
create arrays of `ScalarConstraint` or `VectorConstraint`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @build_constraint(2x >= 1)
ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}(2 x, MathOptInterface.GreaterThan{Float64}(1.0))
```
"""
macro build_constraint(arg)
    function error_fn(str...)
        return _macro_error(:build_constraint, (arg,), __source__, str...)
    end
    if arg isa Symbol
        error_fn(
            "Incomplete constraint specification $arg. " *
            "Are you missing a comparison (<=, >=, or ==)?",
        )
    end
    _, parse_code, build_call = parse_constraint(error_fn, arg)
    return quote
        $parse_code
        $build_call
    end
end

"""
    _constraint_macro(
        args,
        macro_name::Symbol,
        parse_fn::Function,
        source::LineNumberNode,
    )

Returns the code for the macro `@constraint args...` of syntax
```julia
@constraint(model, con, extra_arg, kwargs...)      # single constraint
@constraint(model, ref, con, extra_arg, kwargs...) # group of constraints
```

The expression `con` is parsed by `parse_fn` which returns a `build_constraint`
call code that, when executed, returns an `AbstractConstraint`. The macro
keyword arguments (except the `container` keyword argument which is used to
determine the container type) are added to the `build_constraint` call. The
`extra_arg` is added as terminal positional argument to the `build_constraint`
call along with any keyword arguments (apart from `container` and `base_name`).
The returned value of this call is passed to `add_constraint` which returns a
constraint reference.

`source` is a `LineNumberNode` that should refer to the line that the macro was
called from in the user's code. One way of generating this is via the hidden
variable `__source__`.
"""
function _constraint_macro(
    input_args,
    macro_name::Symbol,
    parse_fn::Function,
    source::LineNumberNode,
)
    error_fn(str...) = _macro_error(macro_name, input_args, source, str...)
    args, kwargs, container = Containers._extract_kw_args(input_args)
    if length(args) < 2 && !isempty(kwargs)
        error_fn(
            "No constraint expression detected. If you are trying to " *
            "construct an equality constraint, use `==` instead of `=`.",
        )
    elseif length(args) < 2
        error_fn("Not enough arguments")
    elseif Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@$(macro_name)s`?")
    end
    model, y, extra = esc(args[1]), args[2], args[3:end]
    # Determine if a reference/container argument was given by the user
    # There are six cases to consider:
    # y                                  | type of y | y.head
    # -----------------------------------+-----------+------------
    # name                               | Symbol    | NA
    # name[1:2]                          | Expr      | :ref
    # name[i = 1:2, j = 1:2; i + j >= 3] | Expr      | :typed_vcat
    # [1:2]                              | Expr      | :vect
    # [i = 1:2, j = 1:2; i + j >= 3]     | Expr      | :vcat
    # a constraint expression            | Expr      | :call or :comparison
    c, x = if y isa Symbol || Meta.isexpr(y, (:vect, :vcat, :ref, :typed_vcat))
        if length(extra) == 0
            error_fn("No constraint expression was given.")
        end
        y, popfirst!(extra)
    else
        nothing, y
    end
    if length(extra) > 1
        error_fn("Cannot specify more than 1 additional positional argument.")
    end
    index_vars, indices = Containers.build_ref_sets(error_fn, c)
    if args[1] in index_vars
        error_fn(
            "Index $(args[1]) is the same symbol as the model. Use a " *
            "different name for the index.",
        )
    end
    is_vectorized, parse_code, build_call = parse_fn(error_fn, x)
    _add_positional_args(build_call, extra)
    _add_kw_args(build_call, kwargs; exclude = [:base_name, :set_string_name])
    base_name = _get_kwarg_value(
        kwargs,
        :base_name;
        default = string(something(Containers._get_name(c), "")),
    )
    set_name_flag = _get_kwarg_value(
        kwargs,
        :set_string_name;
        default = :(set_string_names_on_creation($model)),
    )
    name_expr = Expr(:if, set_name_flag, _name_call(base_name, index_vars), "")
    code = if is_vectorized
        quote
            $parse_code
            # These broadcast calls need to be nested so that the operators
            # are fused. Some broadcasted errors result if you put them on
            # different lines.
            add_constraint.(
                $model,
                model_convert.($model, $build_call),
                $name_expr,
            )
        end
    else
        quote
            $parse_code
            build = model_convert($model, $build_call)
            add_constraint($model, build, $name_expr)
        end
    end
    return _finalize_macro(
        model,
        Containers.container_code(index_vars, indices, code, container),
        source;
        register_name = Containers._get_name(c),
        wrap_let = true,
    )
end

"""
    parse_constraint(error_fn::Function, expr::Expr)

The entry-point for all constraint-related parsing.

## Arguments

 * The `error_fn` function is passed everywhere to provide better error messages
 * `expr` comes from the `@constraint` macro. There are two possibilities:
    * `@constraint(model, expr)`
    * `@constraint(model, name[args], expr)`
   In both cases, `expr` is the main component of the constraint.

## Supported syntax

JuMP currently supports the following `expr` objects:
 * `lhs <= rhs`
 * `lhs == rhs`
 * `lhs >= rhs`
 * `l <= body <= u`
 * `u >= body >= l`
 * `lhs ⟂ rhs`
 * `lhs in rhs`
 * `lhs ∈ rhs`
 * `z --> {constraint}`
 * `!z --> {constraint}`
 * `z <--> {constraint}`
 * `!z <--> {constraint}`
 * `z => {constraint}`
 * `!z => {constraint}`
as well as all broadcasted variants.

## Extensions

The infrastructure behind `parse_constraint` is extendable. See
[`parse_constraint_head`](@ref) and [`parse_constraint_call`](@ref) for details.
"""
function parse_constraint(error_fn::Function, expr::Expr)
    return parse_constraint_head(error_fn, Val(expr.head), expr.args...)
end

"""
    parse_constraint_head(error_fn::Function, ::Val{head}, args...)

Implement this method to intercept the parsing of an expression with head
`head`.

!!! warning
    Extending the constraint macro at parse time is an advanced operation and
    has the potential to interfere with existing JuMP syntax. Please discuss
    with the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev) before
    publishing any code that implements these methods.

## Arguments

 * `error_fn`: a function that accepts a `String` and throws the string as an
   error, along with some descriptive information of the macro from which it was
   thrown.
 * `head`: the `.head` field of the `Expr` to intercept
 * `args...`: the `.args` field of the `Expr`.

## Returns

This function must return:

 * `is_vectorized::Bool`: whether the expression represents a broadcasted
   expression like `x .<= 1`
 * `parse_code::Expr`: an expression containing any setup or rewriting code that
   needs to be called before `build_constraint`
 * `build_code::Expr`: an expression that calls `build_constraint(` or
   `build_constraint.(` depending on `is_vectorized`.

## Existing implementations

JuMP currently implements:

   * `::Val{:call}`, which forwards calls to [`parse_constraint_call`](@ref)
   * `::Val{:comparison}`, which handles the special case of `l <= body <= u`.

See also: [`parse_constraint_call`](@ref), [`build_constraint`](@ref)
"""
function parse_constraint_head(error_fn::Function, ::Val{T}, args...) where {T}
    return error_fn(
        "Unsupported constraint expression: we don't know how to parse " *
        "constraints containing expressions of type :$T.\n\nIf you are " *
        "writing a JuMP extension, implement " *
        "`parse_constraint_head(::Function, ::Val{:$T}, args...)",
    )
end

function parse_constraint_head(
    error_fn::Function,
    ::Val{:call},
    op::Symbol,
    args...,
)
    op, is_vectorized = _check_vectorized(op)
    parse_code, build_call =
        parse_constraint_call(error_fn, is_vectorized, Val(op), args...)
    return is_vectorized, parse_code, build_call
end

function parse_constraint_head(
    error_fn::Function,
    ::Val{:comparison},
    lb,
    lsign::Symbol,
    aff,
    rsign::Symbol,
    ub,
)
    lsign, lvectorized = _check_vectorized(lsign)
    rsign, rvectorized = _check_vectorized(rsign)
    if lvectorized != rvectorized
        error_fn("Operators are inconsistently vectorized.")
    end
    if lsign in (:(<=), :≤) && rsign in (:(<=), :≤)
        # Nothing. What we expect.
    elseif lsign in (:(>=), :≥) && rsign in (:(>=), :≥)
        # Flip lb and ub
        lb, ub = ub, lb
    else
        error_fn(
            "unsupported mix of comparison operators " *
            "`$lb $lsign ... $rsign $ub`.\n\n" *
            "Two-sided rows must of the form `$lb <= ... <= $ub` or " *
            "`$ub >= ... >= $lb`.",
        )
    end
    new_aff, parse_aff = _rewrite_expression(aff)
    new_lb, parse_lb = _rewrite_expression(lb)
    new_ub, parse_ub = _rewrite_expression(ub)
    parse_code = quote
        $parse_aff
        $parse_lb
        $parse_ub
    end
    build_call = if lvectorized
        :(
            build_constraint.(
                $error_fn,
                _desparsify($new_aff),
                _desparsify($new_lb),
                _desparsify($new_ub),
            )
        )
    else
        :(build_constraint($error_fn, $new_aff, $new_lb, $new_ub))
    end
    return lvectorized, parse_code, build_call
end

"""
    operator_to_set(error_fn::Function, ::Val{sense_symbol})

Converts a sense symbol to a set `set` such that
`@constraint(model, func sense_symbol 0)` is equivalent to
`@constraint(model, func in set)` for any `func::AbstractJuMPScalar`.

## Example

Once a custom set is defined you can directly create a JuMP constraint with it:
```jldoctest operator_to_set
julia> struct CustomSet{T} <: MOI.AbstractScalarSet
           value::T
       end

julia> Base.copy(x::CustomSet) = CustomSet(x.value)

julia> model = Model();

julia> @variable(model, x)
x

julia> cref = @constraint(model, x in CustomSet(1.0))
x ∈ CustomSet{Float64}(1.0)
```

However, there might be an appropriate sign that could be used in order to
provide a more convenient syntax:
```jldoctest operator_to_set
julia> JuMP.operator_to_set(::Function, ::Val{:⊰}) = CustomSet(0.0)

julia> MOIU.supports_shift_constant(::Type{<:CustomSet}) = true

julia> MOIU.shift_constant(set::CustomSet, value) = CustomSet(set.value + value)

julia> cref = @constraint(model, x ⊰ 1)
x ∈ CustomSet{Float64}(1.0)
```
Note that the whole function is first moved to the right-hand side, then the
sign is transformed into a set with zero constant and finally the constant is
moved to the set with `MOIU.shift_constant`.
"""
function operator_to_set(error_fn::Function, ::Val{S}) where {S}
    return error_fn("unsupported operator $S")
end

function operator_to_set(error_fn::Function, ::Val{:>})
    return error_fn(
        "unsupported operator `>`.\n\n" *
        "JuMP does not support strict inequalities, use `>=` instead.\n\n" *
        "If you require a strict inequality, you will need to use a " *
        "tolerance. For example, instead of `x > 1`, do `x >= 1 + 1e-4`. " *
        "If the constraint must take integer values, use a tolerance of " *
        "`1.0`. If the constraint may take continuous values, note that this " *
        "work-around can cause numerical issues, and your constraint may not " *
        "hold exactly.",
    )
end

function operator_to_set(error_fn::Function, ::Val{:<})
    return error_fn(
        "unsupported operator `<`.\n\n" *
        "JuMP does not support strict inequalities, use `<=` instead.\n\n" *
        "If you require a strict inequality, you will need to use a " *
        "tolerance. For example, instead of `x < 1`, do `x <= 1 - 1e-4`. " *
        "If the constraint must take integer values, use a tolerance of " *
        "`1.0`. If the constraint may take continuous values, note that this " *
        "work-around can cause numerical issues, and your constraint may not " *
        "hold exactly.",
    )
end

"""
    Nonnegatives()

The JuMP equivalent of the [`MOI.Nonnegatives`](@ref) set, in which the
dimension is inferred from the corresponding function.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, x in Nonnegatives())
[x[1], x[2]] ∈ MathOptInterface.Nonnegatives(2)

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x >= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Nonnegatives(2)
```
"""
struct Nonnegatives end

operator_to_set(::Function, ::Union{Val{:(>=)},Val{:(≥)}}) = Nonnegatives()

"""
    Nonpositives()

The JuMP equivalent of the [`MOI.Nonpositives`](@ref) set, in which the
dimension is inferred from the corresponding function.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, x in Nonpositives())
[x[1], x[2]] ∈ MathOptInterface.Nonpositives(2)

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x <= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Nonpositives(2)
```
"""
struct Nonpositives end

operator_to_set(::Function, ::Union{Val{:(<=)},Val{:(≤)}}) = Nonpositives()

"""
    Zeros()

The JuMP equivalent of the [`MOI.Zeros`](@ref) set, in which the dimension is
inferred from the corresponding function.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, x in Zeros())
[x[1], x[2]] ∈ MathOptInterface.Zeros(2)

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x == b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Zeros(2)
```
"""
struct Zeros end

operator_to_set(::Function, ::Val{:(==)}) = Zeros()

"""
    parse_constraint_call(
        error_fn::Function,
        is_vectorized::Bool,
        ::Val{op},
        args...,
    )

Implement this method to intercept the parsing of a `:call` expression with
operator `op`.

!!! warning
    Extending the constraint macro at parse time is an advanced operation and
    has the potential to interfere with existing JuMP syntax. Please discuss
    with the [developer chatroom](https://gitter.im/JuliaOpt/jump-dev) before
    publishing any code that implements these methods.

## Arguments

 * `error_fn`: a function that accepts a `String` and throws the string as an
   error, along with some descriptive information of the macro from which it was
   thrown.
 * `is_vectorized`: a boolean to indicate if `op` should be broadcast or not
 * `op`: the first element of the `.args` field of the `Expr` to intercept
 * `args...`: the `.args` field of the `Expr`.

## Returns

This function must return:

 * `parse_code::Expr`: an expression containing any setup or rewriting code that
   needs to be called before `build_constraint`
 * `build_code::Expr`: an expression that calls `build_constraint(` or
   `build_constraint.(` depending on `is_vectorized`.

See also: [`parse_constraint_head`](@ref), [`build_constraint`](@ref)
"""
function parse_constraint_call(
    error_fn::Function,
    ::Bool,
    ::Val{T},
    args...,
) where {T}
    return error_fn(
        "Unsupported constraint expression: we don't know how to parse " *
        "constraints containing the operator $T.\n\nIf you are writing a " *
        "JuMP extension, implement " *
        "`parse_constraint_call(::Function, ::Bool, ::Val{$T}, args...)",
    )
end

# `@constraint(model, func in set)`
# `@constraint(model, func ∈ set)`
function parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    ::Union{Val{:in},Val{:∈}},
    func,
    set,
)
    f, parse_code = _rewrite_expression(func)
    build_call = if vectorized
        :(build_constraint.($error_fn, _desparsify($f), Ref($(esc(set)))))
    else
        :(build_constraint($error_fn, $f, $(esc(set))))
    end
    return parse_code, build_call
end

"""
    parse_constraint_call(
        error_fn::Function,
        vectorized::Bool,
        ::Val{op},
        lhs,
        rhs,
    ) where {op}

Fallback handler for binary operators. These might be infix operators like
`@constraint(model, lhs op rhs)`, or normal operators like
`@constraint(model, op(lhs, rhs))`.

In both cases, we rewrite as `lhs - rhs in operator_to_set(error_fn, op)`.

See [`operator_to_set`](@ref) for details.
"""
function parse_constraint_call(
    error_fn::Function,
    vectorized::Bool,
    operator::Val,
    lhs,
    rhs,
)
    func = vectorized ? :($lhs .- $rhs) : :($lhs - $rhs)
    f, parse_code = _rewrite_expression(func)
    set = operator_to_set(error_fn, operator)
    # `_functionize` deals with the pathological case where the `lhs` is a
    # `VariableRef` and the `rhs` is a summation with no terms.
    f = :(_functionize($f))
    build_call = if vectorized
        :(build_constraint.($error_fn, _desparsify($f), Ref($(esc(set)))))
    else
        :(build_constraint($error_fn, $f, $(esc(set))))
    end
    return parse_code, build_call
end

function build_constraint(
    error_fn::Function,
    f,
    set::Nonnegatives,
    args...;
    kwargs...,
)
    return build_constraint(
        error_fn,
        f,
        MOI.GreaterThan(false),
        args...;
        kwargs...,
    )
end

function build_constraint(
    error_fn::Function,
    f,
    set::Nonpositives,
    args...;
    kwargs...,
)
    return build_constraint(
        error_fn,
        f,
        MOI.LessThan(false),
        args...;
        kwargs...,
    )
end

function build_constraint(error_fn::Function, f, set::Zeros, args...; kwargs...)
    return build_constraint(error_fn, f, MOI.EqualTo(false), args...; kwargs...)
end

function build_constraint(
    error_fn::Function,
    f::AbstractVector,
    set::Nonnegatives,
)
    return build_constraint(error_fn, f, MOI.Nonnegatives(length(f)))
end

function build_constraint(
    error_fn::Function,
    f::AbstractVector,
    set::Nonpositives,
)
    return build_constraint(error_fn, f, MOI.Nonpositives(length(f)))
end

function build_constraint(error_fn::Function, f::AbstractVector, set::Zeros)
    return build_constraint(error_fn, f, MOI.Zeros(length(f)))
end

# Generic fallback.
function build_constraint(error_fn::Function, func, set, args...; kwargs...)
    arg_str = join(args, ", ")
    arg_str = isempty(arg_str) ? "" : ", " * arg_str
    kwarg_str = join(Tuple(string(k, " = ", v) for (k, v) in kwargs), ", ")
    kwarg_str = isempty(kwarg_str) ? "" : "; " * kwarg_str
    return error_fn(
        "Unrecognized constraint building format. Tried to invoke " *
        "`build_constraint(error, $(func), $(set)$(arg_str)$(kwarg_str))`, " *
        "but no such method exists. This is due to specifying an unrecognized " *
        "function, constraint set, and/or extra positional/keyword arguments." *
        "\n\nIf you're trying to create a JuMP extension, you need to " *
        "implement `build_constraint` to accomodate these arguments.",
    )
end

function build_constraint(
    error_fn::Function,
    func,
    ::Union{MOI.AbstractScalarSet,MOI.AbstractVectorSet},
)
    return error_fn(
        "Unable to add the constraint because we don't recognize " *
        "$(func) as a valid JuMP function.",
    )
end

function build_constraint(
    ::Function,
    v::AbstractJuMPScalar,
    set::MOI.AbstractScalarSet,
)
    return ScalarConstraint(v, set)
end

function _clear_constant!(expr::Union{GenericAffExpr,GenericQuadExpr})
    offset = constant(expr)
    add_to_expression!(expr, -offset)
    return expr, offset
end

function _clear_constant!(α::Number)
    return zero(α), α
end

function build_constraint(
    ::Function,
    expr::Union{Number,GenericAffExpr,GenericQuadExpr},
    set::MOI.AbstractScalarSet,
)
    if MOI.Utilities.supports_shift_constant(typeof(set))
        expr, offset = _clear_constant!(expr)
        new_set = MOI.Utilities.shift_constant(set, -offset)
        return ScalarConstraint(expr, new_set)
    else
        return ScalarConstraint(expr, set)
    end
end

function build_constraint(
    error_fn::Function,
    ::_MA.Zero,
    set::MOI.AbstractScalarSet,
)
    return build_constraint(error_fn, false, set)
end

function build_constraint(
    ::Function,
    x::AbstractVector{<:Union{Number,AbstractJuMPScalar}},
    set::MOI.AbstractVectorSet,
)
    return VectorConstraint(x, set)
end

function build_constraint(
    error_fn::Function,
    x::AbstractArray,
    set::MOI.AbstractScalarSet,
)
    return error_fn(
        "Unexpected vector in scalar constraint. The left- and right-hand " *
        "sides of the constraint must have the same dimension.",
    )
end

function build_constraint(
    error_fn::Function,
    ::AbstractArray,
    ::AbstractVector,
    ::AbstractVector,
)
    return error_fn(
        "Unexpected vectors in scalar constraint. Did you mean to use the dot ",
        "comparison operators `l .<= f(x) .<= u` instead?",
    )
end

function build_constraint(
    error_fn::Function,
    x::AbstractMatrix,
    set::MOI.AbstractVectorSet,
)
    return error_fn(
        "unexpected matrix in vector constraint. Do you need to flatten the " *
        "matrix into a vector using `vec()`?",
    )
end

function build_constraint(
    error_fn::Function,
    ::Matrix,
    T::MOI.PositiveSemidefiniteConeTriangle,
)
    return error_fn("instead of `$(T)`, use `JuMP.PSDCone()`.")
end

# three-argument build_constraint is used for two-sided constraints.
function build_constraint(
    error_fn::Function,
    func::AbstractJuMPScalar,
    lb::Real,
    ub::Real,
)
    if isnan(lb) || isnan(ub)
        error_fn("Invalid bounds, cannot contain NaN: [$(lb), $(ub)].")
    end
    return build_constraint(
        error_fn,
        func,
        MOI.Interval(
            convert(value_type(variable_ref_type(func)), lb),
            convert(value_type(variable_ref_type(func)), ub),
        ),
    )
end

function build_constraint(
    error_fn::Function,
    ::AbstractJuMPScalar,
    ::Union{AbstractJuMPScalar,Real},
    ::Union{AbstractJuMPScalar,Real},
)
    return error_fn(
        "Interval constraint contains non-constant left- or " *
        "right-hand sides. Reformulate as two separate " *
        "constraints, or move all variables into the central term.",
    )
end

# This method intercepts `@constraint(model, lb <= var <= ub)` and promotes
# `var` to an `AffExpr` to form a `ScalarAffineFunction-in-Interval` instead of
# `VariableIndex-in-Interval`. To create a
# `MOI.VariableIndex`-in-`MOI.Interval`, use
# `@constraint(model, var in MOI.Interval(lb, ub))`. We do this for consistency
# with how one-sided (in)equality constraints are parsed.
function build_constraint(
    error_fn::Function,
    func::AbstractVariableRef,
    lb::Real,
    ub::Real,
)
    return build_constraint(
        error_fn,
        one(value_type(typeof(func))) * func,
        lb,
        ub,
    )
end

function build_constraint(
    ::Function,
    x::AbstractVector{T},
    set::MOI.SOS1,
) where {T<:AbstractJuMPScalar}
    return VectorConstraint(
        x,
        MOI.SOS1{value_type(variable_ref_type(T))}(set.weights),
    )
end

function build_constraint(
    ::Function,
    x::AbstractVector{T},
    set::MOI.SOS2,
) where {T<:AbstractJuMPScalar}
    return VectorConstraint(
        x,
        MOI.SOS2{value_type(variable_ref_type(T))}(set.weights),
    )
end
