#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    @constraint(model, expr, args...; kwargs...)
    @constraint(model, [index_sets...], expr, args...; kwargs...)
    @constraint(model, name, expr, args...; kwargs...)
    @constraint(model, name[index_sets...], expr, args...; kwargs...)

Add a constraint described by the expression `expr`.

The `name` argument is optional. If index sets are passed, a container is built
and the constraint may depend on the indices of the index sets.

The expression `expr` may be one of following forms:

 * `func in set`, constraining the function `func` to belong to the set `set`,
   which is either a [`MOI.AbstractSet`](@ref) or one of the JuMP shortcuts like
   [`SecondOrderCone`](@ref) or [`PSDCone`](@ref)

 * `a <op> b`, where `<op>` is one of `==`, `≥`, `>=`, `≤`, `<=`

 * `l <= f <= u` or `u >= f >= l`, constraining the expression `f` to lie
   between `l` and `u`

 * `f(x) ⟂ x`, which defines a complementarity constraint

 * `z --> {expr}`, which defines an indicator constraint that activates
   when `z` is `1`

 * `!z --> {expr}`, which defines an indicator constraint that activates
   when `z` is `0`

 * `z <--> {expr}`, which defines a reified constraint

 * `expr := rhs`, which defines a Boolean equality constraint

Broadcasted comparison operators like `.==` are also supported for the case when
the left- and right-hand sides of the comparison operator are arrays.

JuMP extensions may additionally provide support for constraint expressions
which are not listed here.

## Keyword arguments

 * `base_name`: sets the name prefix used to generate constraint names. It
   corresponds to the constraint name for scalar constraints, otherwise, the
   constraint names are set to `base_name[...]` for each index `...`.

 * `container = :Auto`: force the container type by passing `container = Array`,
  `container = DenseAxisArray`, `container = SparseAxisArray`, or any another
  container type which is supported by a JuMP extension.

 * `set_string_name::Bool = true`: control whether to set the [`MOI.ConstraintName`](@ref)
   attribute. Passing `set_string_name = false` can improve performance.

Other keyword arguments may be supported by JuMP extensions.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @variable(model, z, Bin);

julia> @constraint(model, x in SecondOrderCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.SecondOrderCone(3)

julia> @constraint(model, [i in 1:3], x[i] == i)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 x[1] = 1
 x[2] = 2
 x[3] = 3

julia> @constraint(model, x .== [1, 2, 3])
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 x[1] = 1
 x[2] = 2
 x[3] = 3

julia> @constraint(model, con_name, 1 <= x[1] + x[2] <= 3)
con_name : x[1] + x[2] ∈ [1, 3]

julia> @constraint(model, con_perp[i in 1:3], x[i] - 1 ⟂ x[i])
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VectorAffineFunction{Float64}, MathOptInterface.Complements}, VectorShape}}:
 con_perp[1] : [x[1] - 1, x[1]] ∈ MathOptInterface.Complements(2)
 con_perp[2] : [x[2] - 1, x[2]] ∈ MathOptInterface.Complements(2)
 con_perp[3] : [x[3] - 1, x[3]] ∈ MathOptInterface.Complements(2)

julia> @constraint(model, z --> {x[1] >= 0})
z --> {x[1] ≥ 0}

julia> @constraint(model, !z --> {2 * x[2] <= 3})
!z --> {2 x[2] ≤ 3}
```
"""
macro constraint(input_args...)
    error_fn = Containers.build_error_fn(:constraint, input_args, __source__)
    args, kwargs = Containers.parse_macro_arguments(error_fn, input_args)
    if length(args) < 2 && !isempty(kwargs)
        error_fn(
            "No constraint expression detected. If you are trying to " *
            "construct an equality constraint, use `==` instead of `=`.",
        )
    elseif length(args) < 2
        error_fn("expected 2 to 4 positional arguments, got $(length(args)).")
    elseif Meta.isexpr(args[2], :block)
        error_fn("Invalid syntax. Did you mean to use `@constraints`?")
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
    c, x = nothing, y
    if y isa Symbol || Meta.isexpr(y, (:vect, :vcat, :ref, :typed_vcat))
        if length(extra) == 0
            error_fn("No constraint expression was given.")
        end
        c, x = y, popfirst!(extra)
    end
    if length(extra) > 1
        error_fn("Cannot specify more than 1 additional positional argument.")
    end
    name, index_vars, indices = Containers.parse_ref_sets(
        error_fn,
        c;
        invalid_index_variables = [args[1]],
    )
    is_vectorized, parse_code, build_call = parse_constraint(error_fn, x)
    Containers.add_additional_args(
        build_call,
        extra,
        kwargs;
        kwarg_exclude = [:base_name, :container, :set_string_name],
    )
    # ; set_string_name
    name_expr = Containers.build_name_expr(name, index_vars, kwargs)
    if name_expr != ""
        set_string_name = if haskey(kwargs, :set_string_name)
            esc(kwargs[:set_string_name])
        else
            :(set_string_names_on_creation($model))
        end
        name_expr = :($set_string_name ? $name_expr : "")
    end
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
        Containers.container_code(index_vars, indices, code, kwargs),
        __source__;
        register_name = name,
        wrap_let = true,
    )
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

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @build_constraint(x .>= 0)
2-element Vector{ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}}:
 ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}(x[1], MathOptInterface.GreaterThan{Float64}(-0.0))
 ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}(x[2], MathOptInterface.GreaterThan{Float64}(-0.0))
```
"""
macro build_constraint(arg)
    error_fn = Containers.build_error_fn(:build_constraint, (arg,), __source__)
    _, parse_code, build_call = parse_constraint(error_fn, arg)
    return quote
        $parse_code
        $build_call
    end
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

function parse_constraint(error_fn::Function, arg)
    return error_fn(
        "Incomplete constraint specification $arg. Are you missing a " *
        "comparison (<=, >=, or ==)?",
    )
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
(for example, a constant vector or matrix). Due to promotion and arithmetic, this
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
    ::Union{Val{:vect},Val{:vcat}},
    args...,
)
    return error_fn(
        """
        Unsupported constraint expression: we don't know how to parse a
        `[ ]` block as a constraint. Have you written:
        ```julia
        @constraint(model, name, [...], ...)
        ```
        instead of:
        ```julia
        @constraint(model, name[...], ...)
        ```""",
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
[x[1], x[2]] ∈ Nonnegatives()

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x >= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ Nonnegatives()
```
"""
struct Nonnegatives end

"""
    GreaterThanZero()

A struct used to intercept when `>=` or `≥` is used in a macro via
[`operator_to_set`](@ref).

This struct is not the same as [`Nonnegatives`](@ref) so that we can disambiguate
`x >= y` and `x - y in Nonnegatives()`.

This struct is not intended for general usage, but it may be useful to some
JuMP extensions.

## Example

```jldoctest
julia> operator_to_set(error, Val(:>=))
GreaterThanZero()
```
"""
struct GreaterThanZero end

operator_to_set(::Function, ::Union{Val{:(>=)},Val{:(≥)}}) = GreaterThanZero()

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
[x[1], x[2]] ∈ Nonpositives()

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x <= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ Nonpositives()
```
"""
struct Nonpositives end

"""
    GreaterThanZero()

A struct used to intercept when `<=` or `≤` is used in a macro via
[`operator_to_set`](@ref).

This struct is not the same as [`Nonpositives`](@ref) so that we can disambiguate
`x <= y` and `x - y in Nonpositives()`.

This struct is not intended for general usage, but it may be useful to some
JuMP extensions.

## Example

```jldoctest
julia> operator_to_set(error, Val(:<=))
LessThanZero()
```
"""
struct LessThanZero end

operator_to_set(::Function, ::Union{Val{:(<=)},Val{:(≤)}}) = LessThanZero()

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
[x[1], x[2]] ∈ Zeros()

julia> A = [1 2; 3 4];

julia> b = [5, 6];

julia> @constraint(model, A * x == b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ Zeros()
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
        :(build_constraint.($error_fn, _desparsify($f), $(esc(set))))
    else
        :(build_constraint($error_fn, $f, $(esc(set))))
    end
    return parse_code, build_call
end

function _functionize(v::V) where {V<:AbstractVariableRef}
    return convert(GenericAffExpr{value_type(V),V}, v)
end

_functionize(v::AbstractArray{<:AbstractVariableRef}) = _functionize.(v)

function _functionize(
    v::LinearAlgebra.Symmetric{V},
) where {V<:AbstractVariableRef}
    return LinearAlgebra.Symmetric(_functionize(v.data))
end

function _functionize(
    v::LinearAlgebra.Hermitian{V},
) where {V<:AbstractVariableRef}
    return LinearAlgebra.Hermitian(_functionize(v.data))
end

_functionize(x) = x

_functionize(::_MA.Zero) = false

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
    ::GreaterThanZero,
    args...;
    kwargs...,
)
    return build_constraint(error_fn, f, Nonnegatives(), args...; kwargs...)
end

function build_constraint(
    error_fn::Function,
    ::Union{Matrix,LinearAlgebra.Symmetric,LinearAlgebra.Hermitian},
    ::GreaterThanZero,
)
    return error_fn(
        """

        The syntax `x >= y` is ambiguous for matrices because we cannot tell if
        you intend a positive semidefinite constraint or an elementwise
        inequality.

        To create a positive semidefinite constraint, pass `PSDCone()` or
        `HermitianPSDCone()`:

        ```julia
        @constraint(model, x >= y, PSDCone())
        ```

        To create an element-wise inequality, pass `Nonnegatives()`, or use
        broadcasting:

        ```julia
        @constraint(model, x >= y, Nonnegatives())
        # or
        @constraint(model, x .>= y)
        ```""",
    )
end

function build_constraint(
    error_fn::Function,
    f,
    ::LessThanZero,
    args...;
    kwargs...,
)
    return build_constraint(error_fn, f, Nonpositives(), args...; kwargs...)
end

function build_constraint(
    error_fn::Function,
    ::Union{Matrix,LinearAlgebra.Symmetric,LinearAlgebra.Hermitian},
    ::LessThanZero,
)
    return error_fn(
        """

        The syntax `x <= y` is ambiguous for matrices because we cannot tell if
        you intend a positive semidefinite constraint or an elementwise
        inequality.

        To create a positive semidefinite constraint, reverse the sense of the
        inequality and pass `PSDCone()` or `HermitianPSDCone()`:

        ```julia
        @constraint(model, y >= x, PSDCone())
        ```

        To create an element-wise inequality, reverse the sense of the
        inequality and pass `Nonnegatives()`, or use broadcasting:

        ```julia
        @constraint(model, y >= x, Nonnegatives())
        # or
        @constraint(model, x .<= y)
        ```""",
    )
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

"""
    build_constraint(error_fn::Function, func, set, args...; kwargs...)

This method should only be implemented by developers creating JuMP extensions.
It should never be called by users of JuMP.
"""
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
