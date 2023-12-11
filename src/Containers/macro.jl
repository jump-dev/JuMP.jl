#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    container_name(x) --> Union{Symbol,Nothing}

Return the name of the container given by the expression `x`.

If the container is anonymous, return `nothing`.

## Examples

```jldoctest
julia> Containers.container_name(:x)
:x

julia> Containers.container_name(nothing)

julia> Containers.container_name(:(y[i in 1:2]))
:y

julia> Containers.container_name(:([i in 1:2]))
```
"""
function container_name(x)
    return error("Expression `$x::$(typeof(x))` cannot be used as a name.")
end


container_name(expr::Union{Symbol,Nothing}) = expr

function container_name(expr::Expr)
    if Meta.isexpr(expr, (:vcat, :vect))
        return nothing
    elseif Meta.isexpr(expr, (:ref, :typed_vcat))
        return container_name(expr.args[1])
    end
    return error("Expression $expr cannot be used as a name.")
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

"""
    parse_macro_arguments(error_fn::Function, args)

Returns a `Tuple{Vector{Any},Dict{Symbol,Any}}` containing the ordered
positional arguments and a dictionary mapping the keyword arguments.

This specially handles the distinction of `@foo(key = value)` and
`@foo(; key = value)` in macros.

Throws an error if mulitple keyword arguments are passed with the same name.
"""
function parse_macro_arguments(error_fn::Function, args)
    pos_args, kw_args = Any[], Dict{Symbol,Any}()
    for arg in _reorder_parameters(args)
        if Meta.isexpr(arg, :(=), 2)
            if haskey(kw_args, arg.args[1])
                error_fn(
                    "the keyword argument `$(arg.args[1])` was given " *
                    "multiple times.",
                )
            end
            kw_args[arg.args[1]] = arg.args[2]
        else
            push!(pos_args, arg)
        end
    end
    return pos_args, kw_args
end

"""
    _explicit_oneto(index_set)

If the `index_set` matches the form of `1:N`, then return
`Base.OneTo(index_set)`.
"""
function _explicit_oneto(index_set)
    s = Meta.isexpr(index_set, :escape) ? index_set.args[1] : index_set
    if Meta.isexpr(s, :call, 3) && s.args[1] == :(:) && s.args[2] == 1
        return :(Base.OneTo($index_set))
    else
        return index_set
    end
end

"""
    _expr_is_splat(expr)

Return `true` if `expr` is a `...` expression (or an `esc`'d one).
"""
function _expr_is_splat(expr)
    if Meta.isexpr(expr, :...)
        return true
    elseif Meta.isexpr(expr, :escape)
        return _expr_is_splat(expr.args[1])
    end
    return false
end

function _parse_index_sets(
    error_fn::Function,
    index_vars,
    index_sets,
    arg::Expr,
)
    index_var, index_set = gensym(), esc(arg)
    if Meta.isexpr(arg, :kw, 2) || Meta.isexpr(arg, :(=), 2)
        # Handle [i=S] and x[i=S]
        index_var, index_set = arg.args[1], esc(arg.args[2])
    elseif Meta.isexpr(arg, :call, 3) &&
           (arg.args[1] === :in || arg.args[1] === :∈)
        # Handle `i in S` and `i ∈ S`
        index_var, index_set = arg.args[2], esc(arg.args[3])
    end
    if index_var in index_vars
        error_fn(
            "The index $(index_var) appears more than once. The " *
            "index associated with each set must be unique.",
        )
    end
    push!(index_vars, index_var)
    push!(index_sets, index_set)
    return
end

function _parse_index_sets(::Function, index_vars, index_sets, arg)
    push!(index_vars, gensym())
    push!(index_sets, esc(arg))
    return
end

"""
    _parse_ref_sets(expr::Expr)

Helper function for macros to construct container objects.

Takes an `Expr` that specifies the container, e.g.,
`:(x[i=1:3,[:red,:blue],k=S; i+k <= 6])`, and returns:

 1. `index_vars`: Names for the index variables, e.g. `[:i, gensym(), :k]`.
    These may also be expressions, like `:((i, j))` from a call like
    `:(x[(i, j) in S])`.
 2. `index_sets`: Sets used for indexing, e.g. `[1:3, [:red,:blue], S]`
 3. `condition`: Expr containing any conditional imposed on indexing, or `:()`
    if none is present
"""
function _parse_ref_sets(error_fn::Function, expr::Expr)
    c = copy(expr)
    index_vars, index_sets, condition = Any[], Any[], :()
    # `:(t[i, j; k])` is a `:ref`, while `:(t[i; j])` is a `:typed_vcat`. In
    # both cases `:t` is the first argument.
    if Meta.isexpr(c, :typed_vcat) || Meta.isexpr(c, :ref)
        popfirst!(c.args)
    end
    if Meta.isexpr(c, :vcat) || Meta.isexpr(c, :typed_vcat)
        # An expression like `t[i; k]` or `[i; k]`. The filtering condition is
        # the second argument.
        if length(c.args) > 2
            error_fn(
                "Unsupported syntax $c: There can be at most one filtering " *
                "condition, which is separated from the indices by a single " *
                "`;`.",
            )
        elseif length(c.args) == 2
            condition = pop!(c.args)
        else
            # expr ends in a trailing `;`, but there is no condition
        end
    elseif Meta.isexpr(c, :ref) || Meta.isexpr(c, :vect)
        # An expression like `t[i, j; k]` or `[i, j; k]`. The filtering
        # condition is a `:parameters` expression in the first argument.
        if Meta.isexpr(c.args[1], :parameters)
            parameters = popfirst!(c.args)
            if length(parameters.args) != 1
                error_fn(
                    "Unsupported syntax $c: There can be at most one " *
                    "filtering condition, which is separated from the " *
                    "indices by a single `;`.",
                )
            end
            condition = parameters.args[1]
        end
    end
    for arg in c.args
        _parse_index_sets(error_fn, index_vars, index_sets, arg)
    end
    return index_vars, index_sets, condition
end

# Catch the case that has no index sets, just a name like `x`.
_parse_ref_sets(::Function, ::Union{Nothing,Symbol}) = (Any[], Any[], :())

depends_on(ex::Expr, s::Symbol) = any(Base.Fix2(depends_on, s), ex.args)

depends_on(ex::Symbol, s::Symbol) = ex == s

# For the case that `ex` might be an iterable literal like `4`.
depends_on(ex, s::Symbol) = false

depends_on(ex, s::QuoteNode) = depends_on(ex, s.value)

depends_on(ex, ::Any) = false

# For the case where the index set is compound, like `[(i, j) in S, k in i:K]`.
depends_on(ex1, ex2::Expr) = any(Base.Fix1(depends_on, ex1), ex2.args)

function _has_dependent_sets(index_vars::Vector{Any}, index_sets::Vector{Any})
    for i in 2:length(index_sets)
        for j in 1:(i-1)
            if depends_on(index_sets[i], index_vars[j])
                return true
            end
        end
    end
    return false
end

"""
    build_ref_sets(error_fn::Function, expr)

Helper function for macros to construct container objects.

!!! warning
    This function is for advanced users implementing JuMP extensions. See
    [`container_code`](@ref) for more details.

## Arguments

 * `error_fn`: a function that takes a `String` and throws an error, potentially
   annotating the input string with extra information such as from which macro
   it was thrown from. Use `error` if you do not want a modified error message.
 * `expr`: an `Expr` that specifies the container, e.g.,
   `:(x[i = 1:3, [:red, :blue], k = S; i + k <= 6])`

## Returns

 1. `index_vars`: a `Vector{Any}` of names for the index variables, e.g.,
    `[:i, gensym(), :k]`. These may also be expressions, like `:((i, j))` from a
    call like `:(x[(i, j) in S])`.
 2. `indices`: an iterator over the indices, for example,
    [`Containers.NestedIterator`](@ref)

## Example

See [`container_code`](@ref) for a worked example.
"""
function build_ref_sets(error_fn::Function, expr)
    index_vars, index_sets, condition = _parse_ref_sets(error_fn, expr)
    if any(_expr_is_splat, index_sets)
        error_fn(
            "cannot use splatting operator `...` in the definition of an " *
            "index set.",
        )
    end
    if !_has_dependent_sets(index_vars, index_sets) && condition == :()
        # Convert any 1:N to Base.OneTo(N)
        new_index_sets = _explicit_oneto.(index_sets)
        indices = :(Containers.vectorized_product($(new_index_sets...)))
        return index_vars, indices
    end
    esc_index_vars = esc.(index_vars)
    indices = Expr(:call, :(Containers.nested))
    for i in 1:length(index_vars)
        push!(
            indices.args,
            :(($(esc_index_vars[1:(i-1)]...),) -> $(index_sets[i])),
        )
    end
    if condition != :()
        f = :(($(esc_index_vars...),) -> $(esc(condition)))
        args = indices.args[2:end]
        indices = :(Containers.nested($(args...); condition = $f))
    end
    return index_vars, indices
end

"""
    container_code(
        index_vars::Vector{Any},
        indices::Expr,
        code,
        requested_container::Union{Symbol,Expr},
    )

Used in macros to construct a call to [`container`](@ref). This should be used
in conjunction with [`build_ref_sets`](@ref).

## Arguments

 * `index_vars::Vector{Any}`: a vector of names for the indices of the
   container. These may also be expressions, like `:((i, j))` from a
   call like `:(x[(i, j) in S])`.
 * `indices::Expr`: an expression that evaluates to an iterator of the indices.
 * `code`: an expression or literal constant for the value to be stored in the
   container as a function of the named `index_vars`.
 * `requested_container`: passed to the third argument of [`container`](@ref).
   For built-in JuMP types, choose one of `:Array`, `:DenseAxisArray`,
   `:SparseAxisArray`, or `:Auto`. For a user-defined container, this expression
   must evaluate to the correct type.

!!! warning
    In most cases, you should `esc(code)` before passing it to `container_code`.

## Example

```jldoctest
julia> macro foo(ref_sets, code)
           index_vars, indices = Containers.build_ref_sets(error, ref_sets)
           return Containers.container_code(
               index_vars,
               indices,
               esc(code),
               :Auto,
            )
       end
@foo (macro with 1 method)

julia> @foo(x[i=1:2, j=["A", "B"]], j^i)
2-dimensional DenseAxisArray{String,2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, ["A", "B"]
And data, a 2×2 Matrix{String}:
 "A"   "B"
 "AA"  "BB"
```
"""
function container_code(
    index_vars::Vector{Any},
    indices::Expr,
    code,
    requested_container::Union{Symbol,Expr},
)
    if isempty(index_vars)
        return code
    end
    esc_index_vars = esc.(index_vars)
    f = :(($(esc_index_vars...),) -> $code)
    # This switch handles the four "built-in" JuMP container types, with a
    # generic fallback for user-defined types.
    container_type = if requested_container == :Auto
        AutoContainerType
    elseif requested_container == :DenseAxisArray
        DenseAxisArray
    elseif requested_container == :SparseAxisArray
        SparseAxisArray
    elseif requested_container == :Array
        Array
    else
        # This is a symbol or expression from outside JuMP, so we need to escape
        # it.
        esc(requested_container)
    end
    return Expr(:call, container, f, indices, container_type, index_vars)
end

"""
    @container([i=..., j=..., ...], expr[, container = :Auto])

Create a container with indices `i`, `j`, ... and values given by `expr` that
may depend on the value of the indices.

    @container(ref[i=..., j=..., ...], expr[, container = :Auto])

Same as above but the container is assigned to the variable of name `ref`.

The type of container can be controlled by the `container` keyword.

!!! note
    When the index set is explicitly given as `1:n` for any expression `n`, it
    is transformed to `Base.OneTo(n)` before being given to [`container`](@ref).

## Example

```jldoctest
julia> Containers.@container([i = 1:3, j = 1:3], i + j)
3×3 Matrix{Int64}:
 2  3  4
 3  4  5
 4  5  6

julia> I = 1:3
1:3

julia> Containers.@container(x[i = I, j = I], i + j);

julia> x
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 1:3
    Dimension 2, 1:3
And data, a 3×3 Matrix{Int64}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.@container([i = 2:3, j = 1:3], i + j)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Base.OneTo(3)
And data, a 2×3 Matrix{Int64}:
 3  4  5
 4  5  6

julia> Containers.@container([i = 1:3, j = 1:3; i <= j], i + j)
SparseAxisArray{Int64, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  2
  [1, 2]  =  3
  [1, 3]  =  4
  [2, 2]  =  4
  [2, 3]  =  5
  [3, 3]  =  6
```
"""
macro container(input_args...)
    args, kw_args = parse_macro_arguments(error, input_args)
    container = get(kw_args, :container, :Auto)
    @assert length(args) == 2
    for key in keys(kw_args)
        @assert key == :container
    end
    index_vars, indices = build_ref_sets(error, args[1])
    code = container_code(index_vars, indices, esc(args[2]), container)
    name = container_name(args[1])
    if name === nothing
        return code
    end
    return :($(esc(name)) = $code)
end
