#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

_get_name(c::Union{Symbol,AbstractString}) = c

function _get_name(c::Expr)
    if isexpr(c, :vcat) || isexpr(c, :vect)
        return Symbol("")  # Anonymous variable
    elseif isexpr(c, :ref) || isexpr(c, :typed_vcat)
        return _get_name(c.args[1])
    end
    return error("Expression $c cannot be used as a name.")
end

"""
    _extract_kw_args(args)

Process the arguments to a macro, separating out the keyword arguments.

Return a tuple of (flat_arguments, keyword arguments, and requested_container),
where `requested_container` is a symbol to be passed to `container_code`.
"""
function _extract_kw_args(args)
    flat_args, kw_args, requested_container = Any[], Any[], :Auto
    for arg in args
        if isexpr(arg, :(=))
            if arg.args[1] == :container
                requested_container = arg.args[2]
            else
                push!(kw_args, arg)
            end
        else
            push!(flat_args, arg)
        end
    end
    return flat_args, kw_args, requested_container
end

"""
    _explicit_oneto(index_set)

If the `index_set` matches the form of `1:N`, then return
`Base.OneTo(index_set)`.
"""
function _explicit_oneto(index_set)
    s = isexpr(index_set, :escape) ? index_set.args[1] : index_set
    if isexpr(s, :call, 3) && s.args[1] == :(:) && s.args[2] == 1
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
    if isexpr(expr, :...)
        return true
    elseif isexpr(expr, :escape)
        return _expr_is_splat(expr.args[1])
    end
    return false
end

function _parse_index_sets(_error::Function, index_vars, index_sets, arg::Expr)
    index_var, index_set = gensym(), esc(arg)
    if isexpr(arg, :kw, 2) || isexpr(arg, :(=), 2)
        # Handle [i=S] and x[i=S]
        index_var, index_set = arg.args[1], esc(arg.args[2])
    elseif isexpr(arg, :call, 3) && (arg.args[1] === :in || arg.args[1] === :∈)
        # Handle `i in S` and `i ∈ S`
        index_var, index_set = arg.args[2], esc(arg.args[3])
    end
    if index_var in index_vars
        _error(
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
function _parse_ref_sets(_error::Function, expr::Expr)
    c = copy(expr)
    index_vars, index_sets, condition = Any[], Any[], :()
    # `:(t[i, j; k])` is a `:ref`, while `:(t[i; j])` is a `:typed_vcat`. In
    # both cases `:t` is the first argument.
    if isexpr(c, :typed_vcat) || isexpr(c, :ref)
        popfirst!(c.args)
    end
    if isexpr(c, :vcat) || isexpr(c, :typed_vcat)
        # An expression like `t[i; k]` or `[i; k]`. The filtering condition is
        # the second argument.
        if length(c.args) > 2
            _error(
                "Unsupported syntax $c: There can be at most one filtering " *
                "condition, which is separated from the indices by a single " *
                "`;`.",
            )
        elseif length(c.args) == 2
            condition = pop!(c.args)
        else
            # expr ends in a trailing `;`, but there is no condition
        end
    elseif isexpr(c, :ref) || isexpr(c, :vect)
        # An expression like `t[i, j; k]` or `[i, j; k]`. The filtering
        # condition is a `:parameters` expression in the first argument.
        if isexpr(c.args[1], :parameters)
            parameters = popfirst!(c.args)
            if length(parameters.args) != 1
                _error(
                    "Unsupported syntax $c: There can be at most one " *
                    "filtering condition, which is separated from the " *
                    "indices by a single `;`.",
                )
            end
            condition = parameters.args[1]
        end
    end
    for arg in c.args
        _parse_index_sets(_error, index_vars, index_sets, arg)
    end
    return index_vars, index_sets, condition
end

# Catch the case that has no index sets, just a name like `x`.
_parse_ref_sets(::Function, ::Symbol) = (Any[], Any[], :())

_depends_on(ex::Expr, s::Symbol) = any(a -> _depends_on(a, s), ex.args)

_depends_on(ex::Symbol, s::Symbol) = ex == s

# For the case that `ex` might be an iterable literal like `4`.
_depends_on(ex, s::Symbol) = false

# For the case where the index set is compound, like `[(i, j) in S, k in i:K]`.
_depends_on(ex1, ex2::Expr) = any(s -> _depends_on(ex1, s), ex2.args)

function _has_dependent_sets(index_vars::Vector{Any}, index_sets::Vector{Any})
    for i in 2:length(index_sets)
        for j in 1:(i-1)
            if _depends_on(index_sets[i], index_vars[j])
                return true
            end
        end
    end
    return false
end

"""
    build_ref_sets(_error::Function, expr)

Helper function for macros to construct container objects.

!!! warning
    This function is for advanced users implementing JuMP extensions. See
    [`container_code`](@ref) for more details.

## Arguments

 * `_error`: a function that takes a `String` and throws an error, potentially
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
function build_ref_sets(_error::Function, expr)
    index_vars, index_sets, condition = _parse_ref_sets(_error, expr)
    if any(_expr_is_splat, index_sets)
        _error(
            "cannot use splatting operator `...` in the definition of an " *
            "index set.",
        )
    end
    if !_has_dependent_sets(index_vars, index_sets) && condition == :()
        # Convert any 1:N to Base.OneTo(N)
        new_index_sets = Containers._explicit_oneto.(index_sets)
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
        requested_container::Union{Symbol,Expr};
        pass_names::Union{Bool,Expr},
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
    requested_container::Union{Symbol,Expr};
    pass_names::Union{Bool,Expr} = false,
)
    if isempty(index_vars)
        return code
    end
    esc_index_vars = esc.(index_vars)
    f = :(($(esc_index_vars...),) -> $code)
    # This switch handles the four "built-in" JuMP container types, with a
    # generic fallback for user-defined types.
    container_type = if requested_container == :Auto
        Containers.AutoContainerType
    elseif requested_container == :DenseAxisArray
        Containers.DenseAxisArray
    elseif requested_container == :SparseAxisArray
        Containers.SparseAxisArray
    elseif requested_container == :Array
        Array
    else
        # This is a symbol or expression from outside JuMP, so we need to escape
        # it.
        pass_names = true
        esc(requested_container)
    end
    return quote
        names = $pass_names ? $index_vars : nothing
        Containers.container($f, $indices, $container_type, names)
    end
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
3×3 Array{Int64,2}:
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
And data, a 3×3 Array{Int64,2}:
 2  3  4
 3  4  5
 4  5  6

julia> Containers.@container([i = 2:3, j = 1:3], i + j)
2-dimensional DenseAxisArray{Int64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, Base.OneTo(3)
And data, a 2×3 Array{Int64,2}:
 3  4  5
 4  5  6

julia> Containers.@container([i = 1:3, j = 1:3; i <= j], i + j)
JuMP.Containers.SparseAxisArray{Int64,2,Tuple{Int64,Int64}} with 6 entries:
  [1, 2]  =  3
  [2, 3]  =  5
  [3, 3]  =  6
  [2, 2]  =  4
  [1, 1]  =  2
  [1, 3]  =  4
```
"""
macro container(args...)
    args, kwargs, requested_container = _extract_kw_args(args)
    @assert length(args) == 2
    enable_keyword_indexing = false
    for kw in kwargs
        @assert Meta.isexpr(kw, :(=), 2)
        @assert kw.args[1] == :enable_keyword_indexing
        enable_keyword_indexing = kw.args[2]
    end
    var, value = args
    index_vars, indices = build_ref_sets(error, var)
    code = container_code(
        index_vars,
        indices,
        esc(value),
        requested_container;
        pass_names = enable_keyword_indexing,
    )
    if isexpr(var, :vect) || isexpr(var, :vcat)
        return code
    else
        name = _get_name(var)
        return :($(esc(name)) = $code)
    end
end
