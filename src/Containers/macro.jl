#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _reorder_parameters(args)
    if isempty(args) || !Meta.isexpr(args[1], :parameters)
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
    parse_macro_arguments(
        error_fn::Function,
        args;
        valid_kwargs::Union{Nothing,Vector{Symbol}} = nothing,
        num_positional_args::Union{Nothing,Int,UnitRange{Int}} = nothing,
    )

Returns a `Tuple{Vector{Any},Dict{Symbol,Any}}` containing the ordered
positional arguments and a dictionary mapping the keyword arguments.

This specially handles the distinction of `@foo(key = value)` and
`@foo(; key = value)` in macros.

An error is thrown if multiple keyword arguments are passed with the same key.

If `valid_kwargs` is a `Vector{Symbol}`, an error is thrown if a keyword is not
in `valid_kwargs`.

If `num_positional_args` is not nothing, an error is thrown if the number of
positional arguments is not in `num_positional_args`.
"""
function parse_macro_arguments(
    error_fn::Function,
    args;
    valid_kwargs::Union{Nothing,Vector{Symbol}} = nothing,
    num_positional_args::Union{Nothing,Int,UnitRange{Int}} = nothing,
)
    pos_args, kwargs = Any[], Dict{Symbol,Any}()
    for arg in _reorder_parameters(args)
        if Meta.isexpr(arg, :(=), 2)
            if haskey(kwargs, arg.args[1])
                error_fn(
                    "the keyword argument `$(arg.args[1])` was given " *
                    "multiple times.",
                )
            elseif valid_kwargs !== nothing && !(arg.args[1] in valid_kwargs)
                error_fn("unsupported keyword argument `$(arg.args[1])`.")
            end
            kwargs[arg.args[1]] = arg.args[2]
        else
            push!(pos_args, arg)
        end
    end
    if num_positional_args isa Int
        n = length(pos_args)
        if n != num_positional_args
            error_fn(
                "expected $num_positional_args positional arguments, got $n.",
            )
        end
    elseif num_positional_args isa UnitRange{Int}
        if !(length(pos_args) in num_positional_args)
            a, b = num_positional_args.start, num_positional_args.stop
            error_fn(
                "expected $a to $b positional arguments, got $(length(pos_args)).",
            )
        end
    end
    return pos_args, kwargs
end

"""
    _explicit_oneto(error_fn, index_set)

If the `index_set` matches the form of `1:N`, then return
`Base.OneTo(index_set)`.
"""
function _explicit_oneto(error_fn, index_set)
    s = Meta.isexpr(index_set, :escape) ? index_set.args[1] : index_set
    index_set = quote
        try
            $index_set
        catch
            $error_fn(
                "unexpected error parsing reference set: ",
                $(Meta.quot(_drop_esc(index_set))),
            )
        end
    end
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

Takes an `Expr` that specifies the container, for example,
`:(x[i=1:3,[:red,:blue],k=S; i+k <= 6])`, and returns:

 1. `index_vars`: Names for the index variables, for example, `[:i, gensym(), :k]`.
    These may also be expressions, like `:((i, j))` from a call like
    `:(x[(i, j) in S])`.
 2. `index_sets`: Sets used for indexing, for example, `[1:3, [:red,:blue], S]`
 3. `condition`: Expr containing any conditional imposed on indexing, or `:()`
    if none is present
"""
function _parse_ref_sets(error_fn::Function, expr::Expr)
    c = copy(expr)
    index_vars, index_sets, condition = Any[], Any[], :()
    # `:(t[i, j; k])` is a `:ref`, while `:(t[i; j])` is a `:typed_vcat`. In
    # both cases `:t` is the first argument.
    if Meta.isexpr(c, :typed_vcat) || Meta.isexpr(c, :ref)
        name = popfirst!(c.args)
        if !(name isa Symbol)
            error_fn(
                "Unsupported syntax: the expression `$name` cannot be used " *
                "as a name.",
            )
        end
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

function _container_name(error_fn::Function, x)
    return error_fn("Expression `$x::$(typeof(x))` cannot be used as a name.")
end

_container_name(::Function, expr::Union{Symbol,Nothing}) = expr

function _container_name(error_fn::Function, expr::Expr)
    if Meta.isexpr(expr, (:vcat, :vect))
        return nothing
    elseif Meta.isexpr(expr, (:ref, :typed_vcat))
        return _container_name(error_fn, expr.args[1])
    end
    return error_fn("Expression $expr cannot be used as a name.")
end

"""
    parse_ref_sets(
        error_fn::Function,
        expr;
        invalid_index_variables::Vector{Symbol} = Symbol[],
    )

Helper function for macros to construct container objects.

!!! warning
    This function is for advanced users implementing JuMP extensions. See
    [`container_code`](@ref) for more details.

## Arguments

 * `error_fn`: a function that takes a `String` and throws an error, potentially
   annotating the input string with extra information such as from which macro
   it was thrown from. Use `error` if you do not want a modified error message.
 * `expr`: an `Expr` that specifies the container, for example,
   `:(x[i = 1:3, [:red, :blue], k = S; i + k <= 6])`

## Returns

 1. `name`: the name of the container, if given, otherwise `nothing`
 2. `index_vars`: a `Vector{Any}` of names for the index variables, for example,
    `[:i, gensym(), :k]`. These may also be expressions, like `:((i, j))` from a
    call like `:(x[(i, j) in S])`.
 3. `indices`: an iterator over the indices, for example,
    [`Containers.NestedIterator`](@ref)

## Example

See [`container_code`](@ref) for a worked example.
"""
function parse_ref_sets(
    error_fn::Function,
    expr::Union{Nothing,Symbol,Expr};
    invalid_index_variables::Vector = Symbol[],
)
    name = _container_name(error_fn, expr)
    index_vars, indices = build_ref_sets(error_fn, expr)
    for name in invalid_index_variables
        if name in index_vars
            error_fn(
                "the index name `$name` conflicts with another variable in " *
                "this scope. Use a different name for the index.",
            )
        end
    end
    return name, index_vars, indices
end

# This method is needed because Julia v1.10 prints LineNumberNode in the string
# representation of an expression.
function _strip_LineNumberNode(x::Expr)
    if Meta.isexpr(x, :block)
        return Expr(:block, filter(!Base.Fix2(isa, LineNumberNode), x.args)...)
    end
    return x
end

_strip_LineNumberNode(x) = x

"""
    build_error_fn(macro_name, args, source)

Return a function that can be used in place of `Base.error`, but which
additionally prints the macro from which it was called.
"""
function build_error_fn(macro_name, args, source)
    str_args = join(_strip_LineNumberNode.(args), ", ")
    msg = "At $(source.file):$(source.line): `@$macro_name($str_args)`: "
    error_fn(str...) = error(msg, str...)
    return error_fn
end

_drop_esc(x) = Meta.isexpr(x, :escape) ? x.args[1] : x

"""
    build_ref_sets(error_fn::Function, expr)

This function is deprecated. Use [`parse_ref_sets`](@ref) instead.
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
        new_index_sets = _explicit_oneto.(error_fn, index_sets)
        indices = :(Containers.vectorized_product($(new_index_sets...)))
        return index_vars, indices
    end
    esc_index_vars = esc.(index_vars)
    indices = Expr(:call, :(Containers.nested))
    for i in 1:length(index_vars)
        push!(
            indices.args,
            quote
                ($(esc_index_vars[1:(i-1)]...),) -> try
                    $(index_sets[i])
                catch
                    $error_fn(
                        "unexpected error parsing reference set: ",
                        $(Meta.quot(_drop_esc(index_sets[i]))),
                    )
                end
            end,
        )
    end
    if condition != :()
        f = quote
            ($(esc_index_vars...),) -> try
                $(esc(condition))
            catch
                $error_fn(
                    "unexpected error parsing condition: ",
                    $(Meta.quot(condition)),
                )
            end
        end
        args = indices.args[2:end]
        indices = :(Containers.nested($(args...); condition = $f))
    end
    return index_vars, indices
end

"""
    add_additional_args(
        call::Expr,
        args::Vector,
        kwargs::Dict{Symbol,Any};
        kwarg_exclude::Vector{Symbol} = Symbol[],
    )

Add the positional arguments `args` to the function call expression `call`,
escaping each argument expression.

This function is able to incorporate additional positional arguments to `call`s
that already have keyword arguments.
"""
function add_additional_args(
    call::Expr,
    args::Vector,
    kwargs::Dict{Symbol,Any};
    kwarg_exclude::Vector{Symbol} = Symbol[],
)
    call_args = call.args
    if Meta.isexpr(call, :.)
        # call is broadcasted
        call_args = call.args[2].args
    end
    # Cache all keyword arguments
    kw_args = filter(Base.Fix2(Meta.isexpr, :kw), call_args)
    # Remove keyword arguments from the end
    filter!(!Base.Fix2(Meta.isexpr, :kw), call_args)
    # Add the new positional arguments
    append!(call_args, esc.(args))
    # Re-add the cached keyword arguments back to the end
    append!(call_args, kw_args)
    for (key, value) in kwargs
        if !(key in kwarg_exclude)
            push!(call_args, esc(Expr(:kw, key, value)))
        end
    end
    return
end

"""
    build_name_expr(
        name::Union{Symbol,Nothing},
        index_vars::Vector,
        kwargs::Dict{Symbol,Any},
    )

Returns an expression for the name of a container element, where `name` and
`index_vars` are the values returned by [`parse_ref_sets`](@ref) and `kwargs`
is the dictionary returned by [`parse_macro_arguments`](@ref).

This assumes that the key in `kwargs` used to over-ride the name choice is
`:base_name`.

## Example

```jldoctest
julia> Containers.build_name_expr(:x, [:i, :j], Dict{Symbol,Any}())
:(string("x", "[", string(\$(Expr(:escape, :i))), ",", string(\$(Expr(:escape, :j))), "]"))

julia> Containers.build_name_expr(nothing, [:i, :j], Dict{Symbol,Any}())
""

julia> Containers.build_name_expr(:y, [:i, :j], Dict{Symbol,Any}(:base_name => "y"))
:(string("y", "[", string(\$(Expr(:escape, :i))), ",", string(\$(Expr(:escape, :j))), "]"))
```
"""
function build_name_expr(
    name::Union{Symbol,Nothing},
    index_vars::Vector,
    kwargs::Dict{Symbol,Any},
)
    base_name = get(kwargs, :base_name, string(something(name, "")))
    if !(base_name isa String)
        base_name = esc(base_name)
    end
    if isempty(index_vars) || base_name == ""
        return base_name
    end
    expr = Expr(:call, :string, base_name, "[")
    for index in index_vars
        # Converting the arguments to strings before concatenating is faster:
        # https://github.com/JuliaLang/julia/issues/29550.
        push!(expr.args, :(string($(esc(index)))))
        push!(expr.args, ",")
    end
    expr.args[end] = "]"
    return expr
end

"""
    container_code(
        index_vars::Vector{Any},
        indices::Expr,
        code,
        requested_container::Union{Symbol,Expr,Dict{Symbol,Any}},
    )

Used in macros to construct a call to [`container`](@ref). This should be used
in conjunction with [`parse_ref_sets`](@ref).

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
   must evaluate to the correct type. You may also pass the `kwargs` dictionary
   from [`parse_macro_arguments`](@ref).

!!! warning
    In most cases, you should `esc(code)` before passing it to `container_code`.

## Example

```jldoctest
julia> macro foo(ref_sets, code)
           name, index_vars, indices =
               Containers.parse_ref_sets(error, ref_sets)
           @assert name !== nothing  # Anonymous container not supported
           container =
               Containers.container_code(index_vars, indices, esc(code), :Auto)
           return quote
               \$(esc(name)) = \$container
           end
       end
@foo (macro with 1 method)

julia> @foo(x[i=1:2, j=["A", "B"]], j^i);

julia> x
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

function container_code(
    index_vars::Vector{Any},
    indices::Expr,
    code,
    kwargs::Dict{Symbol,Any},
)
    container = get(kwargs, :container, :Auto)
    return container_code(index_vars, indices, code, container)
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
    args, kwargs = parse_macro_arguments(
        error,
        input_args;
        num_positional_args = 2,
        valid_kwargs = [:container],
    )
    name, index_vars, indices = parse_ref_sets(error, args[1])
    code = container_code(index_vars, indices, esc(args[2]), kwargs)
    if name === nothing
        return code
    end
    return :($(esc(name)) = $code)
end

"""
    _extract_kw_args(args)

!!! warning
    This function is deprecated. Use [`parse_macro_arguments`](@ref) instead.
"""
function _extract_kw_args(args)
    flat_args, kw_args, requested_container = Any[], Any[], :Auto
    for arg in args
        if Meta.isexpr(arg, :(=))
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
