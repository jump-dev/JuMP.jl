#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta

_get_name(c::Symbol) = c
_get_name(c::Nothing) = ()
_get_name(c::AbstractString) = c
function _get_name(c::Expr)
    if c.head == :string
        return c
    else
        return c.args[1]
    end
end

"""
    _extract_kw_args(args)

Process the arguments to a macro, separating out the keyword arguments.
Return a tuple of (flat_arguments, keyword arguments, and requested_container),
where `requested_container` is a symbol to be passed to `parse_container`.
"""
function _extract_kw_args(args)
    kw_args = filter(x -> isexpr(x, :(=)) && x.args[1] != :container , collect(args))
    flat_args = filter(x->!isexpr(x, :(=)), collect(args))
    requested_container = :Auto
    for kw in args
        if isexpr(kw, :(=)) && kw.args[1] == :container
            requested_container = kw.args[2]
        end
    end
    return flat_args, kw_args, requested_container
end

function _try_parse_idx_set(arg::Expr)
    # [i=1] and x[i=1] parse as Expr(:vect, Expr(:(=), :i, 1)) and
    # Expr(:ref, :x, Expr(:kw, :i, 1)) respectively.
    if arg.head === :kw || arg.head === :(=)
        @assert length(arg.args) == 2
        return true, arg.args[1], arg.args[2]
    elseif isexpr(arg, :call) && arg.args[1] === :in
        return true, arg.args[2], arg.args[3]
    else
        return false, nothing, nothing
    end
end
function _explicit_oneto(index_set)
    s = Meta.isexpr(index_set,:escape) ? index_set.args[1] : index_set
    if Meta.isexpr(s,:call) && length(s.args) == 3 && s.args[1] == :(:) && s.args[2] == 1
        return :(Base.OneTo($index_set))
    else
        return index_set
    end
end

function _expr_is_splat(ex::Expr)
    if ex.head == :(...)
        return true
    elseif ex.head == :escape
        return _expr_is_splat(ex.args[1])
    end
    return false
end
_expr_is_splat(::Any) = false

"""
    _parse_ref_sets(expr::Expr)

Helper function for macros to construct container objects. Takes an `Expr` that
specifies the container, e.g. `:(x[i=1:3,[:red,:blue],k=S; i+k <= 6])`, and
returns:

    1. `idxvars`: Names for the index variables, e.g. `[:i, gensym(), :k]`
    2. `idxsets`: Sets used for indexing, e.g. `[1:3, [:red,:blue], S]`
    3. `condition`: Expr containing any conditional imposed on indexing, or `:()` if none is present
"""
function _parse_ref_sets(_error::Function, expr::Expr)
    c = copy(expr)
    idxvars = Any[]
    idxsets = Any[]
    # `:(t[i,j;k])` is a :ref, while `:(t[i;j])` is a :typed_vcat.
    # In both cases `:t` is the first arg.
    if isexpr(c, :typed_vcat) || isexpr(c, :ref)
        popfirst!(c.args)
    end
    condition = :()
    if isexpr(c, :vcat) || isexpr(c, :typed_vcat)
        # Parameters appear as plain args at the end.
        if length(c.args) > 2
            _error("Unsupported syntax $c.")
        elseif length(c.args) == 2
            condition = pop!(c.args)
        end # else no condition.
    elseif isexpr(c, :ref) || isexpr(c, :vect)
        # Parameters appear at the front.
        if isexpr(c.args[1], :parameters)
            if length(c.args[1].args) != 1
                _error("Invalid syntax: $c. Multiple semicolons are not " *
                       "supported.")
            end
            condition = popfirst!(c.args).args[1]
        end
    end
    if isexpr(c, :vcat) || isexpr(c, :typed_vcat) || isexpr(c, :ref)
        if isexpr(c.args[1], :parameters)
            @assert length(c.args[1].args) == 1
            condition = popfirst!(c.args).args[1]
        end # else no condition.
    end

    for s in c.args
        parse_done = false
        if isa(s, Expr)
            parse_done, idxvar, _idxset = _try_parse_idx_set(s::Expr)
            if parse_done
                idxset = esc(_idxset)
            end
        end
        if !parse_done # No index variable specified
            idxvar = gensym()
            idxset = esc(s)
        end
        push!(idxvars, idxvar)
        push!(idxsets, idxset)
    end
    return idxvars, idxsets, condition
end
_parse_ref_sets(_error::Function, expr) = (Any[], Any[], :())

"""
    _build_ref_sets(_error::Function, expr)

Helper function for macros to construct container objects. Takes an `Expr` that
specifies the container, e.g. `:(x[i=1:3,[:red,:blue],k=S; i+k <= 6])`, and
returns:

    1. `idxvars`: Names for the index variables, e.g. `[:i, gensym(), :k]`
    2. `indices`: Iterators over the indices indexing, e.g.
       `Constainers.NestedIterators((1:3, [:red,:blue], S), (i, ##..., k) -> i + k <= 6)`.
"""
function _build_ref_sets(_error::Function, expr)
    idxvars, idxsets, condition = _parse_ref_sets(_error, expr)
    if any(_expr_is_splat.(idxsets))
        _error("cannot use splatting operator `...` in the definition of an index set.")
    end
    has_dependent = has_dependent_sets(idxvars, idxsets)
    if has_dependent || condition != :()
        esc_idxvars = esc.(idxvars)
        idxfuns = [:(($(esc_idxvars[1:(i - 1)]...),) -> $(idxsets[i])) for i in 1:length(idxvars)]
        if condition == :()
            indices = :(Containers.nested($(idxfuns...)))
        else
            condition_fun = :(($(esc_idxvars...),) -> $(esc(condition)))
            indices = :(Containers.nested($(idxfuns...); condition = $condition_fun))
        end
    else
        indices = :(Containers.vectorized_product($(_explicit_oneto.(idxsets)...)))
    end
    return idxvars, indices
end

function container_code(idxvars, indices, code, requested_container)
    if isempty(idxvars)
        return code
    end
    if !(requested_container in [:Auto, :Array, :DenseAxisArray, :SparseAxisArray])
        # We do this two-step interpolation, first into the string, and then
        # into the expression because interpolating into a string inside an
        # expression has scoping issues.
        error_message = "Invalid container type $requested_container. Must be " *
                        "Auto, Array, DenseAxisArray, or SparseAxisArray."
        return :(error($error_message))
    end
    if requested_container == :DenseAxisArray
        requested_container = :(JuMP.Containers.DenseAxisArray)
    elseif requested_container == :SparseAxisArray
        requested_container = :(JuMP.Containers.SparseAxisArray)
    end
    esc_idxvars = esc.(idxvars)
    func = :(($(esc_idxvars...),) -> $code)
    if requested_container == :Auto
        return :(Containers.container($func, $indices))
    else
        return :(Containers.container($func, $indices, $requested_container))
    end
end
function parse_container(_error, var, value, requested_container)
    idxvars, indices = _build_ref_sets(_error, var)
    return container_code(idxvars, indices, value, requested_container)
end

"""
    @container([i=..., j=..., ...], expr)

Create a container with indices `i`, `j`, ... and values given by `expr` that
may depend on the value of the indices.

    @container(ref[i=..., j=..., ...], expr)

Same as above but the container is assigned to the variable of name `ref`.

The type of container can be controlled by the `container` keyword. See
[Containers in macros](@ref). Note that when the index set is explicitly
given as `1:n` for any expression `n`, it is transformed to `Base.OneTo(n)`
before being given to [`container`](@ref).

## Examples

```jldoctest
julia> Containers.@container([i = 1:3, j = 1:3], i + j)
3×3 Array{Int64,2}:
 2  3  4
 3  4  5
 4  5  6

julia> I = 1:3
1:3

julia> Containers.@container([i = I, j = I], i + j)
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
    args, kw_args, requested_container = _extract_kw_args(args)
    @assert length(args) == 2
    @assert isempty(kw_args)
    var, value = args
    code = parse_container(error, var, esc(value), requested_container)
    anonvar = isexpr(var, :vect) || isexpr(var, :vcat)
    if anonvar
        return code
    else
        name = Containers._get_name(var)
        return :($(esc(name)) = $code)
    end
end
