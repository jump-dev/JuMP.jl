#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

# TODO add doc
depends_on(ex::Expr, s::Symbol) = any(a -> depends_on(a, s), ex.args)
depends_on(ex::Symbol, s::Symbol) = (ex == s)
depends_on(ex, s::Symbol) = false
function depends_on(ex1, ex2)
    @assert isa(ex2, Expr)
    @assert ex2.head == :tuple
    return any(s -> depends_on(ex1, s), ex2.args)
end

function is_dependent(idxvars, idxset, i)
    for (it, idx) in enumerate(idxvars)
        it == i && continue
        depends_on(idxset, idx) && return true
    end
    return false
end

function has_dependent_sets(idxvars, idxsets)
    # check if any index set depends on a previous index var
    for i in 2:length(idxsets)
        for v in idxvars[1:(i-1)]
            if depends_on(idxsets[i], v)
                return true
            end
        end
    end
    return false
end

valid_array_index_set(::Base.OneTo) = true
valid_array_index_set(r::UnitRange) = first(r) == 1
valid_array_index_set(s) = false

"""
    generate_container(T, indexvars, indexsets, requestedtype)

Return a tuple, the first element of which is *code* that generates a container
for objects of type `T` given the index variables, index sets, and
`requestedtype`. `requestedtype` may be one of `:Array`, `:DenseAxisArray`,
`:SparseAxisArray`, or `:Auto`. Return error-producing code if requested type is
incompatible. For the case of `:Auto`, the following rules are used to determine
the appropriate container:

1. If all index sets are either explicit `1:B` objects for any `B` or symbols
   which refer to objects of type `Base.OneTo`, then an `Array` is generated of
   the appropriate size. Types of symbols/expressions are not known at compile
   time, so we defer to type-safe functions to check the `Base.OneTo` condition.

2. If condition (1) does not hold, and the index sets are independent (the index
   variable for one set does not appear in the definition of another), then an
   `DenseAxisArray` is generated of the appropriate size.

3. Otherwise, generate an empty `SparseAxisArray{T,N,NTuple{N,Any}}`.

The second element of the return tuple is a `Bool`, `true` if the container type
automatically checks for duplicate terms in the index sets and `false` otherwise.

### Examples

    generate_container(VariableRef, [:i,:j], [:(1:N), :(1:T)], :Auto)
    # Returns code equivalent to:
    # :(Array{VariableRef}(length(1:N), length(1:T))

    generate_container(VariableRef, [:i,:j], [:(1:N), :(2:T)], :Auto)
    # Returns code equivalent to:
    # :(JuMP.Containers.DenseAxisArray(undef, 1:N, 2:T))

    generate_container(VariableRef, [:i,:j], [:(1:N), :(S)], :Auto)
    # Returns code that generates an Array if S is of type Base.OneTo,
    # otherwise an DenseAxisArray.

    generate_container(VariableRef, [:i,:j], [:(1:N), :(1:j)], :Auto)
    # Returns code equivalent to:
    # :(Containers.SparseAxisArray(Dict{NTuple{N,Any},VariableRef}()))
"""
function generate_container(T, indexvars, indexsets, requestedtype)
    has_dependent = has_dependent_sets(indexvars, indexsets)
    onetosets = falses(length(indexsets))
    for (i, indexset) in enumerate(indexsets)
        s = Meta.isexpr(indexset, :escape) ? indexset.args[1] : indexset
        if Meta.isexpr(s, :call) &&
           length(s.args) == 3 &&
           s.args[1] == :(:) &&
           s.args[2] == 1
            onetosets[i] = true
        end
    end

    if requestedtype == :SparseAxisArray ||
       (requestedtype == :Auto && has_dependent)
        N = length(indexvars)
        @assert N == length(indexsets)
        return :(JuMP.Containers.SparseAxisArray(Dict{NTuple{$N,Any},$T}())),
        false
    end

    sizes = Expr(:tuple, [:(length($rng)) for rng in indexsets]...)
    arrayexpr = :(Array{$T}(undef, $sizes...))

    if requestedtype == :Array
        has_dependent && return :(error(
            "Unable to create requested Array because index sets are dependent.",
        )),
        true
        if all(onetosets)
            return arrayexpr, true
        else
            # all sets must be one-based intervals
            condition = Expr(:&&)
            for (i, indexset) in enumerate(indexsets)
                if !onetosets[i]
                    push!(
                        condition.args,
                        Expr(
                            :call,
                            :(JuMP.Containers.valid_array_index_set),
                            indexset,
                        ),
                    )
                end
            end
            return :(
                $condition ||
                    error("Index set for array is not one-based interval.");
                $arrayexpr
            ),
            true
        end
    end

    axisexpr = :(JuMP.Containers.DenseAxisArray{$T}(undef))
    append!(axisexpr.args, indexsets)

    if requestedtype == :DenseAxisArray
        has_dependent && return :(error(
            "Unable to create requested DenseAxisArray because index sets are dependent.",
        )),
        true
        return axisexpr, true
    end

    @assert requestedtype == :Auto
    if all(onetosets) # definitely Array
        return arrayexpr, true
    end

    # Fallback, we have to decide once we know the types of the symbolic sets
    condition = Expr(:&&)
    for (i, indexset) in enumerate(indexsets)
        if !onetosets[i]
            push!(condition.args, Expr(:call, :isa, indexset, :(Base.OneTo)))
        end
    end

    # TODO: check type stability
    return :(
        if $condition
            $arrayexpr
        else
            $axisexpr
        end
    ), true
end
