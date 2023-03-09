#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    struct NestedIterator{T}
        iterators::T # Tuple of functions
        condition::Function
    end

Iterators over the tuples that are produced by a nested for loop.

Construct a `NestedIterator` using [`nested`](@ref).

## Example

```jldoctest nested_iterator_docstring
julia> iterators = (() -> 1:2, (i,) -> ["A", "B"]);

julia> condition = (i, j) -> isodd(i) || j == "B";

julia> x = Containers.NestedIterator(iterators, condition);

julia> for (i, j) in x
           println((i, j))
       end
(1, "A")
(1, "B")
(2, "B")
```
is the same as
```jldoctest nested_iterator_docstring
julia> for i in iterators[1]()
           for j in iterators[2](i)
               if condition(i, j)
                   println((i, j))
               end
           end
       end
(1, "A")
(1, "B")
(2, "B")
```
"""
struct NestedIterator{T,C}
    iterators::T # Tuple of functions
    condition::C
end

"""
    nested(iterators...; condition = (args...) -> true)

Create a [`NestedIterator`](@ref).

## Example

```jldoctest
julia> iterator = Containers.nested(
           () -> 1:2,
           (i,) -> ["A", "B"];
           condition = (i, j) -> isodd(i) || j == "B",
       );

julia> collect(iterator)
3-element Vector{Tuple{Int64, String}}:
 (1, "A")
 (1, "B")
 (2, "B")
```
"""
function nested(iterators...; condition = (args...) -> true)
    return NestedIterator(iterators, condition)
end

Base.IteratorSize(::Type{<:NestedIterator}) = Base.SizeUnknown()

Base.IteratorEltype(::Type{<:NestedIterator}) = Base.EltypeUnknown()

function _next_iterate(
    iterators,
    condition,
    elems,
    states,
    iterator,
    elem_state,
)
    while true
        if elem_state === nothing
            return
        end
        elem, state = elem_state
        elems_states = _first_iterate(
            Base.tail(iterators),
            condition,
            (elems..., elem),
            (states..., (iterator, state, elem)),
        )
        if elems_states !== nothing
            return elems_states
        end
        # This could be written as a recursive function where we call
        # `_next_iterate` here with this new value of `_next_iterate` instead of
        # the `while` loop`.
        #
        # However, if there are too many consecutive elements for which
        # `condition`is `false` for the last iterator, this will result in a
        # `StackOverflow`. See https://github.com/jump-dev/JuMP.jl/issues/2335
        elem_state = iterate(iterator, state)
    end
end

function _first_iterate(::Tuple{}, condition, elems, states)
    if condition(elems...)
        return elems, states
    end
    return
end

function _first_iterate(iterators, condition, elems, states)
    iterator = iterators[1](elems...)
    return _next_iterate(
        iterators,
        condition,
        elems,
        states,
        iterator,
        iterate(iterator),
    )
end

_tail_iterate(::Tuple{}, condition, elems, states, prev_states) = nothing

function _tail_iterate(iterators, condition, elems, states, prev_states)
    next = _tail_iterate(
        Base.tail(iterators),
        condition,
        (elems..., states[1][3]),
        Base.tail(states),
        (prev_states..., states[1]),
    )
    if next !== nothing
        return next
    end
    iterator = states[1][1]
    return _next_iterate(
        iterators,
        condition,
        elems,
        prev_states,
        iterator,
        iterate(iterator, states[1][2]),
    )
end

function Base.iterate(it::NestedIterator)
    return _first_iterate(it.iterators, it.condition, tuple(), tuple())
end

function Base.iterate(it::NestedIterator, states)
    return _tail_iterate(it.iterators, it.condition, tuple(), states, tuple())
end

function _eltype_or_any(::NestedIterator{<:Tuple{Vararg{Any,N}}}) where {N}
    return NTuple{N,Any}
end
