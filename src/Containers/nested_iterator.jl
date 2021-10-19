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

If `length(iterators) == 3`:
```julia
x = NestedIterator(iterators, condition)
for (i1, i2, i3) in x
    # produces (i1, i2, i3)
end
```
is the same as
```julia
for i1 in iterators[1]()
    for i2 in iterator[2](i1)
        for i3 in iterator[3](i1, i2)
            if condition(i1, i2, i3)
                # produces (i1, i2, i3)
            end
        end
    end
end
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

```julia
nested(1:2, ["A", "B"]; condition = (i, j) -> isodd(i) || j == "B")
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
