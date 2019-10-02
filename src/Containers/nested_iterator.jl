#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    struct NestedIterator{T}
        iterators::T # Tuple of functions
        condition::Function
    end

Iterators over the tuples that are produced by a nested for loop.
For instance, if `length(iterators) == 3` , this corresponds to the tuples
`(i1, i2, i3)` produced by:
```
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
struct NestedIterator{T, C}
    iterators::T # Tuple of functions
    condition::C
end
function nested(iterators...; condition = (args...) -> true)
    return NestedIterator(iterators, condition)
end
Base.IteratorSize(::Type{<:NestedIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:NestedIterator}) = Base.EltypeUnknown()
function next_iterate(iterators, condition, elems, states, iterator, elem_state)
    if elem_state === nothing
        return nothing
    end
    elem, state = elem_state
    elems_states = first_iterate(
        Base.tail(iterators), condition, (elems..., elem),
        (states..., (iterator, state, elem)))
    if elems_states !== nothing
        return elems_states
    end
    return next_iterate(iterators, condition, elems, states, iterator, iterate(iterator, state))
end
function first_iterate(::Tuple{}, condition, elems, states)
    if condition(elems...)
        return elems, states
    else
        return nothing
    end
end
function first_iterate(iterators, condition, elems, states)
    iterator = iterators[1](elems...)
    return next_iterate(iterators, condition, elems, states, iterator, iterate(iterator))
end
tail_iterate(::Tuple{}, condition, elems, states, prev_states) = nothing
function tail_iterate(iterators, condition, elems, states, prev_states)
    next = tail_iterate(Base.tail(iterators), condition, (elems..., states[1][3]), Base.tail(states), (prev_states..., states[1]))
    if next !== nothing
        return next
    end
    iterator = states[1][1]
    next_iterate(iterators, condition, elems, prev_states, iterator, iterate(iterator, states[1][2]))
end
Base.iterate(it::NestedIterator) = first_iterate(it.iterators, it.condition, tuple(), tuple())
Base.iterate(it::NestedIterator, states) = tail_iterate(it.iterators, it.condition, tuple(), states, tuple())
