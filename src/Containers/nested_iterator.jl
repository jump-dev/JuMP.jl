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
struct NestedIterator{T}
    iterators::T # Tuple of functions
    condition::Function
end
function nested(iterators...; condition = (args...) -> true)
    return NestedIterator(iterators, condition)
end
Base.IteratorSize(::Type{<:NestedIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:NestedIterator}) = Base.EltypeUnknown()
function next_iterate(it::NestedIterator, i, elems, states, iterator, elem_state)
    if elem_state === nothing
        return nothing
    end
    elem, state = elem_state
    elems_states = first_iterate(
        it, i + 1, (elems..., elem),
        (states..., (iterator, state, elem)))
    if elems_states !== nothing
        return elems_states
    end
    return next_iterate(it, i, elems, states, iterator, iterate(iterator, state))
end
function first_iterate(it::NestedIterator, i, elems, states)
    if i > length(it.iterators)
        if it.condition(elems...)
            return elems, states
        else
            return nothing
        end
    end
    iterator = it.iterators[i](elems...)
    return next_iterate(it, i, elems, states, iterator, iterate(iterator))
end
function tail_iterate(it::NestedIterator, i, elems, states)
    if i > length(it.iterators)
        return nothing
    end
    next = tail_iterate(it, i + 1, (elems..., states[i][3]), states)
    if next !== nothing
        return next
    end
    iterator = states[i][1]
    next_iterate(it, i, elems, states[1:(i - 1)], iterator, iterate(iterator, states[i][2]))
end
Base.iterate(it::NestedIterator) = first_iterate(it, 1, tuple(), tuple())
Base.iterate(it::NestedIterator, states) = tail_iterate(it, 1, tuple(), states)
