#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# With `Iterators.ProductIterator`, everything works except when the
# `IteratorSize` of an iterator is not `Base.HasShape{1}` or `Base.HasLength`
# which is notably the case for scalars for it is `Base.HasShape{0}` and
# multidimensional arrays for which it is `Base.HasShape{N}` where `N` is the
# dimension. For instance:
# ```julia
# julia> collect(Iterators.product(2, 3))
# 0-dimensional Array{Tuple{Int64,Int64},0}:
# (2, 3)
# ```
# while we could like it to be a `3`-dimensional array of size `(1, 1)`.
# When the user does `@container([2, 3], 1)`, a `DenseAxisArray` of size
# `(1, 1)`. Another example:
# ```julia
# julia> collect(Iterators.product([1 2; 3 4]))
# 2×2 Array{Tuple{Int64},2}:
# (1,)  (2,)
# (3,)  (4,)
# ```
# while we need the size to be `(4,)`, not `(2, 2)` when the user does
# `@container([i = [1, 2; 3 4]], i^2)`.
# Long story short, we want to tried everything as a interator without shape
# while `Iterators.ProductIterator` does care about preserving the shape
# when doing the cartesian product.
"""
    struct VectorizedProductIterator{T}
        prod::Iterators.ProductIterator{T}
    end

Cartesian product of the iterators `prod.iterators`. It is the same iterator as
`Base.Iterators.ProductIterator` except that it is independent of the
`IteratorSize` of the elements of `prod.iterators`.
For instance:
* `size(Iterators.product(1, 2))` is `tuple()` while
  `size(VectorizedProductIterator(1, 2))` is `(1, 1)`.
* `size(Iterators.product(ones(2, 3)))` is `(2, 3)` while
  `size(VectorizedProductIterator(ones(2, 3)))` is `(1, 1)`.
"""
struct VectorizedProductIterator{T}
    prod::Iterators.ProductIterator{T}
end
function vectorized_product(iterators...)
    return VectorizedProductIterator(Iterators.product(iterators...))
end
function Base.IteratorSize(::Type{<:VectorizedProductIterator{<:Tuple{Vararg{Any, N}}}}) where N
    return Base.HasShape{N}()
end
Base.IteratorEltype(::Type{<:VectorizedProductIterator}) = Base.EltypeUnknown()
Base.size(it::VectorizedProductIterator) = _prod_size(it.prod.iterators)
_prod_size(::Tuple{}) = ()
_prod_size(t::Tuple) = (length(t[1]), _prod_size(Base.tail(t))...)
Base.axes(it::VectorizedProductIterator) = _prod_indices(it.prod.iterators)
_prod_indices(::Tuple{}) = ()
_prod_indices(t::Tuple) = (Base.OneTo(length(t[1])), _prod_indices(Base.tail(t))...)
Base.ndims(it::VectorizedProductIterator) = length(axes(it))
Base.length(it::VectorizedProductIterator) = prod(size(it))
Base.iterate(it::VectorizedProductIterator, args...) = iterate(it.prod, args...)
