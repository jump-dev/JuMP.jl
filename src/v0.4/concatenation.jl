#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import Base: promote_rule, promote_type, cat_t, hcat, vcat, hvcat

immutable DummyJuMPArray end

Base.promote_rule{T<:OneIndexedArray,S<:OneIndexedArray}(::Type{T},::Type{S}) = DummyJuMPArray
Base.promote_rule{T<:OneIndexedArray,S<:Union(AbstractArray,Number,JuMPTypes)}(::Type{T},::Type{S}) = DummyJuMPArray

_tofull(x) = x
_tofull(x::OneIndexedArray) = x.innerArray

function Base.cat_t(catdims, ::Type{DummyJuMPArray}, X...)
    Y = map(_tofull, X)
    T = promote_type(map(x->isa(x,AbstractArray) ? eltype(x) : typeof(x), Y)...)
    @assert T != DummyJuMPArray # else we hit an infinite recursion below...
    cat_t(catdims, T, Y...)
end

Base.hcat(X::OneIndexedArray...) = hcat([_tofull(x) for x in X]...)
Base.vcat(X::OneIndexedArray...) = vcat([_tofull(x) for x in X]...)
Base.hvcat(rows::Tuple{Vararg{Int}}, X::OneIndexedArray...) = hvcat(rows, [_tofull(x) for x in X]...)
