#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import Base: promote_rule, promote_type, cat_t, hcat, vcat, hvcat, cat

immutable DummyJuMPArray end

Base.promote_rule{T<:OneIndexedArray,S<:OneIndexedArray}(::Type{T},::Type{S}) = DummyJuMPArray
Base.promote_rule{T<:OneIndexedArray,S<:Union(AbstractArray,Number,JuMPTypes)}(::Type{T},::Type{S}) = DummyJuMPArray

# The following methods are needed on v0.3 since T(x) does not fall back to convert(T,x)
AffExpr(v::Variable) = AffExpr(Variable[v], Float64[1], 0.0)
AffExpr(aff::AffExpr) = aff
QuadExpr(v::Variable) = QuadExpr(Variable[], Variable[], Float64[], AffExpr(v))
QuadExpr(aff::AffExpr) = QuadExpr(Variable[], Variable[], Float64[], aff)
QuadExpr(q::QuadExpr) = q

_tofull(x) = x
_tofull(x::OneIndexedArray) = x.innerArray

Base.hcat(X::OneIndexedArray...) = hcat([_tofull(x) for x in X]...)
Base.vcat(X::OneIndexedArray...) = vcat([_tofull(x) for x in X]...)
Base.hvcat(rows::(Int...), X::OneIndexedArray...) = hvcat(rows, [_tofull(x) for x in X]...)

typealias BaseTypes Union(Number,AbstractArray)
typealias CatTypes Union(Number,JuMPTypes,AbstractArray)
typealias CatTypesFull Union(CatTypes,OneIndexedArray)

Base.vcat(X::CatTypes...) = Base.cat(1, X...)
Base.hcat(X::CatTypes...) = Base.cat(2, X...)

Base.vcat(X::Variable...) = Variable[ X[i] for i=1:length(X) ]

###############################################################################
# copied from Base at https://github.com/JuliaLang/julia/blob/0c24dca65c031820b91721139f0291068086955c/base/abstractarray.jl#L613-L668
# This version is exactly identical to what's in Base so we (hopefully) don't
# mess with anything.
# The Julia language is licensed under the MIT License. More details about the
# licensing can be found at http://julialang.org/license
function Base.cat(catdim::Integer, X::BaseTypes...)
    nargs = length(X)
    dimsX = map((a->isa(a,AbstractArray) ? size(a) : (1,)), X)
    ndimsX = map((a->isa(a,AbstractArray) ? ndims(a) : 1), X)
    d_max = maximum(ndimsX)

    if catdim > d_max + 1
        for i=1:nargs
            if dimsX[1] != dimsX[i]
                error("all inputs must have same dimensions when concatenating along a higher dimension");
            end
        end
    elseif nargs >= 2
        for d=1:d_max
            if d == catdim; continue; end
            len = d <= ndimsX[1] ? dimsX[1][d] : 1
            for i = 2:nargs
                if len != (d <= ndimsX[i] ? dimsX[i][d] : 1)
                    error("mismatch in dimension ", d)
                end
            end
        end
    end

    cat_ranges = [ catdim <= ndimsX[i] ? dimsX[i][catdim] : 1 for i=1:nargs ]

    function compute_dims(d)
        if d == catdim
            if catdim <= d_max
                return sum(cat_ranges)
            else
                return nargs
            end
        else
            if d <= ndimsX[1]
                return dimsX[1][d]
            else
                return 1
            end
        end
    end

    ndimsC = max(catdim, d_max)
    dimsC = ntuple(ndimsC, compute_dims)::(Int...)
    typeC = promote_type(map(x->isa(x,AbstractArray) ? eltype(x) : typeof(x), X)...)
    C = similar(isa(X[1],AbstractArray) ? full(X[1]) : [X[1]], typeC, dimsC)

    range = 1
    for k=1:nargs
        nextrange = range+cat_ranges[k]
        cat_one = [ i != catdim ? (1:dimsC[i]) : (range:nextrange-1) for i=1:ndimsC ]
        C[cat_one...] = X[k]
        range = nextrange
    end
    return C
end
# End copy
###############################################################################

###############################################################################
# copied from Base at https://github.com/JuliaLang/julia/blob/0c24dca65c031820b91721139f0291068086955c/base/abstractarray.jl#L613-L668
function Base.cat(catdim::Integer, X::CatTypes...)
    nargs = length(X)
    dimsX = map((a->isa(a,AbstractArray) ? size(a) : (1,)), X)
    ndimsX = map((a->isa(a,AbstractArray) ? ndims(a) : 1), X)
    d_max = maximum(ndimsX)

    if catdim > d_max + 1
        for i=1:nargs
            if dimsX[1] != dimsX[i]
                error("all inputs must have same dimensions when concatenating along a higher dimension");
            end
        end
    elseif nargs >= 2
        for d=1:d_max
            if d == catdim; continue; end
            len = d <= ndimsX[1] ? dimsX[1][d] : 1
            for i = 2:nargs
                if len != (d <= ndimsX[i] ? dimsX[i][d] : 1)
                    error("mismatch in dimension ", d)
                end
            end
        end
    end

    cat_ranges = [ catdim <= ndimsX[i] ? dimsX[i][catdim] : 1 for i=1:nargs ]

    function compute_dims(d)
        if d == catdim
            if catdim <= d_max
                return sum(cat_ranges)
            else
                return nargs
            end
        else
            if d <= ndimsX[1]
                return dimsX[1][d]
            else
                return 1
            end
        end
    end

    ndimsC = max(catdim, d_max)
    dimsC = ntuple(ndimsC, compute_dims)::(Int...)
    typeC = promote_type(map(x->isa(x,AbstractArray) ? eltype(x) : typeof(x), X)...)
    # add typeC typed vcat to avoid stack overflow with [x] for x::Variable
    C = similar(isa(X[1],AbstractArray) ? full(X[1]) : typeC[X[1]], typeC, dimsC)

    range = 1
    for k=1:nargs
        nextrange = range+cat_ranges[k]
        cat_one = [ i != catdim ? (1:dimsC[i]) : (range:nextrange-1) for i=1:ndimsC ]
        C[cat_one...] = X[k]
        range = nextrange
    end
    return C
end
# End copy
###############################################################################

Base.vcat(X::CatTypesFull...) = cat(1, [_tofull(x) for x in X]...)
Base.hcat(X::CatTypesFull...) = cat(2, [_tofull(x) for x in X]...)
Base.cat(catdim::Integer, X::CatTypesFull...) = cat(catdim, [_tofull(x) for x in X]...)
