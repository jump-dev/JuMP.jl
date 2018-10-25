#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

struct JuMPArray{T,N,NT} <: JuMPContainer{T,N}
    innerArray::Array{T,N}
    indexsets::NT
    lookup::NTuple{N,Any}
    meta::Dict{Symbol,Any}
end

@generated function JuMPArray(innerArray::Array{T,N}, indexsets::NTuple{N,Any}) where {T,N}
    dicttuple = Expr(:tuple)
    for i in 1:N
        inner = quote
            idxset = indexsets[$i]
            ret = Dict{eltype(idxset), Int}()
        end
        tupelem = indexsets.parameters[i]
        if !(tupelem == UnitRange{Int} || tupelem == StepRange{Int})
            inner = quote
                $inner
                cnt = 1
                for x in idxset
                    ret[x] = cnt
                    cnt += 1
                end
                ret
            end
        end
        push!(dicttuple.args, inner)
    end
    :(JuMPArray(innerArray, indexsets, $dicttuple, Dict{Symbol,Any}()))
end

# Note that this is needed because automatic broadcasting of addition operators is deprecated in 1.0.
# The explicitly broadcasted version ought to work in 0.6, but causes a mysterious generated function
# related error because of the presence of this function in `@generated Base.setindex!`.
if VERSION < v"0.7-"
    _to_cartesian_range_sub(I, rng) = I - (first(rng) - 1)
else
    _to_cartesian_range_sub(I, rng) = I .- (first(rng) .- 1)
end

function _to_cartesian(d,NT,idx...)
    indexing = Any[]
    for (i,S) in enumerate(NT.parameters)
        idxtype = idx[1][i]
        if S == UnitRange{Int}
            if idxtype == Colon
                # special stuff
                push!(indexing, Colon())
            elseif idxtype <: AbstractRange
                push!(indexing, quote
                    rng = d.indexsets[$i]
                    I = idx[$i]
                    _to_cartesian_range_sub(I, rng)
                end)
            else
                push!(indexing, quote
                    rng = d.indexsets[$i]
                    I = idx[$i]
                    first(rng) <= I <= last(rng) || error("Failed attempt to index JuMPArray along dimension $($i): $I ∉ $(d.indexsets[$i])")
                    _to_cartesian_range_sub(I, rng)
                end)
            end
        elseif S == StepRange{Int}
            if idx[1][i] == Colon
                push!(indexing, Colon())
            else
                push!(indexing, quote
                    rng = $(d.indexsets[i])
                    I = idx[$i]
                    first(rng) <= I <= last(rng) || error("Failed attempt to index JuMPArray along dimension $($i): $I ∉ $(d.indexsets[$i])")
                    dv, rv = divrem(I - first(rng), step(rng))
                    rv == 0 || error("Failed attempt to index JuMPArray along dimension $($i): $I ∉ $(d.indexsets[$i])")
                    dv + 1
                end)
            end
        else
            push!(indexing, quote
                if !haskey(d.lookup[$i],idx[$i])
                    error("Failed attempt to index JuMPArray along dimension $($i): $(idx[$i]) ∉ $(d.indexsets[$i])")
                end
                d.lookup[$i][idx[$i]]::Int
            end)
        end
    end
    indexing
end

Base.getindex(d::JuMPArray, ::Colon) = d.innerArray[:]

@generated function Base.getindex(d::JuMPArray{T,N,NT}, idx...) where {T,N,NT}
    if N != length(idx)
        error("Indexed into a JuMPArray with $(length(idx)) indices (expected $N indices)")
    end
    Expr(:call, :getindex, :(d.innerArray), _to_cartesian(d,NT,idx)...)
end

@generated function Base.setindex!(d::JuMPArray{T,N,NT}, v, idx...) where {T,N,NT}
    if N != length(idx)
        error("Indexed into a JuMPArray with $(length(idx)) indices (expected $N indices)")
    end
    Expr(:call, :setindex!, :(d.innerArray), :v, _to_cartesian(d,NT,idx)...)
end
