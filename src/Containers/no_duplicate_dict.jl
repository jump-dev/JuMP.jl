#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    struct NoDuplicateDict{K, V} <: AbstractDict{K, V}
        dict::Dict{K, V}
    end

Same as `Dict{K, V}` but errors if constructed from an iterator with duplicate
keys.
"""
struct NoDuplicateDict{K, V} <: AbstractDict{K, V}
    dict::Dict{K, V}
    NoDuplicateDict{K, V}() where {K, V} = new{K, V}(Dict{K, V}())
end

# Implementation of the `AbstractDict` API.
Base.empty(::NoDuplicateDict, ::Type{K}, ::Type{V}) where {K, V} = NoDuplicateDict{K, V}()
Base.iterate(d::NoDuplicateDict, args...) = iterate(d.dict, args...)
Base.length(d::NoDuplicateDict) = length(d.dict)
Base.haskey(dict::NoDuplicateDict, key) = haskey(dict.dict, key)
Base.getindex(dict::NoDuplicateDict, key) = getindex(dict.dict, key)
function Base.setindex!(dict::NoDuplicateDict, value, key)
    if haskey(dict, key)
        error("Repeated index ", key, ". Index sets must have unique elements.")
    end
    setindex!(dict.dict, value, key)
end
function NoDuplicateDict{K, V}(it) where {K, V}
    dict = NoDuplicateDict{K, V}()
    for (k, v) in it
        dict[k] = v
    end
    return dict
end
function NoDuplicateDict(it)
    return Base.dict_with_eltype((K, V) -> NoDuplicateDict{K, V}, it, eltype(it))
end
