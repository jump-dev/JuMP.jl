const ArrayIndices{N} = Base.Iterators.ProductIterator{NTuple{N, Base.OneTo{Int}}}
container(f::Function, indices) = container(f, indices, default_container(indices))
default_container(::ArrayIndices) = Array
function container(f::Function, indices::ArrayIndices, ::Type{Array})
    return map(I -> f(I...), indices)
end
default_container(::Base.Iterators.ProductIterator) = DenseAxisArray
function container(f::Function, indices::Base.Iterators.ProductIterator,
                   ::Type{DenseAxisArray})
    return DenseAxisArray(map(I -> f(I...), indices), indices.iterators...)
end
default_container(::NestedIterator) = SparseAxisArray
function container(f::Function, indices,
                   ::Type{SparseAxisArray})
    mappings = map(I -> I => f(I...), indices)
    data = Dict(mappings)
    if length(mappings) != length(data)
        unique_indices = Set()
        duplicate = nothing
        for index in indices
            if index in unique_indices
                duplicate = index
                break
            end
            push!(unique_indices, index)
        end
        # TODO compute idx
        error("Repeated index ", duplicate, ". Index sets must have unique elements.")
    end
    return SparseAxisArray(Dict(data))
end
