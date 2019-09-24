const ArrayIndices{N} = Iterators.ProductIterator{NTuple{N, Base.OneTo{Int}}}
container(f::Function, indices) = container(f, indices, default_container(indices))
default_container(::ArrayIndices) = Array
function container(f::Function, indices::ArrayIndices, ::Type{Array})
    return map(I -> f(I...), indices)
end
function _oneto(indices)
    if indices isa UnitRange{Int} && indices == 1:length(indices)
        return Base.OneTo(length(indices))
    end
    error("Index set for array is not one-based interval.")
end
function container(f::Function, indices::Iterators.ProductIterator,
                   ::Type{Array})
    container(f, Iterators.ProductIterator(_oneto.(indices.iterators)), Array)
end
default_container(::Iterators.ProductIterator) = DenseAxisArray
function container(f::Function, indices::Iterators.ProductIterator,
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
