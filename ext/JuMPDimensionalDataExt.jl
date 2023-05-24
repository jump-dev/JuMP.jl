module JuMPDimensionalDataExt
using DimensionalData,JuMP

function JuMP.Containers.container(
    f::Function, 
    indices::JuMP.Containers.VectorizedProductIterator,
    c::Type{DimensionalData.DimArray},
    names,
)   
    for i in names
        if i isa Symbol && !Base.isidentifier(i) #a symbol was created automatically by JuMP.
            throw(ArgumentError("Cannot build DimArray. invalid indices."))
        end
    end
    dims = NamedTuple(i => j for (i, j) in zip(names, indices.prod.iterators))
    return DimensionalData.DimArray(map(i -> f(i...), indices), dims)
end

#error on nested iterators, TODO?
function JuMP.Containers.container(
    f::Function, 
    indices::JuMP.Containers.NestedIterator,
    c::Type{DimensionalData.DimArray},
    names,
)
    throw(MethodError(JuMP.Containers.container,(f,indices,c)))
end

end #module
