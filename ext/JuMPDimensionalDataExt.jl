#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module JuMPDimensionalDataExt

import DimensionalData
import JuMP

function JuMP.Containers.container(
    f::F, 
    indices::JuMP.Containers.VectorizedProductIterator,
    c::Type{DimensionalData.DimArray},
    names,
) where {F<:Function}
    for i in names
        if i isa Symbol && !Base.isidentifier(i)
            # A symbol was created automatically by JuMP.
            throw(ArgumentError("Cannot build DimArray. invalid indices."))
        end
    end
    dims = NamedTuple(i => j for (i, j) in zip(names, indices.prod.iterators))
    return DimensionalData.DimArray(map(i -> f(i...), indices), dims)
end

# TODO: imporve error for nested iterators
function JuMP.Containers.container(
    f::Function, 
    indices::JuMP.Containers.NestedIterator,
    c::Type{DimensionalData.DimArray},
    names,
)
    return throw(MethodError(JuMP.Containers.container, (f, indices, c)))
end

end #module
