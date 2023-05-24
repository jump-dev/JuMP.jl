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
    names::AbstractVector,
) where {F<:Function}
    dims = NamedTuple(i => j for (i, j) in zip(names, indices.prod.iterators))
    return DimensionalData.DimArray(map(i -> f(i...), indices), dims)
end

function JuMP.Containers.container(
    ::Function,
    ::JuMP.Containers.NestedIterator,
    ::Type{DimensionalData.DimArray},
    ::AbstractVector,
)
    return error(
        "Unable to create a `DimArray` because the container does not form " *
        "a dense rectangular array",
    )
end

end #module
