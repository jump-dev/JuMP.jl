#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module JuMPDataFramesExt

import DataFrames
import JuMP

function JuMP.Containers.container(
    f::F,
    indices::Union{
        JuMP.Containers.VectorizedProductIterator,
        JuMP.Containers.NestedIterator,
    },
    ::Type{DataFrames.DataFrame},
    names::AbstractVector,
) where {F<:Function}
    rows = vec(collect(indices))
    data = [name => getindex.(rows, i) for (i, name) in enumerate(names)]
    df = DataFrames.DataFrame(data)
    df.value = map(i -> f(i...), rows)
    return df
end

end #module
