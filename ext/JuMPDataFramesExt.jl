#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module JuMPDataFramesExt

import DataFrames
import JuMP

function JuMP.Containers.container(
    f::Function,
    indices,
    ::Type{DataFrames.DataFrame},
    names::AbstractVector,
)
    rows = vec(collect(indices))
    df = DataFrames.DataFrame(NamedTuple{tuple(names...)}(arg) for arg in rows)
    df.value = [f(arg...) for arg in rows]
    return df
end

end #module
