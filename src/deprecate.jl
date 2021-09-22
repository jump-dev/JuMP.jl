#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function Base.getproperty(
    x::Union{NonlinearExpression,NonlinearParameter,NLPEvaluator},
    key::Symbol,
)
    if key == :m
        @warn("`.m` is deprecated in favor of `.model`.")
        return getfield(x, :model)
    end
    return getfield(x, key)
end
