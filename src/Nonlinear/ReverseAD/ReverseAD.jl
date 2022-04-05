#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module ReverseAD

import ForwardDiff
import MathOptInterface
import ..Nonlinear
import SparseArrays

const MOI = MathOptInterface

include("Coloring/Coloring.jl")
include("graph_tools.jl")
include("types.jl")
include("utils.jl")

include("reverse_mode.jl")
include("forward_over_reverse.jl")
include("mathoptinterface_api.jl")

end  # module
