module ReverseDiffSparse

import Calculus
using DualNumbers
using Base.Meta

# package code goes here
include("types.jl")
include("revmode.jl")
include("hessian.jl")

end # module
