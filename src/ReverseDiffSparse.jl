module ReverseDiffSparse

import Calculus
using DualNumbers
using Base.Meta
if Pkg.installed("ArrayViews") != nothing
    eval(Expr(:import,:ArrayViews))
    const subarr = ArrayViews.view
else
    const subarr = Base.sub
end

# package code goes here
include("types.jl")
include("revmode.jl")
include("hessian.jl")
include("coloring.jl")
include("exprlist.jl")

end # module
