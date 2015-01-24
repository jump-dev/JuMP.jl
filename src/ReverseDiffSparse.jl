module ReverseDiffSparse

using Compat
import Calculus
using DualNumbers
using Base.Meta
# Override basic math functions to return NaN instead of throwing errors.
# This is what NLP solvers expect, and
# sometimes the results aren't needed anyway,
# because the code may compute derivatives wrt constants.
import NaNMath: sin, cos, tan, asin, acos, acosh, atanh, log, log2, log10, lgamma, log1p, pow

if isdir(Pkg.dir("ArrayViews"))
    eval(Expr(:import,:ArrayViews))
    const subarr = ArrayViews.view
else
    const subarr = Base.sub
end

issum(s::Symbol) = (s == :sum) || (s == :∑) || (s == :Σ)
isprod(s::Symbol) = (s == :prod) || (s == :∏)

# package code goes here
include("types.jl")
include("revmode.jl")
include("hessian.jl")
include("coloring.jl")
include("exprlist.jl")
include("export.jl")

end # module
