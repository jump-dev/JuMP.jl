module ReverseDiffSparse

using Compat
import Calculus
import DualNumbers: Dual, Dual4, epsilon, epsilon1, epsilon2, epsilon3, epsilon4
using Base.Meta
# Override basic math functions to return NaN instead of throwing errors.
# This is what NLP solvers expect, and
# sometimes the results aren't needed anyway,
# because the code may compute derivatives wrt constants.
import NaNMath: sin, cos, tan, asin, acos, acosh, atanh, log, log2, log10, lgamma, log1p, pow

const subarr = ArrayViews.view

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
