module _Derivatives

using LinearAlgebra
using SparseArrays

using Base.Meta
using ForwardDiff
import MathOptInterface
const MOI = MathOptInterface
using ..JuMP

# Needed for some of the univariate derivatives in univariate_expressions.jl.
using SpecialFunctions

const TAG = :jump_tag

# Override basic math functions to return NaN instead of throwing errors.
# This is what NLP solvers expect, and
# sometimes the results aren't needed anyway,
# because the code may compute derivatives wrt constants.
import NaNMath:
    sin,
    cos,
    tan,
    asin,
    acos,
    acosh,
    atanh,
    log,
    log2,
    log10,
    lgamma,
    log1p,
    pow,
    sqrt

include("univariate_expressions.jl")
include("types.jl")
include("conversion.jl")
include("coloring.jl")
export Coloring
include("linearity.jl")
include("sparsity.jl")
include("forward.jl")
include("reverse.jl")
include("subexpressions.jl")

end # module
