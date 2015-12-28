module ReverseDiffSparse2

using Base.Meta
import DualNumbers: Dual, epsilon
import Calculus
import Lazy

include("types.jl")
include("conversion.jl")
include("linearity.jl")
include("sparsity.jl")
include("forward.jl")
include("reverse.jl")
include("coloring.jl")
export Coloring
include("mpb_wrapper.jl")

end # module
