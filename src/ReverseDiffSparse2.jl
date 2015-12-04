module ReverseDiffSparse2

using Base.Meta
import Calculus
import Lazy

include("types.jl")
include("conversion.jl")
include("linearity.jl")
include("sparsity.jl")
include("forward.jl")
include("reverse.jl")
include("coloring.jl")

end # module
