abstract type AbstractVectorSet end

"""
    moi_set(s::AbstractVectorSet, dim::Int)

Returns the MOI set of dimension `dim` corresponding to the JuMP set `s`.
"""
function moi_set end

# Used in `@constraint model f in s`
function build_constraint(_error::Function, f::AbstractVector,
                         s::AbstractVectorSet)
    return build_constraint(_error, f, moi_set(s, length(f)))
end

struct SecondOrderCone <: AbstractVectorSet end
moi_set(::SecondOrderCone, dim::Int) = MOI.SecondOrderCone(dim)
struct RotatedSecondOrderCone <: AbstractVectorSet end
moi_set(::RotatedSecondOrderCone, dim::Int) = MOI.RotatedSecondOrderCone(dim)
