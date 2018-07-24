abstract type AbstractVectorSet end

"""
    moi_set(s::AbstractVectorSet, dim::Int)

Returns the MOI set of dimension `dim` corresponding to the JuMP set `s`.
"""
function moi_set end

# Used in `@constraint model f in s`
function buildconstraint(_error::Function, f::AbstractVector,
                         s::AbstractVectorSet)
    return buildconstraint(_error, f, moi_set(s, length(f)))
end

struct SOCone <: AbstractVectorSet end
moi_set(::SOCone, dim::Int) = MOI.SecondOrderCone(dim)
struct RSOCone <: AbstractVectorSet end
moi_set(::RSOCone, dim::Int) = MOI.RotatedSecondOrderCone(dim)
