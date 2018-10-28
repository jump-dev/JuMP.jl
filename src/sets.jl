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

"""
    SecondOrderCone

Second order cone object that can be used to constrain the euclidean norm of a
vector `x` to be less than or euqal to a nonnegative scalar `t`. This is a
shortcut for the `MathOptInterface.SecondOrderCone`.

## Examples

The following constrains ``\\|(x-1, x-2)\\|_2 \\le t`` and ``t \\ge 0``:
```jldoctest; setup = :(using JuMP)
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, t)
t

julia> @constraint(model, [t, x-1, x-2] in SecondOrderCone())
[t, x - 1, x - 2] ∈ MathOptInterface.SecondOrderCone(3)
```
"""
struct SecondOrderCone <: AbstractVectorSet end
moi_set(::SecondOrderCone, dim::Int) = MOI.SecondOrderCone(dim)

"""
    RotatedSecondOrderCone

Rotated second order cone object that can be used to constrain the square of the
euclidean norm of a vector `x` to be less than or equal to ``2tu`` where `t` and
`u` are nonnegative scalars. This is a shortcut for the
`MathOptInterface.RotatedSecondOrderCone`.

## Examples

The following constrains ``\\|(x-1, x-2)\\|_2 \\le 2tx`` and ``t, x \\ge 0``:
```jldoctest; setup = :(using JuMP)
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, t)
t

julia> @constraint(model, [t, x, x-1, x-2] in RotatedSecondOrderCone())
[t, x, x - 1, x - 2] ∈ MathOptInterface.RotatedSecondOrderCone(4)
```
"""
struct RotatedSecondOrderCone <: AbstractVectorSet end
moi_set(::RotatedSecondOrderCone, dim::Int) = MOI.RotatedSecondOrderCone(dim)
