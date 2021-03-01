abstract type AbstractVectorSet end

# Used in `@constraint(model, [1:n] in s)`
function build_variable(
    _error::Function,
    variables::Vector{<:ScalarVariable},
    set::AbstractVectorSet,
)
    return VariablesConstrainedOnCreation(
        variables,
        moi_set(set, length(variables)),
    )
end

# Used in `@constraint(model, func in set)`
function build_constraint(
    _error::Function,
    func::AbstractVector,
    set::AbstractVectorSet,
)
    return build_constraint(_error, func, moi_set(set, length(func)))
end

"""
    SecondOrderCone

Second order cone object that can be used to constrain the euclidean norm of a
vector `x` to be less than or equal to a nonnegative scalar `t`. This is a
shortcut for the `MOI.SecondOrderCone`.

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
`MOI.RotatedSecondOrderCone`.

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

# Deprecation for JuMP v0.18 -> JuMP v0.19 transition
function LinearAlgebra.norm(::AbstractVector{<:AbstractJuMPScalar})
    return error(
        "JuMP no longer performs automatic transformation of `norm()` ",
        "expressions into second-order cone constraints. They should now ",
        "be expressed using the SecondOrderCone() set. For example, ",
        "`@constraint(model, norm(x) <= t)` should now be written as ",
        "`@constraint(model, [t; x] in SecondOrderCone())`",
    )
end

"""
    SOS1

SOS1 (Special Ordered Sets type 1) object than can be used to constrain a
vector `x` to a set where at most 1 variable can take a non-zero value, all
others being at 0.
The `weights`, when specified, induce an ordering of the variables; as such,
they should be unique values. The *k*th element in the set corresponds to the
*k*th weight in `weights`. See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm)
for a description of SOS constraints and their potential uses.
This is a shortcut for the `MathOptInterface.SOS1` set.
"""
struct SOS1 <: AbstractVectorSet
    weights::Vector{Float64}
    function SOS1(weights::AbstractVector = Float64[])
        return new(convert(Vector{Float64}, weights))
    end
end
function moi_set(set::SOS1, dim::Int)
    if length(set.weights) == 0
        return MOI.SOS1(collect(1:1.0:dim))
    elseif length(set.weights) == dim
        return MOI.SOS1(set.weights)
    else
        error("Weight vector in SOS1 is not of length $(dim).")
    end
end

"""
    SOS2

SOS1 (Special Ordered Sets type 2) object than can be used to constrain a
vector `x` to a set where at most 2 variables can take a non-zero value, all
others being at 0. In addition, if two are non-zero these must be consecutive
in their ordering.
The `weights` induce an ordering of the variables; as such, they should be unique
values. The *k*th element in the set corresponds to the *k*th weight in `weights`.
See [here](http://lpsolve.sourceforge.net/5.5/SOS.htm) for a description of SOS
constraints and their potential uses.
This is a shortcut for the `MathOptInterface.SOS2` set.
"""
struct SOS2 <: AbstractVectorSet
    weights::Vector{Float64}
    function SOS2(weights::AbstractVector = Float64[])
        return new(convert(Vector{Float64}, weights))
    end
end
function moi_set(set::SOS2, dim::Int)
    if length(set.weights) == 0
        return MOI.SOS2(collect(1:1.0:dim))
    elseif length(set.weights) == dim
        return MOI.SOS2(set.weights)
    else
        error("Weight vector in SOS2 is not of length $(dim).")
    end
end
