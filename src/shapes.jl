#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See https://github.com/jump-dev/JuMP.jl
#############################################################################

"""
    AbstractShape

Abstract vectorizable shape. Given a flat vector form of an object of shape
`shape`, the original object can be obtained by [`reshape_vector`](@ref).
"""
abstract type AbstractShape end

"""
    dual_shape(shape::AbstractShape)::AbstractShape

Returns the shape of the dual space of the space of objects of shape `shape`. By
default, the `dual_shape` of a shape is itself. See the examples section below
for an example for which this is not the case.

## Example

Consider polynomial constraints for which the dual is moment constraints and
moment constraints for which the dual is polynomial constraints. Shapes for
polynomials can be defined as follows:
```julia
struct Polynomial
    coefficients::Vector{Float64}
    monomials::Vector{Monomial}
end
struct PolynomialShape <: AbstractShape
    monomials::Vector{Monomial}
end
JuMP.reshape_vector(x::Vector, shape::PolynomialShape) = Polynomial(x, shape.monomials)
```
and a shape for moments can be defined as follows:
```julia
struct Moments
    coefficients::Vector{Float64}
    monomials::Vector{Monomial}
end
struct MomentsShape <: AbstractShape
    monomials::Vector{Monomial}
end
JuMP.reshape_vector(x::Vector, shape::MomentsShape) = Moments(x, shape.monomials)
```
Then `dual_shape` allows the definition of the shape of the dual of polynomial
and moment constraints:
```julia
dual_shape(shape::PolynomialShape) = MomentsShape(shape.monomials)
dual_shape(shape::MomentsShape) = PolynomialShape(shape.monomials)
```
"""
dual_shape(shape::AbstractShape) = shape

"""
    reshape_set(vectorized_set::MOI.AbstractSet, shape::AbstractShape)

Return a set in its original shape `shape` given its vectorized form
`vectorized_form`.

## Example

Given a [`SymmetricMatrixShape`](@ref) of vectorized form
`[1, 2, 3] in MOI.PositiveSemidefinieConeTriangle(2)`, the
following code returns the set of the original constraint
`Symmetric(Matrix[1 2; 2 3]) in PSDCone()`:

```jldoctest
julia> reshape_set(MOI.PositiveSemidefiniteConeTriangle(2), SymmetricMatrixShape(2))
PSDCone()
```
"""
function reshape_set end

"""
    reshape_vector(vectorized_form::Vector, shape::AbstractShape)

Return an object in its original shape `shape` given its vectorized form
`vectorized_form`.

## Example

Given a [`SymmetricMatrixShape`](@ref) of vectorized form `[1, 2, 3]`, the
following code returns the matrix `Symmetric(Matrix[1 2; 2 3])`:

```jldoctest
julia> reshape_vector([1, 2, 3], SymmetricMatrixShape(2))
2×2 LinearAlgebra.Symmetric{Int64, Matrix{Int64}}:
 1  2
 2  3
```
"""
function reshape_vector end

"""
    shape(c::AbstractConstraint)::AbstractShape

Return the shape of the constraint `c`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> c = @constraint(model, x[2] <= 1);

julia> shape(constraint_object(c))
ScalarShape()

julia> d = @constraint(model, x in SOS1());

julia> shape(constraint_object(d))
VectorShape()
```
"""
function shape end

"""
    ScalarShape()

An [`AbstractShape`](@ref) that represents scalar constraints.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> c = @constraint(model, x[2] <= 1);

julia> shape(constraint_object(c))
ScalarShape()
```
"""
struct ScalarShape <: AbstractShape end

reshape_vector(α, ::ScalarShape) = α

"""
    VectorShape()

An [`AbstractShape`](@ref) that represents vector-valued constraints.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> c = @constraint(model, x in SOS1());

julia> shape(constraint_object(c))
VectorShape()
```
"""
struct VectorShape <: AbstractShape end

reshape_vector(vectorized_form, ::VectorShape) = vectorized_form

vectorize(x, ::VectorShape) = x

"""
    ArrayShape{N}(dims::NTuple{N,Int}) where {N}

An [`AbstractShape`](@ref) that represents array-valued constraints.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2, 1:3]);

julia> c = @constraint(model, x >= 0, Nonnegatives())
[x[1,1]  x[1,2]  x[1,3]
 x[2,1]  x[2,2]  x[2,3]] ∈ Nonnegatives()

julia> shape(constraint_object(c))
ArrayShape{2}((2, 3))
```
"""
struct ArrayShape{N} <: AbstractShape
    dims::NTuple{N,Int}
end

reshape_vector(x, shape::ArrayShape) = reshape(x, shape.dims)
reshape_vector(::Nothing, shape::ArrayShape) = nothing

vectorize(x, ::ArrayShape) = vec(x)
