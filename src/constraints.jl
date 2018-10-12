#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

"""
    AbstractShape

Abstract vectorizable shape. Given a flat vector form of an object of shape
`shape`, the original object can be obtained by [`reshape`](@ref).
"""
abstract type AbstractShape end

"""
    dual_shape(shape::AbstractShape)::AbstractShape

Returns the shape of the dual space of the space of objects of shape `shape`. By
default, the `dual_shape` of a shape is itself. See the examples section below
for an example for which this is not the case.

## Examples

Consider polynomial constraints for which the dual is moment constraints and
moment constraints for which the dual is polynomial constraints. Shapes for
polynomials can be defined as follows:
```julia
struct Polynomial
    coefficients::Vector{Float64}
    monomials::Vector{Monomial}
end
struct PolynomialShape <: JuMP.AbstractShape
    monomials::Vector{Monomial}
end
JuMP.reshape(x::Vector, shape::PolynomialShape) = Polynomial(x, shape.monomials)
```
and a shape for moments can be defined as follows:
```julia
struct Moments
    coefficients::Vector{Float64}
    monomials::Vector{Monomial}
end
struct MomentsShape <: JuMP.AbstractShape
    monomials::Vector{Monomial}
end
JuMP.reshape(x::Vector, shape::MomentsShape) = Moments(x, shape.monomials)
```
The `dual_shape` allows to define the shape of the dual of polynomial and moment
constraints:
```julia
dual_shape(shape::PolynomialShape) = MomentsShape(shape.monomials)
dual_shape(shape::MomentsShape) = PolynomialShape(shape.monomials)
```
"""
dual_shape(shape::AbstractShape) = shape

"""
    reshape(vectorized_form::Vector, shape::AbstractShape)

Return an object in its original shape `shape` given its vectorized form
`vectorized_form`.

## Examples

Given a [`SymmetricMatrixShape`](@ref) of vectorized form `[1, 2, 3]`, the
following code returns the matrix `Symmetric(Matrix[1 2; 2 3])`:
```julia
reshape([1, 2, 3], SymmetricMatrixShape(2))
```
"""
function reshape end

"""
    shape(c::AbstractConstraint)::AbstractShape

Return the shape of the constraint `c`.
"""
function shape end

"""
    ScalarShape

Shape of scalar constraints.
"""
struct ScalarShape <: AbstractShape end
reshape(α, ::ScalarShape) = α

"""
    VectorShape

Vector for which the vectorized form corresponds exactly to the vector given.
"""
struct VectorShape <: AbstractShape end
reshape(vectorized_form, ::VectorShape) = vectorized_form

##########################################################################
# Constraint Reference
# Holds the index of a constraint in a Model.
# TODO: Rename "m" field (breaks style guidelines).
struct ConstraintRef{M <: AbstractModel, C, Shape <: AbstractShape}
    m::M
    index::C
    shape::Shape
end

if VERSION >= v"0.7-"
    Base.broadcastable(cref::ConstraintRef) = Ref(cref)
end

"""
    name(v::ConstraintRef)

Get a constraint's name.
"""
name(cr::ConstraintRef{Model,<:MOICON}) = MOI.get(cr.m, MOI.ConstraintName(), cr)

set_name(cr::ConstraintRef{Model,<:MOICON}, s::String) = MOI.set(cr.m, MOI.ConstraintName(), cr, s)

"""
    delete(model::Model, constraint_ref::ConstraintRef)

Delete the constraint associated with `constraint_ref` from the model `model`.
"""
function delete(model::Model, constraint_ref::ConstraintRef{Model})
    if model !== constraint_ref.m
        error("The constraint reference you are trying to delete does not " *
              "belong to the model.")
    end
    MOI.delete(model.moi_backend, index(constraint_ref))
end

"""
    is_valid(model::Model, constraint_ref::ConstraintRef{Model})

Return `true` if `constraint_ref` refers to a valid constraint in `model`.
"""
function is_valid(model::Model, constraint_ref::ConstraintRef{Model})
    return (model === constraint_ref.m &&
            MOI.is_valid(model.moi_backend, constraint_ref.index))
end

#############################################################################
# AbstractConstraint
# Abstract base type for all constraint types
abstract type AbstractConstraint end

struct ScalarConstraint{F <: AbstractJuMPScalar,
                        S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
end

moi_function_and_set(c::ScalarConstraint) = (moi_function(c.func), c.set)
shape(::ScalarConstraint) = ScalarShape()
function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractScalarFunction, SetType <: MOI.AbstractScalarSet}
    model = ref.m
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return ScalarConstraint(jump_function(model, f), s)
end

struct VectorConstraint{F <: AbstractJuMPScalar,
                        S <: MOI.AbstractVectorSet,
                        Shape <: AbstractShape} <: AbstractConstraint
    func::Vector{F}
    set::S
    shape::Shape
end
function VectorConstraint(func::Vector{<:AbstractJuMPScalar},
                          set::MOI.AbstractVectorSet)
    VectorConstraint(func, set, VectorShape())
end

moi_function_and_set(c::VectorConstraint) = (moi_function(c.func), c.set)
shape(c::VectorConstraint) = c.shape
function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractVectorFunction, SetType <: MOI.AbstractVectorSet}
    model = ref.m
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return VectorConstraint(jump_function(model, f), s, ref.shape)
end

"""
    add_constraint(m::Model, c::AbstractConstraint, name::String="")

Add a constraint `c` to `Model m` and sets its name.
"""
function add_constraint(m::Model, c::AbstractConstraint, name::String="")
    f, s = moi_function_and_set(c)
    if !MOI.supports_constraint(m.moi_backend, typeof(f), typeof(s))
        if m.moi_backend isa MOI.Bridges.LazyBridgeOptimizer
            bridge_message = " and there are no bridges that can reformulate it into supported constraints."
        else
            bridge_message = ", try using `bridge_constraints=true` in the `JuMP.Model` constructor if you believe the constraint can be reformulated to constraints supported by the solver."
        end
        error("Constraints of type $(typeof(f))-in-$(typeof(s)) are not supported by the solver" * bridge_message)
    end
    cindex = MOI.add_constraint(m.moi_backend, f, s)
    cref = ConstraintRef(m, cindex, shape(c))
    if !isempty(name)
        set_name(cref, name)
    end
    return cref
end

"""
    set_rhs(constraint, value)

Set the right-hand side of a constraint. Supported constraints are
 - ScalarAffineFunction -in- LessThan
 - ScalarAffineFunction -in- GreaterThan
 - ScalarAffineFunction -in- EqualTo.

Interval sets are *not* supported.

Note that prior to this step, JuMP will move all variable terms onto the
left-hand side, and all constant terms onto the right-hand side. For example,
given a constraint `2x + 1 <= 2`, `JuMP.set_constant(c, 3)` will create the
constraint `2x <= 3`, not `2x + 1 <= 3`.
"""
function set_rhs(constraint::ConstraintRef{Model, MOICON{F, SetType}},
                 value) where {F <: MOI.ScalarAffineFunction,
                               SetType <: Union{MOI.LessThan, MOI.GreaterThan,
                                                MOI.EqualTo}}
    MOI.set(constraint.m.moi_backend, MOI.ConstraintSet(), index(constraint),
            SetType(value))
end

"""
    set_coefficient(constraint::ConstraitRef, variable::VariableRef, value)

Set the coefficient of `variable` in the constraint `constraint` to `value`.

Note that prior to this step, JuMP will aggregate multiple terms containing the
same variable. For example, given a constraint `2x + 3x <= 2`,
`JuMP.set_coefficient(c, x, 4)` will create the constraint `4x <= 2`.
"""
function set_coefficient(constraint::ConstraintRef{Model, MOICON{F}},
                         variable, value) where {T, F <: Union{
                             MOI.ScalarAffineFunction{T},
                             MOI.ScalarQuadraticFunction{T}}}
    MOI.modify(con_ref.m.moi_backend, index(con_ref),
        MOI.ScalarCoefficientChange(index(variable), convert(T, value)))
end
