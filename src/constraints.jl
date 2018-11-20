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
struct ConstraintRef{M <: AbstractModel, C, Shape <: AbstractShape}
    model::M
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
name(cr::ConstraintRef{Model,<:MOICON}) = MOI.get(cr.model, MOI.ConstraintName(), cr)

set_name(cr::ConstraintRef{Model,<:MOICON}, s::String) = MOI.set(cr.model, MOI.ConstraintName(), cr, s)

"""
    delete(model::Model, constraint_ref::ConstraintRef)

Delete the constraint associated with `constraint_ref` from the model `model`.
"""
function delete(model::Model, constraint_ref::ConstraintRef{Model})
    if model !== constraint_ref.model
        error("The constraint reference you are trying to delete does not " *
              "belong to the model.")
    end
    MOI.delete(backend(model), index(constraint_ref))
end

"""
    is_valid(model::Model, constraint_ref::ConstraintRef{Model})

Return `true` if `constraint_ref` refers to a valid constraint in `model`.
"""
function is_valid(model::Model, constraint_ref::ConstraintRef{Model})
    return (model === constraint_ref.model &&
            MOI.is_valid(backend(model), constraint_ref.index))
end

#############################################################################
# AbstractConstraint
# Abstract base type for all constraint types
abstract type AbstractConstraint end

"""
    jump_function(constraint::JuMP.AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `JuMP.AbstractJuMPScalar` or `Vector{JuMP.AbstractJuMPScalar}`.
"""
function jump_function end

"""
    moi_function(constraint::JuMP.AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `MathOptInterface.AbstractFunction`.
"""
function moi_function(constraint::JuMP.AbstractConstraint)
    return moi_function(jump_function(constraint))
end

"""
    moi_set(constraint::AbstractConstraint)

Return the set of the constraint `constraint` in the function-in-set form as a
`MathOptInterface.AbstractSet`.

    moi_set(s::AbstractVectorSet, dim::Int)

Returns the MOI set of dimension `dim` corresponding to the JuMP set `s`.
"""
function moi_set end

struct ScalarConstraint{F <: AbstractJuMPScalar,
                        S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
end

jump_function(constraint::ScalarConstraint) = constraint.func
moi_set(constraint::ScalarConstraint) = constraint.set
shape(::ScalarConstraint) = ScalarShape()
function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractScalarFunction, SetType <: MOI.AbstractScalarSet}
    model = ref.model
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

jump_function(constraint::VectorConstraint) = constraint.func
moi_set(constraint::VectorConstraint) = constraint.set
shape(c::VectorConstraint) = c.shape
function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractVectorFunction, SetType <: MOI.AbstractVectorSet}
    model = ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return VectorConstraint(jump_function(model, f), s, ref.shape)
end

function moi_add_constraint(model::MOI.ModelLike, f::MOI.AbstractFunction,
                            s::MOI.AbstractSet)
    if !MOI.supports_constraint(model, typeof(f), typeof(s))
        if moi_mode(model) == Direct
            bridge_message = "."
        elseif moi_bridge_constraints(model)
            bridge_message = " and there are no bridges that can reformulate it into supported constraints."
        else
            bridge_message = ", try using `bridge_constraints=true` in the `JuMP.Model` constructor if you believe the constraint can be reformulated to constraints supported by the solver."
        end
        error("Constraints of type $(typeof(f))-in-$(typeof(s)) are not supported by the solver" * bridge_message)
    end
    return MOI.add_constraint(model, f, s)
end

"""
    add_constraint(model::Model, c::AbstractConstraint, name::String="")

Add a constraint `c` to `Model model` and sets its name.
"""
function add_constraint(model::Model, c::AbstractConstraint, name::String="")
    # The type of backend(model) is unknown so we directly redirect to another
    # function.
    cindex = moi_add_constraint(backend(model), moi_function(c), moi_set(c))
    cref = ConstraintRef(model, cindex, shape(c))
    if !isempty(name)
        set_name(cref, name)
    end
    return cref
end

"""
    set_coefficient(constraint::ConstraintRef, variable::VariableRef, value)

Set the coefficient of `variable` in the constraint `constraint` to `value`.

Note that prior to this step, JuMP will aggregate multiple terms containing the
same variable. For example, given a constraint `2x + 3x <= 2`,
`JuMP.set_coefficient(c, x, 4)` will create the constraint `4x <= 2`.


```jldoctest; setup = :(using JuMP), filter=r"≤|<="
model = Model()
@variable(model, x)
@constraint(model, con, 2x + 3x <= 2)
JuMP.set_coefficient(con, x, 4)
con

# output

con : 4 x <= 2.0
```
"""
function set_coefficient(constraint::ConstraintRef{Model, MOICON{F, S}},
                         variable, value) where {S, T, F <: Union{
                             MOI.ScalarAffineFunction{T},
                             MOI.ScalarQuadraticFunction{T}}}
    MOI.modify(backend(constraint.model), index(constraint),
        MOI.ScalarCoefficientChange(index(variable), convert(T, value)))
    return
end

"""
    shadow_price(constraint::ConstraintRef)

The change in the objective from an infinitesimal relaxation of the constraint.
This value is computed from [`dual`](@ref) and can be queried only when
`has_duals` is `true` and the objective sense is `MinSense` or `MaxSense`
(not `FeasibilitySense`). For linear constraints, the shadow prices differ at
most in sign from the `dual` value depending on the objective sense.

## Notes

- The function simply translates signs from `dual` and does not validate
  the conditions needed to guarantee the sensitivity interpretation of the
  shadow price. The caller is responsible, e.g., for checking whether the solver
  converged to an optimal primal-dual pair or a proof of infeasibility.
- The computation is based on the current objective sense of the model. If this
  has changed since the last solve, the results will be incorrect.
- Relaxation of equality constraints (and hence the shadow price) is defined
  based on which sense of the equality constraint is active.
"""
function shadow_price(constraint::ConstraintRef{Model, <:MOICON})
    error("The shadow price is not defined or not implemented for this type " *
          "of constraint.")
end

# Internal function.
function shadow_price_less_than_(dual_value, sense::MOI.OptimizationSense)
    # When minimizing, the shadow price is nonpositive and when maximizing the
    # shadow price is nonnegative (because relaxing a constraint can only
    # improve the objective). By MOI convention, a feasible dual on a LessThan
    # set is nonpositive, so we flip the sign when maximizing.
    if sense == MOI.MaxSense
        return -dual_value
    elseif sense == MOI.MinSense
        return dual_value
    else
        error("The shadow price is not available because the objective sense " *
              "$sense is not minimization or maximization.")
    end
end

# Internal function.
function shadow_price_greater_than_(dual_value, sense::MOI.OptimizationSense)
    # By MOI convention, a feasible dual on a GreaterThan set is nonnegative,
    # so we flip the sign when minimizing. (See comment in the method above).
    if sense == MOI.MaxSense
        return dual_value
    elseif sense == MOI.MinSense
        return -dual_value
    else
        error("The shadow price is not available because the objective sense " *
              "$sense is not minimization or maximization.")
    end
end

function shadow_price(constraint::ConstraintRef{Model, MOICON{F, S}}
                      ) where {S <: MOI.LessThan, F}
    model = constraint.model
    if !has_duals(model)
        error("The shadow price is not available because no dual result is " *
              "available.")
    end
    return shadow_price_less_than_(dual(constraint),
                                   objective_sense(model))
end

function shadow_price(constraint::ConstraintRef{Model, MOICON{F, S}}
                      ) where {S <: MOI.GreaterThan, F}
    model = constraint.model
    if !has_duals(model)
        error("The shadow price is not available because no dual result is " *
              "available.")
    end
    return shadow_price_greater_than_(dual(constraint),
                                      objective_sense(model))
end

function shadow_price(constraint::ConstraintRef{Model, MOICON{F, S}}
                      ) where {S <: MOI.EqualTo, F}
    model = constraint.model
    if !has_duals(model)
        error("The shadow price is not available because no dual result is " *
              "available.")
    end
    sense = objective_sense(model)
    dual_val = dual(constraint)
    if dual_val > 0
        # Treat the equality constraint as if it were a GreaterThan constraint.
        return shadow_price_greater_than_(dual_val, sense)
    else
        # Treat the equality constraint as if it were a LessThan constraint.
        return shadow_price_less_than_(dual_val, sense)
    end
end
