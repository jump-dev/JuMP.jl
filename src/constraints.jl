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
`shape`, the original object can be obtained by [`reshape_result`](@ref).
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
struct PolynomialShape <: AbstractShape
    monomials::Vector{Monomial}
end
JuMP.reshape_result(x::Vector, shape::PolynomialShape) = Polynomial(x, shape.monomials)
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
JuMP.reshape_result(x::Vector, shape::MomentsShape) = Moments(x, shape.monomials)
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
    reshape_result(vectorized_form::Vector, shape::AbstractShape)

Return an object in its original shape `shape` given its vectorized form
`vectorized_form`.

## Examples

Given a [`SymmetricMatrixShape`](@ref) of vectorized form `[1, 2, 3]`, the
following code returns the matrix `Symmetric(Matrix[1 2; 2 3])`:
```julia
reshape_result([1, 2, 3], SymmetricMatrixShape(2))
```
"""
function reshape_result end

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
reshape_result(α, ::ScalarShape) = α

"""
    VectorShape

Vector for which the vectorized form corresponds exactly to the vector given.
"""
struct VectorShape <: AbstractShape end
reshape_result(vectorized_form, ::VectorShape) = vectorized_form

"""
    ConstraintRef

Holds a reference to the model and the corresponding MOI.ConstraintIndex.
"""
struct ConstraintRef{M <: AbstractModel, C, Shape <: AbstractShape}
    model::M
    index::C
    shape::Shape
end

"""
    struct ConstraintNotOwned{C <: ConstraintRef} <: Exception
        constraint_ref::C
    end

The constraint `constraint_ref` was used in a model different to
`owner_model(constraint_ref)`.
"""
struct ConstraintNotOwned{C <: ConstraintRef} <: Exception
    constraint_ref::C
end

"""
    owner_model(cref::ConstraintRef)

Returns the model to which `cref` belongs.
"""
owner_model(cref::ConstraintRef) = cref.model

"""
    check_belongs_to_model(cref::ConstraintRef, model::AbstractModel)

Throw `ConstraintNotOwned` if `owner_model(cref)` is not `model`.
"""
function check_belongs_to_model(cref::ConstraintRef, model::AbstractModel)
    if owner_model(cref) !== model
        throw(ConstraintNotOwned(cref))
    end
end

Base.broadcastable(cref::ConstraintRef) = Ref(cref)

"""
    name(v::ConstraintRef)

Get a constraint's name attribute.
"""
function name(cr::ConstraintRef{Model,<:_MOICON})
    return MOI.get(cr.model, MOI.ConstraintName(), cr)::String
end

"""
    set_name(v::ConstraintRef, s::AbstractString)

Set a constraint's name attribute.
"""
set_name(cr::ConstraintRef{Model,<:_MOICON}, s::String) = MOI.set(cr.model, MOI.ConstraintName(), cr, s)

"""
    constraint_by_name(model::AbstractModel,
                       name::String)::Union{ConstraintRef, Nothing}

Returns the reference of the constraint with name attribute `name` or `Nothing`
if no constraint has this name attribute. Throws an error if several
constraints have `name` as their name attribute.

    constraint_by_name(model::AbstractModel,
                       name::String,
                       F::Type{<:Union{AbstractJuMPScalar,
                                       Vector{<:AbstractJuMPScalar},
                                       MOI.AbstactFunction}},
                       S::Type{<:MOI.AbstractSet})::Union{ConstraintRef, Nothing}

Similar to the method above, except that it throws an error if the constraint is
not an `F`-in-`S` contraint where `F` is either the JuMP or MOI type of the
function, and `S` is the MOI type of the set. This method is recommended if you
know the type of the function and set since its returned type can be inferred
while for the method above (i.e. without `F` and `S`), the exact return type of
the constraint index cannot be inferred.

```jldoctest objective_function; setup = :(using JuMP), filter = r"Stacktrace:.*"s
julia> using JuMP

julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x)
x

julia> @constraint(model, con, x^2 == 1)
con : x² = 1.0

julia> constraint_by_name(model, "kon")

julia> constraint_by_name(model, "con")
con : x² = 1.0

julia> constraint_by_name(model, "con", AffExpr, MOI.EqualTo{Float64})

julia> constraint_by_name(model, "con", QuadExpr, MOI.EqualTo{Float64})
con : x² = 1.0
```
"""
function constraint_by_name end

function constraint_by_name(model::Model, name::String)
    index = MOI.get(backend(model), MOI.ConstraintIndex, name)
    if index isa Nothing
        return nothing
    else
        return constraint_ref_with_index(model, index)
    end
end

function constraint_by_name(model::Model, name::String,
                            F::Type{<:MOI.AbstractFunction},
                            S::Type{<:MOI.AbstractSet})
    index = MOI.get(backend(model), MOI.ConstraintIndex{F, S}, name)
    if index isa Nothing
        return nothing
    else
        return constraint_ref_with_index(model, index)
    end
end
function constraint_by_name(model::Model, name::String,
                            F::Type{<:Union{ScalarType,
                                            Vector{ScalarType}}},
                            S::Type) where ScalarType <: AbstractJuMPScalar
    return constraint_by_name(model, name, moi_function_type(F), S)
end

# Creates a ConstraintRef with default shape
function constraint_ref_with_index(model::AbstractModel,
                        index::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction,
                                                   <:MOI.AbstractScalarSet})
    return ConstraintRef(model, index, ScalarShape())
end
function constraint_ref_with_index(model::AbstractModel,
                        index::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction,
                                                   <:MOI.AbstractVectorSet})
    return ConstraintRef(model, index, VectorShape())
end

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
"""
    abstract type AbstractConstraint

An abstract base type for all constraint types. `AbstractConstraint`s store the
function and set directly, unlike [`ConstraintRef`](@ref)s that are merely
references to constraints stored in a model. `AbstractConstraint`s do not need
to be attached to a model.
"""
abstract type AbstractConstraint end

"""
    struct BridgeableConstraint{C, B} <: AbstractConstraint
        constraint::C
        bridge_type::B
    end

Constraint `constraint` that can be bridged by the bridge of type `bridge_type`.
Adding this constraint to a model is equivalent to
```julia
add_bridge(model, bridge_type)
add_constraint(model, constraint)
```

## Examples

Given a new scalar set type `CustomSet` with a bridge `CustomBridge` that can
bridge `F`-in-`CustomSet` constraints, when the user does
```julia
model = Model()
@variable(model, x)
@constraint(model, x + 1 in CustomSet())
optimize!(model)
```
with an optimizer that does not support `F`-in-`CustomSet` constraints, the
constraint will not be bridge unless he manually calls `add_bridge(model,
CustomBridge)`. In order to automatically add the `CustomBridge` to any model to
which an `F`-in-`CustomSet` is added, simply add the following method:
```julia
function JuMP.build_constraint(_error::Function, func::AbstractJuMPScalar,
                               set::CustomSet)
    constraint = ScalarConstraint(func, set)
    return JuMP.BridgeableConstraint(constraint, CustomBridge)
end
```

### Note

JuMP extensions should extend `JuMP.build_constraint` only if they also defined
`CustomSet`, for three
reasons:
1. It is problematic if multiple extensions overload the same JuMP method.
2. A missing method will not inform the users that they forgot to load the
   extension module defining the `build_constraint` method.
3. Defining a method where neither the function nor any of the argument types
   are defined in the package is called [*type piracy*](https://docs.julialang.org/en/v1/manual/style-guide/index.html#Avoid-type-piracy-1)
   and is discouraged in the Julia style guide.
```
"""
struct BridgeableConstraint{C, B} <: AbstractConstraint
    constraint::C
    bridge_type::B
end

function add_constraint(model::Model, c::BridgeableConstraint, name::String="")
    add_bridge(model, c.bridge_type)
    return add_constraint(model, c.constraint, name)
end

"""
    jump_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `AbstractJuMPScalar` or `Vector{AbstractJuMPScalar}`.
"""
function jump_function end

"""
    moi_function(constraint::AbstractConstraint)

Return the function of the constraint `constraint` in the function-in-set form
as a `MathOptInterface.AbstractFunction`.
"""
function moi_function(constraint::AbstractConstraint)
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

"""
    constraint_object(ref::ConstraintRef)

Return the underlying constraint data for the constraint referenced by `ref`.
"""
function constraint_object end

"""
    struct ScalarConstraint

The data for a scalar constraint. The `func` field containts a JuMP object
representing the function and the `set` field contains the MOI set.
See also the [documentation](@ref Constraints) on JuMP's representation of
constraints for more background.
"""
struct ScalarConstraint{F <: AbstractJuMPScalar,
                        S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::F
    set::S
end

jump_function(constraint::ScalarConstraint) = constraint.func
moi_set(constraint::ScalarConstraint) = constraint.set
shape(::ScalarConstraint) = ScalarShape()

function constraint_object(ref::ConstraintRef{Model, _MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractScalarFunction, SetType <: MOI.AbstractScalarSet}
    model = ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return ScalarConstraint(jump_function(model, f), s)
end
function check_belongs_to_model(c::ScalarConstraint, model)
    check_belongs_to_model(c.func, model)
end

"""
    struct VectorConstraint

The data for a vector constraint. The `func` field containts a JuMP object
representing the function and the `set` field contains the MOI set. The
`shape` field contains an [`AbstractShape`](@ref) matching the form in which
the constraint was constructed (e.g., by using matrices or flat vectors).
See also the [documentation](@ref Constraints) on JuMP's representation of
constraints.
"""
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
function constraint_object(ref::ConstraintRef{Model, _MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.AbstractVectorFunction, SetType <: MOI.AbstractVectorSet}
    model = ref.model
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return VectorConstraint(jump_function(model, f), s, ref.shape)
end
function check_belongs_to_model(c::VectorConstraint, model)
    for func in c.func
        check_belongs_to_model(func, model)
    end
end

function moi_add_constraint(model::MOI.ModelLike, f::MOI.AbstractFunction,
                            s::MOI.AbstractSet)
    if !MOI.supports_constraint(model, typeof(f), typeof(s))
        if moi_mode(model) == DIRECT
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
    check_belongs_to_model(c, model)
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
`set_coefficient(c, x, 4)` will create the constraint `4x <= 2`.


```jldoctest; setup = :(using JuMP), filter=r"≤|<="
model = Model()
@variable(model, x)
@constraint(model, con, 2x + 3x <= 2)
set_coefficient(con, x, 4)
con

# output

con : 4 x <= 2.0
```
"""
function set_coefficient(constraint::ConstraintRef{Model, _MOICON{F, S}},
                         variable, value) where {S, T, F <: Union{
                             MOI.ScalarAffineFunction{T},
                             MOI.ScalarQuadraticFunction{T}}}
    MOI.modify(backend(constraint.model), index(constraint),
        MOI.ScalarCoefficientChange(index(variable), convert(T, value)))
    return
end

"""
    value(cref::ConstraintRef)

Get the primal value of this constraint in the result returned by a solver. That
is, if `cref` is the reference of a constraint `func`-in-`set`, it returns the
value of `func` evaluated at the value of the variables (given by
[`value(::VariableRef)`](@ref)).
Use [`has_values`](@ref) to check if a result exists before asking for values.

## Note

For scalar contraints, the constant is moved to the `set` so it is not taken
into account in the primal value of the constraint. For instance, the constraint
`@constraint(model, 2x + 3y + 1 == 5)` is transformed into
`2x + 3y`-in-`MOI.EqualTo(4)` so the value returned by this function is the
evaluation of `2x + 3y`.
```
"""
function value(cref::ConstraintRef{Model, <:_MOICON})
    return reshape_result(_constraint_primal(cref), cref.shape)
end

# Returns the value of MOI.ConstraintPrimal in a type-stable way
function _constraint_primal(
    cref::ConstraintRef{Model, <:_MOICON{
        <:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}})::Float64
    return MOI.get(cref.model, MOI.ConstraintPrimal(), cref)
end
function _constraint_primal(
    cref::ConstraintRef{Model, <:_MOICON{
        <:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}})::Vector{Float64}
    return MOI.get(cref.model, MOI.ConstraintPrimal(), cref)
end

"""
    has_duals(model::Model)

Return true if the solver has a dual solution available to query, otherwise
return false.

See also [`dual`](@ref) and [`shadow_price`](@ref).
"""
has_duals(model::Model) = dual_status(model) != MOI.NO_SOLUTION

"""
    dual(cref::ConstraintRef)

Get the dual value of this constraint in the result returned by a solver.
Use `has_dual` to check if a result exists before asking for values.
See also [`shadow_price`](@ref).
"""
function dual(cref::ConstraintRef{Model, <:_MOICON})
    return reshape_result(_constraint_dual(cref), dual_shape(cref.shape))
end

# Returns the value of MOI.ConstraintPrimal in a type-stable way
function _constraint_dual(
    cref::ConstraintRef{Model, <:_MOICON{
        <:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}})::Float64
    return MOI.get(cref.model, MOI.ConstraintDual(), cref)
end
function _constraint_dual(
    cref::ConstraintRef{Model, <:_MOICON{
        <:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}})::Vector{Float64}
    return MOI.get(cref.model, MOI.ConstraintDual(), cref)
end


"""
    shadow_price(constraint::ConstraintRef)

The change in the objective from an infinitesimal relaxation of the constraint.
This value is computed from [`dual`](@ref) and can be queried only when
`has_duals` is `true` and the objective sense is `MIN_SENSE` or `MAX_SENSE`
(not `FEASIBILITY_SENSE`). For linear constraints, the shadow prices differ at
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
function shadow_price(constraint::ConstraintRef{Model, <:_MOICON})
    error("The shadow price is not defined or not implemented for this type " *
          "of constraint.")
end

# Internal function.
function shadow_price_less_than_(dual_value, sense::MOI.OptimizationSense)
    # When minimizing, the shadow price is nonpositive and when maximizing the
    # shadow price is nonnegative (because relaxing a constraint can only
    # improve the objective). By MOI convention, a feasible dual on a LessThan
    # set is nonpositive, so we flip the sign when maximizing.
    if sense == MOI.MAX_SENSE
        return -dual_value
    elseif sense == MOI.MIN_SENSE
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
    if sense == MOI.MAX_SENSE
        return dual_value
    elseif sense == MOI.MIN_SENSE
        return -dual_value
    else
        error("The shadow price is not available because the objective sense " *
              "$sense is not minimization or maximization.")
    end
end

function shadow_price(constraint::ConstraintRef{Model, _MOICON{F, S}}
                      ) where {S <: MOI.LessThan, F}
    model = constraint.model
    if !has_duals(model)
        error("The shadow price is not available because no dual result is " *
              "available.")
    end
    return shadow_price_less_than_(dual(constraint),
                                   objective_sense(model))
end

function shadow_price(constraint::ConstraintRef{Model, _MOICON{F, S}}
                      ) where {S <: MOI.GreaterThan, F}
    model = constraint.model
    if !has_duals(model)
        error("The shadow price is not available because no dual result is " *
              "available.")
    end
    return shadow_price_greater_than_(dual(constraint),
                                      objective_sense(model))
end

function shadow_price(constraint::ConstraintRef{Model, _MOICON{F, S}}
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

function _error_if_not_concrete_type(t)
    if !isconcretetype(t)
        error("`$t` is not a concrete type. Did you miss a type parameter?")
    end
    return
end

"""
    all_constraints(model::Model, function_type, set_type)::Vector{VariableRef}

Return a list of all constraints currently in the model where the function
has type `function_type` and the set has type `set_type`. The constraints are
ordered by creation time.

See also [`list_of_constraint_types`](@ref).

# Example
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Bin);

julia> @constraint(model, 2x <= 1);

julia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}},ScalarShape},1}:
 x ≥ 0.0

julia> all_constraints(model, VariableRef, MOI.ZeroOne)
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.ZeroOne},ScalarShape},1}:
 x binary

julia> all_constraints(model, AffExpr, MOI.LessThan{Float64})
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 2 x ≤ 1.0
```
"""
function all_constraints(model::Model,
                         function_type::Type{<:AbstractJuMPScalar},
                         set_type::Type{<:MOI.AbstractSet})
    _error_if_not_concrete_type(function_type)
    _error_if_not_concrete_type(set_type)
    # TODO: Support JuMP's set helpers like SecondOrderCone().
    f_type = moi_function_type(function_type)
    constraint_ref_type = ConstraintRef{Model, _MOICON{f_type, set_type},
                                        ScalarShape}
    result = constraint_ref_type[]
    for idx in MOI.get(model, MOI.ListOfConstraintIndices{f_type, set_type}())
        push!(result, ConstraintRef(model, idx, ScalarShape()))
    end
    return result
end

# TODO: Support vector function types. This is blocked by not having the shape
# information available.

"""
    list_of_constraint_types(model::Model)

Return a list of tuples of the form `(F, S)` where `F` is a JuMP function type
and `S` is an MOI set type such that `all_constraints(model, F, S)` returns
a nonempty list.

# Example
```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0, Bin);

julia> @constraint(model, 2x <= 1);

julia> list_of_constraint_types(model)
3-element Array{Tuple{DataType,DataType},1}:
 (VariableRef, MathOptInterface.ZeroOne)
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
"""
function list_of_constraint_types(model::Model)
    list = MOI.get(
        model, MOI.ListOfConstraints())::Vector{Tuple{DataType, DataType}}
    return Tuple{DataType, DataType}[(jump_function_type(model, f), s)
                                     for (f,s) in list]
end
