```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" ∈ | in "]
```

# Constraints

This page explains how to write various types of constraints in JuMP. Before
reading further, please make sure you are familiar with JuMP models, and JuMP
[Variables](@ref). If you want to add nonlinear constraints, read
[Nonlinear Modeling](@ref) instead.

JuMP is based on the MathOptInterface API. Because of this, JuMP thinks of a
constraint as a *function* belonging to a *set*. For example, instead of
thinking  about a constraint ``a^\top x \le b`` as a *less-than-or-equal-to*
constraint, JuMP thinks about this as the *scalar affine* function ``a^\top x``
belonging to the *less-than* set ``(-\infty, b]``. Thus, instead of a
*less-than-or-equal-to* constraint, we consider this constraint to be a  *scalar
affine -in- less than* constraint. More generally, we use the shorthand
*function-in-set* to refer to constraints composed of different types of
functions and sets. In the rest of this page, we will introduce the different
types of functions and sets that JuMP knows about as needed. You can also more
details about this *function-in-set* concept in the MathOptInterface
documentation.

## The `@constraint` macro

Constraints are added to a JuMP model using the [`@constraint`](@ref) macro. It
is similar to the [`@variable`](@ref) macro. Here is an example of how to add
the constraint ``2x \le 1`` to a JuMP model:
```jldoctest con1; setup = :(model = JuMP.Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0
```
Wasn't that easy! Let's unpack what happened, because just like
[`@variable`](@ref) there are a few subtle things going on.
 1. The mathematical constraint ``2x \le 1`` was added to the model.
 2. A JuMP called `con` was created that acts as a reference to the
    constraint.
 3. JuMP set the *name* of the constraint to `"con"`.

Just like the Julia variables created in [`@variable`](@ref), `con` can be bound
to a different value. For example:
```jldoctest con1
julia> con
con : 2 x <= 1.0

julia> con = 1
1

julia> con
1
```
However, the reference can be retrieved by querying the model using the symbolic
name:
```jldoctest con1
julia> con = model[:con]
con : 2 x <= 1.0

julia> con
con : 2 x <= 1.0
```
Because the named variables and constraints are stored in the same namespace,
creating a constraint with the same name as a variable or an existing constraint
will result in an error. To overcome this limitation, it is possible to create
anonymous constraints, just like it is possible to create
[Anonymous JuMP variables](@ref). This is done by dropping the second argument
to [`@constraint`](@ref):
```jldoctest con1
julia> con = @constraint(model, 2x <= 1)
2 x <= 1.0
```

It is also possible use different comparison operators (e.g., `>=` and `==`) to
create the following types of constraints:
```jldoctest con1
julia> @constraint(model, 2x >= 1)
2 x >= 1.0

julia> @constraint(model, 2x == 1)
2 x == 1.0

julia> @constraint(model, 1 <= 2x <= 3)
2 x ∈ [1.0, 3.0]
```

## Constraint containers

So far, we've adding constraints one-by-one. However, just like
[Variable containers](@ref), JuMP provides a mechanism for building groups of
constraints in *containers*. Three types of constraint containers are supported:
`Array`s, `JuMPArray`s, and `Dict`ionaries. We explain each of these in the
following.

### [Arrays](@id constraint_arrays)

The syntax for adding `Array`s of JuMP constraints is very similar to the
[syntax for creating `Array`s of JuMP variables](@ref Arrays):
```jldoctest constraint_arrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:3], i * x <= i + 1)
3-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 con[1] : x <= 2.0
 con[2] : 2 x <= 3.0
 con[3] : 3 x <= 4.0
```
These arrays can be accessed and sliced like normal Julia arrays:
```jldoctest constraint_arrays
julia> con[1]
con[1] : x <= 2.0

julia> con[2:3]
2-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 con[2] : 2 x <= 3.0
 con[3] : 3 x <= 4.0
```

Anonymous containers can also be constructed by dropping the name (e.g. `con`)
before the square brackets:
```jldoctest constraint_arrays
julia> @constraint(model, [i = 1:2], i * x <= i + 1)
2-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 x <= 2.0
 2 x <= 3.0
```

Just like [`@variable`](@ref), JuMP will form an `Array` of constraints when it
can determine at compile time that the indices are one-based integer ranges.
Therefore `con[1:b]` will work, but `con[a:b]` will throw an error. If JuMP
cannot determine that the indices are one-based integer ranges, JuMP will
constraint a `JuMPArray` instead.

### JuMPArrays

The syntax for constructing a `JuMPArray` of constraints is very similar to the
[syntax for constructing](@ref variable_jump_arrays) a JuMPArray of variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 2:3], i * x <= j + 1)
2-dimensional JuMPArray{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,2,...} with index sets:
    Dimension 1, 1:2
    Dimension 2, 2:3
And data, a 2×2 Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,2}:
 con[1,2] : x <= 3.0    con[1,3] : x <= 4.0
 con[2,2] : 2 x <= 3.0  con[2,3] : 2 x <= 4.0
```

### Dictionaries

The syntax for constructing a dictionary of constraints is very similar to the
[syntax for constructing](@ref variable_dictionaries) a dictionary of variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)
Dict{Any,ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C} with 2 entries:
  (1, 2) => con[1,2] : x <= 3.0
  (2, 1) => con[2,1] : 2 x <= 2.0
```

### Forcing the container type

When creating a container of constraints, JuMP will attempt to choose the
tightest container type that can store the constraints. Thus, it will prefer
to create an Array before a JuMPArray, and a JuMPArray before a dictionary.
However, because this happens at compile time, it does not always make the best
choice. To illustrate this, consider the following example:
```jldoctest con_force_container; setup=:(model=Model(); @variable(model, x))
julia> A = 1:2
1:2

julia> @constraint(model, con[i=A], i * x <= i + 1)
1-dimensional JuMPArray{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 con[1] : x <= 2.0
 con[2] : 2 x <= 3.0
```

Since the value (and type) of `A` is unknown at compile time, JuMP is unable to
infer that `A` is a one-based integer range. Therefore, JuMP creates a
`JuMPArray`, even though it could store these constraints in a standard
one-dimensional `Array`.

We can share our knowledge that it is possible to store these constraints as
an array by setting the `container` keyword:
```jldoctest con_force_container
julia> @constraint(model, con2[i=A], i * x <= i + 1, container=Array)
2-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 con2[1] : x <= 2.0
 con2[2] : 2 x <= 3.0
```
JuMP now creates an Array of constraints instead of a JuMPArray. Note that
choosing an invalid container type will throw an error.

## Quadratic constraints

`@constraint(model, x^2 + y^2 <= z^2)`

## Constraints of a collection of variables

Cones, SOS1, SO2, etc



DRAFT: Describe how constraints are represented (link to MOI docs). Constraints
are very similar to variables in (1) how names work (2) how attributes work, and
(3) the macro syntax for constructing them. They're a bit different because
they're parameterized by function-set type. Describe constraints vs.
`ConstraintRefs`. Describe `JuMP.constraint_object`. How to delete constraints.
How to modify constraints by setting attributes and `MOI.modifyconstraint!`.
Describe semidefinite constraints and symmetry handling. Refer to NLP docs for
nonlinear constraints.

## Sets

As mentioned in the documentation of the [`@constraint`](@ref) and
[`@SDconstraint`](@ref) macros, the following sets can be used to create
constraints in addition to [any MOI set](http://www.juliaopt.org/MathOptInterface.jl/v0.6.2/apireference.html#Sets-1).

## Constraint modifications

`JuMP.fix`

`JuMP.set_coefficient(constraint, variable, value)`

## Duals

## Reference

```@docs
@constraint
@SDconstraint
SecondOrderCone
RotatedSecondOrderCone
PSDCone
JuMP.set_coefficient
JuMP.dual
JuMP.shadow_price
```

## Constructing constraints without adding them to the model

For advanced use cases.

```@docs
JuMP.@build_constraint
```
