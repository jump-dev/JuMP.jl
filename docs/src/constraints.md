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
[Variables](@ref).

JuMP is based on the MathOptInterface API. Because of this, JuMP thinks of a
constraint as a *function* belonging to a *set*. For example, instead of
thinking  about a constraint $a^\top x \le b$ as a *less-than-or-equal-to*
constraint, JuMP thinks about this as the *scalar affine* function $a^\top x$
belonging to the *less-than* set $(-\infty, b]$. Thus, instead of a
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
the constraint `2x <= 1` to a JuMP model:
```jldoctest con1; setup = :(model = JuMP.Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0
```
Wasn't that easy! Let's unpack what happened, because just like
[`@variable`](@ref) there are a few subtle things going on.
 1. The mathematical constraint $2x \le 1$ was added to the model.
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

### Arrays

```jldoctest constraint_arrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2], i * x <= i + 1)
2-element Array{ConstraintRef{Model,C,Shape} where Shape<:JuMP.AbstractShape where C,1}:
 con[1] : x <= 2.0
 con[2] : 2 x <= 3.0
```

### JuMPArrays

### Dictionaries

DRAFT: Describe how constraints are represented (link to MOI docs). Constraints
are very similar to variables in (1) how names work (2) how attributes work, and
(3) the macro syntax for constructing them. They're a bit different because
they're parameterized by function-set type. Describe constraints vs.
`ConstraintRefs`. Describe `JuMP.constraint_object`. How to delete constraints.
How to modify constraints by setting attributes and `MOI.modifyconstraint!`.
Describe semidefinite constraints and symmetry handling. Refer to NLP docs for
nonlinear constraints.

```@docs
@constraint
@SDconstraint
```

## Sets

As mentioned in the documentation of the [`@constraint`](@ref) and
[`@SDconstraint`](@ref) macros, the following sets can be used to create
constraints in addition to [any MOI set](http://www.juliaopt.org/MathOptInterface.jl/v0.6.2/apireference.html#Sets-1).

```@docs
SecondOrderCone
RotatedSecondOrderCone
PSDCone
```

## Constraint modifications

`JuMP.set_coefficient(constraint, variable, value)`

```@docs
JuMP.set_coefficient
```

## Duals

```@docs
JuMP.dual
JuMP.shadow_price
```

## Constructing constraints without adding them to the model

For advanced use cases.

```@docs
JuMP.@build_constraint
```
