```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    const MOI = JuMP.MathOptInterface
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in "]
```

# Constraints

This page explains how to write various types of constraints in JuMP. Before
reading further, please make sure you are familiar with JuMP models, and JuMP
[Variables](@ref). If you want to add nonlinear constraints, read
[Nonlinear Modeling](@ref) instead.

JuMP is based on the MathOptInterface API. Because of this, JuMP thinks of a
constraint as a *function* belonging to a *set*. For example, instead of
thinking about a constraint ``a^\top x \le b`` as a *less-than-or-equal-to*
constraint, JuMP thinks about this as the *scalar affine* function ``a^\top x``
belonging to the *less-than* set ``(-\infty, b]``. Thus, instead of a
*less-than-or-equal-to* constraint, we consider this constraint to be a *scalar
affine -in- less than* constraint. More generally, we use the shorthand
*function-in-set* to refer to constraints composed of different types of
functions and sets. In the rest of this page, we will introduce the different
types of functions and sets that JuMP knows about as needed. You can read more
details about this *function-in-set* concept in the MathOptInterface
documentation.

!!! note
    Throughout this page (and these docs), we use `MOI` as a shorthand for the
    `MathOptInterface` module. This can be created by including the following
    lines after `using JuMP` in your code.
    ```julia
    using MathOptInterface
    const MOI = MathOptInterface
    ```

## The `@constraint` macro

Constraints are added to a JuMP model using the [`@constraint`](@ref) macro. It
is similar to the [`@variable`](@ref) macro. Here is an example of how to add
the constraint ``2x \le 1`` to a JuMP model:
```jldoctest con1; setup = :(model = Model(); @variable(model, x))
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
2 x = 1.0

julia> @constraint(model, 1 <= 2x <= 3)
2 x ∈ [1.0, 3.0]
```

Note that JuMP normalizes the constraints given by the user by moving all of the
terms containing variables to the left-hand side, and all of the constant terms
to the right-hand side. Thus, we get:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, 2x + 1 <= 4x + 4)
-2 x <= 3.0
```

## [Duality](@id constraint_duality)

JuMP adopts the notion of duality from MathOptInterface. Roughly speaking, the
dual variable associated with a constraint can be viewed as the decrease in the
value of the objective given an infinitesimal relaxation of the constraint. For
example, in a linear constraint, this relaxation is equivalent to `2x <= 1 + δ`
and `2x >= 1 - δ`. If the constraint is an equality constraint, it depends on
which direction is binding.

Notably, this definition is independent of the objective sense. That is, the
dual associated with a constraint is the same when the objective is a
maximization and when it is a minimization. **This is different to traditional
duality in linear programming.**

The dual value associated with a constraint can be accessed using the
[`JuMP.dual`](@ref) function. For example:

```jldoctest
julia> model = Model()
A JuMP Model

julia> @variable(model, x)
x

julia> @constraint(model, con, x <= 1)
con : x <= 1.0
```

```@meta
DocTestSetup = quote
    using JuMP
    const MOI = JuMP.MathOptInterface
    model = Model(
        with_optimizer(
            MOI.Utilities.MockOptimizer,
            JuMP.JuMPMOIModel{Float64}(),
            eval_objective_value = false,
            eval_variable_constraint_dual = false));
    @variable(model, x);
    @constraint(model, con, x <= 1);
    @objective(model, Max, -2x);
    JuMP.optimize!(model);
    mock = JuMP.caching_optimizer(model).optimizer;
    MOI.set(mock, MOI.ConstraintDual(), JuMP.optimizer_index(con), -2.0)
end
```

```jldoctest con_duality
julia> @objective(model, Min, -2x)
-2 x

julia> JuMP.optimize!(model)

julia> JuMP.dual(con)
-2.0

julia> @objective(model, Max, 2x)
2 x

julia> JuMP.optimize!(model)

julia> JuMP.dual(con)
-2.0
```

To help users who may be less familiar with conic duality, JuMP provides the
[`JuMP.shadow_price`](@ref) function which returns a value that can be
interpreted as the improvement in the objective given a 1-unit change in the
right-hand side of the constraint. [`JuMP.shadow_price`](@ref) can only be used
on linear constraints with a `<=`, `>=`, or `==` comparison operator.

In the example above, `JuMP.dual(con)` returned `-2.0` regardless of the
optimization sense. However, in the second case when the optimization sense is
`Max`, [`JuMP.shadow_price`](@ref) returns:
```jldoctest con_duality
julia> JuMP.shadow_price(con)
2.0
```

To query the dual variables associated a variable bounds, first obtain a
constraint reference using one of [`JuMP.UpperBoundRef`](@ref),
[`JuMP.LowerBoundRef`](@ref), or [`JuMP.FixRef`](@ref), and then call
[`JuMP.dual`](@ref) on the returned constraint reference.

```@meta
DocTestSetup = quote
    using JuMP
    const MOI = JuMP.MathOptInterface
end
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
Therefore `con[1:b]` will work, but `con[a:b]` will not. If JuMP cannot
determine that the indices are one-based integer ranges (e.g., in the case of
`con[a:b]`), JuMP will create a `JuMPArray` instead.

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
tightest container type that can store the constraints. However, because this
happens at compile time, it does not always make the best choice. Just like in
[`@variable`](@ref), we can force the type of container using the `container`
keyword. For syntax and the reason behind this, take a look at the
[variable docs](@ref variable_forcing).

## Vectorized constraints

We can also add constraints to JuMP using vectorized linear algebra. For
example:

```jldoctest con_vector; setup=:(model = Model())
julia> @variable(model, x[i=1:2])
2-element Array{VariableRef,1}:
 x[1]
 x[2]

julia> A = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4

julia> b = [5, 6]
2-element Array{Int64,1}:
 5
 6

julia> @constraint(model, con, A * x .== b)
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}},JuMP.ScalarShape},1}:
 x[1] + 2 x[2] == 5.0
 3 x[1] + 4 x[2] == 6.0
```

!!! note
    Make sure to use the broadcasting syntax `.==` (or `.>=` and `.<=`). If you
    use a standard comparison, an error will be thrown.

Instead of adding an array of `ScalarAffineFunction-in-EqualTo` constraints, we
can instead construct a `VectorAffineFunction-in-Nonnegatives` constraint as
follows:
```jldoctest con_vector
julia> @constraint(model, A * x - b in MOI.Nonnegatives(2))
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] in MathOptInterface.Nonnegatives(2)
```

In addition to the `Nonnegatives` set, MathOptInterface defines a number of
other vector-valued sets such as `Nonpositives`. See the [MathOptInterface
documentation](http://www.juliaopt.org/MathOptInterface.jl/dev/apireference/#Sets-1)
for more information.

Note also that for the first time we have used an explicit *function-in-set*
description of the constraint. Read more about this below in the [Function-Set
pairs](@ref) section of this documentation.

## Constraints on a single variable

In [Variables](@ref), we saw how to modify the variable bounds, as well as add
binary and integer restrictions to the domain of each variable. This can also be
achieved using the [`@constraint`](@ref) macro. For example, `MOI.ZeroOne()`
restricts the domain to ``\{0, 1\}:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.ZeroOne())
x in MathOptInterface.ZeroOne()
```
and `MOI.Integer()` restricts to the domain to the integers ``\mathbb{Z}``:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Integer())
x in MathOptInterface.Integer()
```

JuMP also supports modeling semi-continuous variables, whose domain is ``\{0\} ∪
[l, u]``, using the `MOI.Semicontinuous` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semicontinuous(1.5, 3.5))
x in MathOptInterface.Semicontinuous{Float64}(1.5, 3.5)
```
as well as semi-integer variables, whose domain is ``{0} ∪ {l, l+1, \dots, u}``,
using the `MOI.Semiinteger` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semiinteger(1.0, 3.0))
x in MathOptInterface.Semiinteger{Float64}(1.0, 3.0)
```

## Constraints on a collection of variables

In addition to constraining the domain of a single variable, JuMP supports
placing constraints of a subset of the variables. We already saw an example of
this in the [Quadratic constraints](@ref) section when we constrained a vector
of variables to belong to the second order cone.

In a special ordered set of type I (often denoted SOS-I), at most one variable
can take a non-zero value. We can construct SOS-I constraints using the
`MOI.SOS1` set:
```jldoctest con_sos; setup=:(model = Model())
julia> @variable(model, x[1:3])
3-element Array{VariableRef,1}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in MOI.SOS1([1.0, 2.0, 3.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([1.0, 2.0, 3.0])
```
Note that we have to pass `MOI.SOS1` a *weight* vector. This vector implies an
ordering on the variables. If the decision variables are related and have a
physical ordering (e.g., they correspond to the size of a factory to be built,
and the SOS-I constraint enforces that only one factory can be built), then the
weight vector, although not used directly in the constraint, can help the solver
make better decision in the solution process.

This ordering is more important in a special ordered set of type II (SOS-II), in
which at most two values can be non-zero, and if there are two non-zeros, they
must be consecutive according to the ordering. For example, in the following
constraint, the possible non-zero pairs are (`x[1]` and `x[3]`) and (`x[2]` and
`x[3]`):
```jldoctest con_sos
julia> @constraint(model, x in MOI.SOS2([3.0, 1.0, 2.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])
```

## Quadratic constraints

In addition to affine functions, JuMP also supports constraints with quadratic
terms. (For more general nonlinear functions, read [Nonlinear Modeling](@ref).)
For example:
```jldoctest con_quadratic; setup=:(model=Model())
julia> @variable(model, x[i=1:2])
2-element Array{VariableRef,1}:
 x[1]
 x[2]

julia> @variable(model, t)
t

julia> @constraint(model, x[1]^2 + x[2]^2 <= t^2)
x[1]² + x[2]² - t² <= 0.0
```
Note that this quadratic constraint is equivalent to a second order cone
constraint where `||x[1]^2 + x[2]^2||\_2 ≤ t` and `t ≥ 0`. Instead of writing
out the quadratic expansion, we can pass JuMP the constraint in
*function*-in-*set* form. To do so, we need to define the function and the set.

The function is a vector of variables:
```jldoctest con_quadratic
julia> [t, x[1], x[2]]
3-element Array{VariableRef,1}:
 t
 x[1]
 x[2]
```
Note that the variable `t` comes first, followed by the `x` arguments. The set
is an instance of [`JuMP.SecondOrderCone`](@ref): `JuMP.SecondOrderCone()`.
Thus, we can add the second order cone constraint as follows:
```jldoctest con_quadratic
julia> @constraint(model, [t, x[1], x[2]] in JuMP.SecondOrderCone())
[t, x[1], x[2]] in MathOptInterface.SecondOrderCone(3)
```

JuMP also supports the [`RotatedSecondOrderCone`](@ref) which requires the
addition of a perspective variable `u`. The rotated second order cone
constraints the variables `t`, `u`, and `x` such that: `||x[1]^2 + x[2]^2||\_2 ≤
t × u` and `t, u ≥ 0`. It can be added as follows:
```jldoctest con_quadratic
julia> @variable(model, u)
u

julia> @constraint(model, [t, u, x[1], x[2]] in JuMP.RotatedSecondOrderCone())
[t, u, x[1], x[2]] in MathOptInterface.RotatedSecondOrderCone(4)
```

In addition to the second order cone and rotated second order cone,
MathOptInterface defines a number of other conic sets such as the exponential
and power cones. See the [MathOptInterface documentation](http://www.juliaopt.org/MathOptInterface.jl/dev/apireference/#Sets-1)
for more information.

## Semidefinite constraints

TODO: discuss [`@SDconstraint`] and [`PSDCone`].

## Constraint modifications

A common paradigm, especially in linear programming, is to repeatedly solve a
model with different coefficients. Most often, this involves changing the
"right-hand side" of a linear constraint. This presents a challenge for JuMP
because it leads to ambiguities. For example, what is the right-hand side term
of `@constraint(model, 2x + 1 <= x - 3)`?

To avoid these ambiguities, JuMP includes the ability to *fix* variables to a
value using the [`JuMP.fix`](@ref) function. Fixing a variable sets its lower
and upper bound to the same value. Thus, changes in the right-hand side can be
simulated by adding a dummy variable and fixing it to different values. Here is
an example:

```jldoctest con_fix; setup = :(model = Model(); @variable(model, x))
julia> @variable(model, rhs_term)
rhs_term

julia> @constraint(model, con, 2x <= rhs_term)
con : 2 x - rhs_term <= 0.0

julia> JuMP.fix(rhs_term, 1.0)
```
!!! note
    Even though `rhs_term` is fixed, it is still a decision variable. Thus, `x *
    rhs_term` is bilinear.

It is also possible to modify the scalar coefficients (but notably *not* the
quadratic coefficients) using the [`JuMP.set_coefficient`](@ref) function. Here
is an example:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> JuMP.set_coefficient(con, x, 3)

julia> con
con : 3 x <= 1.0
```

## Constraint deletion

Constraints can be deleted from a model using [`JuMP.delete`](@ref). Just like
variable references, it is possible to check if a constraint reference is valid
using [`JuMP.is_valid`](@ref). Here is an example of deleting a constraint:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> JuMP.is_valid(model, con)
true

julia> JuMP.delete(model, con)

julia> JuMP.is_valid(model, con)
false
```

## Function-Set pairs

DRAFT: Describe how constraints are represented (link to MOI docs). Constraints
are very similar to variables in (1) how names work (2) how attributes work, and
(3) the macro syntax for constructing them. They're a bit different because
they're parameterized by function-set type. Describe constraints vs.
`ConstraintRefs`. Describe `JuMP.constraint_object`. How to delete constraints.
How to modify constraints by setting attributes and `MOI.modifyconstraint!`.
Describe semidefinite constraints and symmetry handling. Refer to NLP docs for
nonlinear constraints.


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
JuMP.fix
JuMP.delete
JuMP.is_valid
JuMP.LowerBoundRef
JuMP.UpperBoundRef
JuMP.FixRef
```

## Constructing constraints without adding them to the model

For advanced use cases.

```@docs
JuMP.@build_constraint
```
