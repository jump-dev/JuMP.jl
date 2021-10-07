```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Constraints](@id jump_constraints)

This page explains how to write various types of constraints in JuMP. Before
reading further, please make sure you are familiar with JuMP models, and JuMP
[Variables](@ref jump_variables). For nonlinear constraints, see
[Nonlinear Modeling](@ref) instead.

JuMP is based on the MathOptInterface (MOI) API. Because of this, JuMP thinks of a
constraint as the restriction that the output of a *function* belongs to a
*set*. For example, instead of representing a constraint ``a^\top x \le b`` as a
*less-than-or-equal-to* constraint, JuMP models this as the *scalar affine*
function ``a^\top x`` belonging to the *less-than* set ``(-\infty, b]``. Thus,
instead of a *less-than-or-equal-to* constraint, we consider this constraint to
be a *scalar affine -in- less than* constraint. More generally, we use the
shorthand *function-in-set* to refer to constraints composed of different types
of functions and sets. In the rest of this page, we will introduce the different
types of functions and sets that JuMP knows about as needed. You can read more
details about this *function-in-set* concept in the MOI documentation.

!!! note
    The examples use `MOI` as an alias for the `MathOptInterface` module. This
    alias is defined by `using JuMP`. You may also define it in your code by
    ```julia
    import MathOptInterface
    const MOI = MathOptInterface
    ```

## The `@constraint` macro

Constraints are added to a JuMP model using the [`@constraint`](@ref) macro.
Here is an example of how to add the constraint ``2x \le 1`` to a JuMP model:
```jldoctest con1; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0
```
Wasn't that easy! Let's unpack what happened, because just like
[`@variable`](@ref) there are a few subtle things going on.
 1. The mathematical constraint ``2x \le 1`` was added to the model.
 2. A Julia variable called `con` was created that is a reference to the
    constraint.
 3. This Julia variable was stored in `model` and can be accessed by
    `model[:con]`.
 4. JuMP set the name attribute (the one that is shown when printing) of the
    constraint to `"con"`.

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

Note that JuMP normalizes the constraints by moving all of the terms containing
variables to the left-hand side, and all of the constant terms to the right-hand
side. Thus, we get:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, 2x + 1 <= 4x + 4)
-2 x <= 3.0
```

## The `@constraints` macro

Like [`@variables`](@ref variables), there is a "plural" version of the
[`@constraint`](@ref) macro:
```jldoctest; setup=:(model=Model(); @variable(model, x))
julia> @constraints(model, begin
           2x <=  1
            x >= -1
       end)

julia> print(model)
Feasibility
Subject to
 x ≥ -1.0
 2 x ≤ 1.0
```

## [Duality](@id constraint_duality)

JuMP adopts the notion of [conic duality from MOI](https://jump.dev/MathOptInterface.jl/v0.9.1/apimanual/#Duals-1).
For linear programs, a feasible dual on a `>=` constraint is nonnegative and a
feasible dual on a `<=` constraint is nonpositive. If the constraint is an
equality constraint, it depends on which direction is binding.

!!! note
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.** If you have a linear program, and you
    want the textbook definition, you probably want to use [`shadow_price`](@ref)
    and [`reduced_cost`](@ref) instead.

The dual value associated with a constraint in the most recent solution can be
accessed using the [`dual`](@ref) function. You can use the [`has_duals`](@ref)
function to check whether the model has a dual solution available to query.
For example:

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, x <= 1)
con : x <= 1.0

julia> has_duals(model)
false
```
```@meta
DocTestSetup = quote
    using JuMP
    model = Model(() -> MOIU.MockOptimizer(
                            MOIU.Model{Float64}(),
                            eval_objective_value = false,
                            eval_variable_constraint_dual = false));
    @variable(model, x);
    @constraint(model, con, x <= 1);
    @objective(model, Max, -2x);
    optimize!(model);
    mock = unsafe_backend(model);
    MOI.set(mock, MOI.TerminationStatus(), MOI.OPTIMAL)
    MOI.set(mock, MOI.DualStatus(), MOI.FEASIBLE_POINT)
    MOI.set(mock, MOI.ConstraintDual(), optimizer_index(con), -2.0)
end
```

```jldoctest con_duality
julia> @objective(model, Min, -2x)
-2 x

julia> optimize!(model)

julia> has_duals(model)
true

julia> dual(con)
-2.0

julia> @objective(model, Max, 2x)
2 x

julia> optimize!(model)

julia> dual(con)
-2.0
```

To help users who may be less familiar with conic duality, JuMP provides the
[`shadow_price`](@ref) function which returns a value that can be
interpreted as the improvement in the objective in response to an infinitesimal
relaxation (on the scale of one unit) in the right-hand side of the constraint.
[`shadow_price`](@ref) can be used only on linear constraints with a `<=`,
`>=`, or `==` comparison operator.

In the example above, `dual(con)` returned `-2.0` regardless of the
optimization sense. However, in the second case when the optimization sense is
`Max`, [`shadow_price`](@ref) returns:
```jldoctest con_duality
julia> shadow_price(con)
2.0
```

To query the dual variables associated with a variable bound, first obtain a
constraint reference using one of [`UpperBoundRef`](@ref),
[`LowerBoundRef`](@ref), or [`FixRef`](@ref), and then call [`dual`](@ref) on
the returned constraint reference. The [`reduced_cost`](@ref) function may
simplify this process as it returns the shadow price of an active bound of
a variable (or zero, if no active bound exists).

```@meta
DocTestSetup = quote
    using JuMP
end
```

## Constraint names

The name, i.e. the value of the `MOI.ConstraintName` attribute, of a constraint
can be obtained by [`name(::JuMP.ConstraintRef)`](@ref) and set by
[`set_name(::JuMP.ConstraintRef, ::String)`](@ref).

```jldoctest
julia> model = Model(); @variable(model, x);

julia> @constraint(model, con, x <= 1)
con : x <= 1.0

julia> name(con)
"con"

julia> set_name(con, "my_con_name")

julia> con
my_con_name : x <= 1.0
```

Specify a constraint name in the macro via `base_name`:
```jldoctest constraint_name_vector
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

Note that names apply to each element of the container, not to the container of
constraints:
```jldoctest constraint_name_vector
julia> name(con[1])
"my_con[1]"

julia> set_name(con[1], "c")

julia> con
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

### Retrieve a constraint by name

Retrieve a constraint from a model using [`constraint_by_name`](@ref):
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "c")
c : x ≤ 1.0
```
If the name is not present, `nothing` will be returned:
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "bad_name")
```

You can only look up invididual constraints using [`constraint_by_name`](@ref).
Something like this will not work:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0

julia> constraint_by_name(model, "my_con")
```

To look up a collection of constraints, do not use [`constraint_by_name`](@ref).
Instead, register them using the `model[:key] = value` syntax:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> model[:con] = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0

julia> model[:con]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1.0
 my_con[2] : x ≤ 2.0
```

## Start Values

Provide a starting value (also called warmstart) for a constraint's dual using
[`set_dual_start_value`](@ref).

The start value of a constraint's dual can be queried using [`dual_start_value`](@ref).
If no start value has been set, [`dual_start_value`](@ref) will return `nothing`.

```jldoctest constraint_dual_start; setup=:(model=Model())
julia> @variable(model, x)
x

julia> @constraint(model, con, x >= 10)
con : x ≥ 10.0

julia> dual_start_value(con)

julia> set_dual_start_value(con, 2)

julia> dual_start_value(con)
2.0
```

A vector constraint will require a vector warmstart:

```jldoctest constraint_dual_start_vector; setup=:(model=Model())
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, con, x in SecondOrderCone())
con : [x[1], x[2], x[3]] in MathOptInterface.SecondOrderCone(3)

julia> dual_start_value(con)

julia> set_dual_start_value(con, [1.0, 2.0, 3.0])

julia> dual_start_value(con)
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```

To take the dual solution from the last solve and use it as the starting point
for a new solve, use:

```julia
for (F, S) in list_of_constraint_types(model)
    for con in all_constraints(model, F, S)
        set_dual_start_value(con, dual(con))
    end
end
```

!!! note
    Some constraints might not have well defined duals, hence one might need to
    filter `(F, S)` pairs.

## Constraint containers

So far, we've added constraints one-by-one. However, just like
[Variable containers](@ref), JuMP provides a mechanism for building groups of
constraints compactly. References to these groups of constraints are returned in
*containers*. Three types of constraint containers are supported: `Array`s,
`DenseAxisArray`s, and `SparseAxisArray`s. We explain each of these in the
following.

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### [Arrays](@id constraint_arrays)

One way of adding a group of constraints compactly is the following:
```jldoctest constraint_arrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:3], i * x <= i + 1)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1] : x ≤ 2.0
 con[2] : 2 x ≤ 3.0
 con[3] : 3 x ≤ 4.0
```
JuMP returns references to the three constraints in an `Array` that is bound to
the Julia variable `con`. This array can be accessed and sliced as you would
with any Julia array:
```jldoctest constraint_arrays
julia> con[1]
con[1] : x <= 2.0

julia> con[2:3]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[2] : 2 x ≤ 3.0
 con[3] : 3 x ≤ 4.0
```

Anonymous containers can also be constructed by dropping the name (e.g. `con`)
before the square brackets:
```jldoctest constraint_arrays
julia> @constraint(model, [i = 1:2], i * x <= i + 1)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x ≤ 2.0
 2 x ≤ 3.0
```

Just like [`@variable`](@ref), JuMP will form an `Array` of constraints when it
can determine at parse time that the indices are one-based integer ranges.
Therefore `con[1:b]` will create an `Array`, but `con[a:b]` will not. A special
case is `con[Base.OneTo(n)]` which will produce an `Array`. If JuMP cannot
determine that the indices are one-based integer ranges (e.g., in the case of
`con[a:b]`), JuMP will create a `DenseAxisArray` instead.

### DenseAxisArrays

The syntax for constructing a [`DenseAxisArray`](@ref Containers.DenseAxisArray)
of constraints is very similar to the
[syntax for constructing](@ref variable_jump_arrays) a `DenseAxisArray` of
variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 2:3], i * x <= j + 1)
2-dimensional DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 2:3
And data, a 2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1,2] : x ≤ 3.0    con[1,3] : x ≤ 4.0
 con[2,2] : 2 x ≤ 3.0  con[2,3] : 2 x ≤ 4.0
```

### SparseAxisArrays

The syntax for constructing a
[`SparseAxisArray`](@ref Containers.SparseAxisArray) of constraints is very
similar to the [syntax for constructing](@ref variable_sparseaxisarrays) a
`SparseAxisArray` of variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 2, Tuple{Int64, Int64}} with 2 entries:
  [1, 2]  =  con[1,2] : x ≤ 3.0
  [2, 1]  =  con[2,1] : 2 x ≤ 2.0
```

### Forcing the container type

When creating a container of constraints, JuMP will attempt to choose the
tightest container type that can store the constraints. However, because this
happens at parse time, it does not always make the best choice. Just like in
[`@variable`](@ref), we can force the type of container using the `container`
keyword. For syntax and the reason behind this, take a look at the
[variable docs](@ref variable_forcing).

## Vectorized constraints

We can also add constraints to JuMP using vectorized linear algebra. For
example:

```jldoctest con_vector; setup=:(model = Model())
julia> @variable(model, x[i=1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> A = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> b = [5, 6]
2-element Vector{Int64}:
 5
 6

julia> @constraint(model, con, A * x .== b)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 x[1] + 2 x[2] = 5.0
 3 x[1] + 4 x[2] = 6.0
```

!!! note
    Make sure to use [Julia's dot syntax](https://docs.julialang.org/en/v1/manual/functions/index.html#man-vectorized-1)
    in front of the comparison operators (e.g. `.==`, `.>=`, and `.<=`). If you
    use a comparison without the dot, an error will be thrown.

Instead of adding an array of `ScalarAffineFunction-in-EqualTo` constraints, we
can instead construct a `VectorAffineFunction-in-Nonnegatives` constraint as
follows:
```jldoctest con_vector
julia> @constraint(model, A * x - b in MOI.Nonnegatives(2))
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] in MathOptInterface.Nonnegatives(2)
```

In addition to the `Nonnegatives` set, MOI defines a number of
other vector-valued sets such as `Nonpositives`. See the
[MOI documentation](https://jump.dev/MathOptInterface.jl/v0.9.1/apireference/#Sets-1)
for more information.

Note also that for the first time we have used an explicit *function-in-set*
description of the constraint. Read more about this representation for
constraints in the
[MOI documentation](https://jump.dev/MathOptInterface.jl/v0.9.1/apimanual/#Constraints-by-function-set-pairs-1).

## Constraints on a single variable

In [Variables](@ref jump_variables), we saw how to modify the variable bounds,
as well as add binary and integer restrictions to the domain of each variable.
This can also be achieved using the [`@constraint`](@ref) macro. For example,
`MOI.ZeroOne()` restricts the domain to ``\{0, 1\}``:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.ZeroOne())
x binary
```
and `MOI.Integer()` restricts to the domain to the integers ``\mathbb{Z}``:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Integer())
x integer
```

JuMP also supports modeling semi-continuous variables, whose domain is ``\{0\} ∪
[l, u]``, using the `MOI.Semicontinuous` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semicontinuous(1.5, 3.5))
x in MathOptInterface.Semicontinuous{Float64}(1.5, 3.5)
```
as well as semi-integer variables, whose domain is ``\{0\} ∪ \{l, l+1, \dots, u\}``,
using the `MOI.Semiinteger` set:
```jldoctest; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, x in MOI.Semiinteger(1.0, 3.0))
x in MathOptInterface.Semiinteger{Float64}(1.0, 3.0)
```

## [Quadratic constraints](@id quad_constraints)

In addition to affine functions, JuMP also supports constraints with quadratic
terms. (For more general nonlinear functions, see [Nonlinear Modeling](@ref).)
For example:
```jldoctest con_quadratic; setup=:(model=Model())
julia> @variable(model, x[i=1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @variable(model, t >= 0)
t

julia> @constraint(model, x[1]^2 + x[2]^2 <= t^2)
x[1]² + x[2]² - t² <= 0.0
```
Note that this quadratic constraint (including the lower bound on `t`) is
equivalent to a second order cone constraint where `||x[1]^2 + x[2]^2||\_2 ≤ t`
and `t ≥ 0`. Instead of writing out the quadratic expansion, we can pass JuMP
the constraint in *function*-in-*set* form. To do so, we need to define the
function and the set.

The function is a vector of variables:
```jldoctest con_quadratic
julia> [t, x[1], x[2]]
3-element Vector{VariableRef}:
 t
 x[1]
 x[2]
```
Note that the variable `t` comes first, followed by the `x` arguments. The set
is an instance of [`SecondOrderCone`](@ref): `SecondOrderCone()`.
Thus, we can add the second order cone constraint as follows:
```jldoctest con_quadratic
julia> @constraint(model, [t, x[1], x[2]] in SecondOrderCone())
[t, x[1], x[2]] in MathOptInterface.SecondOrderCone(3)
```

JuMP also supports the [`RotatedSecondOrderCone`](@ref) which requires the
addition of a perspective variable `u`. The rotated second order cone
constraints the variables `t`, `u`, and `x` such that: `||x[1]^2 + x[2]^2||\_2 ≤
t × u` and `t, u ≥ 0`. It can be added as follows:
```jldoctest con_quadratic
julia> @variable(model, u)
u

julia> @constraint(model, [t, u, x[1], x[2]] in RotatedSecondOrderCone())
[t, u, x[1], x[2]] in MathOptInterface.RotatedSecondOrderCone(4)
```

In addition to the second order cone and rotated second order cone,
MOI defines a number of other conic sets such as the exponential
and power cones. See the [MathOptInterface documentation](https://jump.dev/MathOptInterface.jl/v0.9.1/apireference/#Sets-1)
for more information.

## Constraints on a collection of variables

In addition to constraining the domain of a single variable, JuMP supports
placing constraints of a subset of the variables. We already saw an example of
this in the [Quadratic constraints](@ref quad_constraints) section when we
constrained a vector of variables to belong to the second order cone.

In a special ordered set of type I (often denoted SOS-I), at most one variable
can take a non-zero value. We can construct SOS-I constraints using the
`MOI.SOS1` set:
```jldoctest con_sos; setup=:(model = Model())
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
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
make a better decision in the solution process.

This ordering is more important in a special ordered set of type II (SOS-II), in
which at most two values can be non-zero, and if there are two non-zeros, they
must be consecutive according to the ordering. For example, in the following
constraint, the possible non-zero pairs are (`x[1]` and `x[3]`) and (`x[2]` and
`x[3]`):
```jldoctest con_sos
julia> @constraint(model, x in MOI.SOS2([3.0, 1.0, 2.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])
```

## Indicator constraints

JuMP provides a special syntax for creating indicator constraints, that is,
enforce a constraint to hold depending on the value of a binary variable.
In order to constrain the constraint `x + y <= 1` to hold when a binary
variable `a` is one, use the following syntax:
```jldoctest indicator; setup=:(model = Model())
julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @variable(model, a, Bin)
a

julia> @constraint(model, a => {x + y <= 1})
a => {x + y ≤ 1.0}
```
If instead the constraint should hold when `a` is zero, simply add a `!` or `¬`
before the binary variable.
```jldoctest indicator
julia> @constraint(model, !a => {x + y <= 1})
!a => {x + y ≤ 1.0}
```

## Semidefinite constraints

Use [`PSDCone`](ref) to constrain a matrix to be symmetric positive
semidefinite (PSD). For example,
```jldoctest con_psd; setup=:(model = Model())
julia> @variable(model, X[1:2, 1:2])
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[2,1]  X[2,2]

julia> @constraint(model, X >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[2,1]  X[2,2]] ∈ PSDCone()
```

The inequality `X >= Y` between two square matrices `X` and `Y` is understood as
constraining `X - Y` to be symmetric positive semidefinite.
```jldoctest con_psd
julia> Y = [1 2; 2 1]
2×2 Matrix{Int64}:
 1  2
 2  1

julia> @constraint(model, X >= Y, PSDCone())
[X[1,1] - 1  X[1,2] - 2;
 X[2,1] - 2  X[2,2] - 1] ∈ PSDCone()
```

!!! tip
    `@constraint(model, X >= Y, Set())` is short-hand for
    `@constraint(model, X - Y in Set())`. Therefore, the following calls are
    equivalent:
     * `@constraint(model, X >= Y, PSDCone())`
     * `@constraint(model, Y <= X, PSDCone())`
     * `@constraint(model, X - Y in PSDCone())`
    This also works for any vector-valued cone, so if `x` and `y` are vectors of
    length 2, you can write `@constraint(model, x >= y, MOI.Nonnegatives(2))`
    instead of `@constraint(model, x - y in MOI.Nonnegatives(2))`.

!!! warning
    Non-zero constants are not supported in this syntax:
    ```jldoctest con_psd
    julia> @constraint(model, X >= 1, PSDCone())
    ERROR: Operation `sub_mul!` between `Matrix{VariableRef}` and `Int64` is not allowed. You should use broadcast.
    Stacktrace:
    [...]
    ```
    Use instead:
    ```jldoctest con_psd
    julia> @constraint(model, X .- 1 >= 0, PSDCone())
    [X[1,1] - 1  X[1,2] - 1;
     X[2,1] - 1  X[2,2] - 1] ∈ PSDCone()
    ```

### Symmetry

Solvers supporting PSD constraints usually expect to be given a matrix that
is *symbolically* symmetric, that is, for which the expression in corresponding
off-diagonal entries are the same. In our example, the expressions of entries
`(1, 2)` and `(2, 1)` are respectively `X[1,2] - 2` and `X[2,1] - 2` which are
different.

To bridge the gap between the constraint modeled and what the solver
expects, solvers may add an equality constraint `X[1,2] - 2 == X[2,1] - 2` to
force symmetry. Use `LinearAlgebra.Symmetric` to explicitly tell the solver that
the matrix is symmetric:
```jldoctest con_psd
julia> import LinearAlgebra

julia> Z = [X[1, 1] X[1, 2]; X[1, 2] X[2, 2]]
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[1,2]  X[2,2]

julia> @constraint(model, LinearAlgebra.Symmetric(Z) >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[1,2]  X[2,2]] ∈ PSDCone()
```

Note that the lower triangular entries are silently ignored even if they are
different so use it with caution:
```jldoctest con_psd
julia> @constraint(model, LinearAlgebra.Symmetric(X) >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[1,2]  X[2,2]] ∈ PSDCone()
```
(Note the `(2, 1)` element of the constraint is `X[1,2]`, not `X[2,1]`.)

## Modify a constraint

### Modifying a constant term (Option 1)

Use [`set_normalized_rhs`](@ref) to modify the right-hand side (constant)
term of a constraint. Use [`normalized_rhs`](@ref) to query the right-hand
side term.

```jldoctest con_fix; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> set_normalized_rhs(con, 3)

julia> con
con : 2 x <= 3.0

julia> normalized_rhs(con)
3.0
```

!!! note
    JuMP normalizes constraints into a standard form by moving all constant
    terms onto the right-hand side of the constraint.
    ```julia
    @constraint(model, 2x - 1 <= 2)
    ```
    will be normalized to
    ```julia
    @constraint(model, 2x <= 3)
    ```
    [`set_normalized_rhs`](@ref) sets the right-hand side term of the
    normalized constraint.

### Modifying a constant term (Option 2)

If constraints are complicated, e.g., they are composed of a number of
components, each of which has a constant term, then it may be difficult to
calculate what the right-hand side term should be in the standard form.

For this situation, JuMP includes the ability to *fix* variables to a
value using the [`fix`](@ref) function. Fixing a variable sets its lower
and upper bound to the same value. Thus, changes in a constant term can be
simulated by adding a dummy variable and fixing it to different values. Here is
an example:

```jldoctest con_fix; setup = :(model = Model(); @variable(model, x))
julia> @variable(model, const_term)
const_term

julia> @constraint(model, con, 2x <= const_term + 1)
con : 2 x - const_term <= 1.0

julia> fix(const_term, 1.0)
```
The constraint `con` is now equivalent to `2x <= 2`.

!!! note
    Even though `const_term` is fixed, it is still a decision variable. Thus,
    `const_term * x` is bilinear. Fixed variables are not replaced with
    constants when communicating the problem to a solver.

Another option is to use [`add_to_function_constant`](@ref). The constant given
is added to the function of a `func`-in-`set` constraint. In the following
example, adding `2` to the function has the effect of removing `2` to the
right-hand side:
```jldoctest con_add; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> add_to_function_constant(con, 2)

julia> con
con : 2 x <= -1.0

julia> normalized_rhs(con)
-1.0
```

In the case of interval constraints, the constant is removed in each bounds.
```jldoctest con_add_interval; setup = :(model = Model(); @variable(model, x))
julia> @constraint(model, con, 0 <= 2x + 1 <= 2)
con : 2 x ∈ [-1.0, 1.0]

julia> add_to_function_constant(con, 3)

julia> con
con : 2 x ∈ [-4.0, -2.0]
```

### Modifying a variable coefficient

To modify the coefficients for a linear term in a constraint (but
notably not yet the coefficients on a quadratic term), use
[`set_normalized_coefficient`](@ref). To query
the current coefficient, use [`normalized_coefficient`](@ref).
```jldoctest; setup = :(model = Model(); @variable(model, x[1:2]))
julia> @constraint(model, con, 2x[1] + x[2] <= 1)
con : 2 x[1] + x[2] ≤ 1.0

julia> set_normalized_coefficient(con, x[2], 0)

julia> con
con : 2 x[1] ≤ 1.0

julia> normalized_coefficient(con, x[2])
0.0
```

!!! note
    JuMP normalizes constraints into a standard form by moving all terms
    involving variables onto the left-hand side of the constraint.
    ```julia
    @constraint(model, 2x <= 1 - x)
    ```
    will be normalized to
    ```julia
    @constraint(model, 3x <= 1)
    ```
    [`set_normalized_coefficient`](@ref) sets the coefficient of the
    normalized constraint.

## Delete a constraint

Use [`delete`](@ref) to delete a constraint from a model. Use [`is_valid`](@ref)
to check if a constraint belongs to a model and has not been deleted.

```jldoctest constraints_delete; setup = :(model=Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> is_valid(model, con)
true

julia> delete(model, con)

julia> is_valid(model, con)
false
```

Deleting a constraint does not unregister the symbolic reference from the model.
Therefore, creating a new constraint of the same name will throw an error:
```jldoctest constraints_delete
julia> @constraint(model, con, 2x <= 1)
ERROR: An object of name con is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :con)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:con]`.
[...]
```

After calling [`delete`](@ref), call [`unregister`](@ref) to remove the symbolic
reference:
```jldoctest constraints_delete
julia> unregister(model, :con)

julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0
```

!!! info
    [`delete`](@ref) does not automatically [`unregister`](@ref) because we do
    not distinguish between names that are automatically registered by JuMP
    macros, and names that are manually registered by the user by setting values
    in [`object_dictionary`](@ref). In addition, deleting a constraint and then
    adding a new constraint of the same name is an easy way to introduce bugs
    into your code.

## Accessing constraints from a model

You can query the types of constraints currently present in the model by calling
[`list_of_constraint_types`](@ref). Then, given a function and set type, use
[`num_constraints`](@ref) to access the number of constraints of this type and
[`all_constraints`](@ref) to access a list of their references. Then use
[`constraint_object`](@ref) to get an instance of an
[`AbstractConstraint`](@ref) object, either [`ScalarConstraint`](@ref) or
[`VectorConstraint`](@ref), that stores the constraint data.

```jldoctest
julia> model = Model();

julia> @variable(model, x[i=1:2] >= i, Int);

julia> @constraint(model, x[1] + x[2] <= 1);

julia> list_of_constraint_types(model)
3-element Vector{Tuple{Type, Type}}:
 (AffExpr, MathOptInterface.LessThan{Float64})
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (VariableRef, MathOptInterface.Integer)

julia> num_constraints(model, VariableRef, MOI.Integer)
2

julia> all_constraints(model, VariableRef, MOI.Integer)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex, MathOptInterface.Integer}, ScalarShape}}:
 x[1] integer
 x[2] integer

julia> num_constraints(model, VariableRef, MOI.GreaterThan{Float64})
2

julia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex, MathOptInterface.GreaterThan{Float64}}, ScalarShape}}:
 x[1] ≥ 1.0
 x[2] ≥ 2.0

julia> num_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
1

julia> less_than_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
1-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x[1] + x[2] ≤ 1.0

julia> con = constraint_object(less_than_constraints[1])
ScalarConstraint{AffExpr, MathOptInterface.LessThan{Float64}}(x[1] + x[2], MathOptInterface.LessThan{Float64}(1.0))

julia> con.func
x[1] + x[2]

julia> con.set
MathOptInterface.LessThan{Float64}(1.0)
```

## Complementarity constraints

A mixed complementarity constraint `F(x) ⟂ x` consists of finding `x` in the
interval `[lb, ub]`, such that the following holds:

- `F(x) == 0` if `lb < x < ub`
- `F(x) >= 0` if `lb == x`
- `F(x) <= 0` if `x == ub`

For more information, see the [`MOI.Complements` documentation](https://jump.dev/MathOptInterface.jl/v0.9/apireference/#MathOptInterface.Complements).

JuMP supports mixed complementarity constraints via `complements(F(x), x)` or
`F(x) ⟂ x` in the [`@constraint`](@ref) macro. The interval set `[lb, ub]` is
obtained from the variable bounds on `x`.

For example, to define the problem `2x - 1 ⟂ x` with `x ∈ [0, ∞)`, do:
```jldoctest complementarity; setup=:(model=Model())
julia> @variable(model, x >= 0)
x

julia> @constraint(model, 2x - 1 ⟂ x)
[2 x - 1, x] ∈ MathOptInterface.Complements(2)
```
This problem has a unique solution at `x = 0.5`.

The perp operator `⟂` can be entered in most editors (and the Julia REPL) by
typing `\perp<tab>`.

An alternative approach that does not require the `⟂` symbol uses the
`complements` function as follows:
```jldoctest complementarity
julia> @constraint(model, complements(2x - 1, x))
[2 x - 1, x] ∈ MathOptInterface.Complements(2)
```

In both cases, the mapping `F(x)` is supplied as the first argument, and the
matching variable `x` is supplied as the second.

Vector-valued complementarity constraints are also supported:
```jldoctest complementarity
julia> @variable(model, -2 <= y[1:2] <= 2)
2-element Vector{VariableRef}:
 y[1]
 y[2]

julia> M = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> q = [5, 6]
2-element Vector{Int64}:
 5
 6

julia> @constraint(model, M * y + q ⟂ y)
[y[1] + 2 y[2] + 5, 3 y[1] + 4 y[2] + 6, y[1], y[2]] ∈ MathOptInterface.Complements(4)
```

## Special Ordered Sets (SOS1 and SOS2)

### Type 1

In a Special Ordered Set of Type 1 (often denoted SOS-I or SOS1), at most one
element can take a non-zero value.

Construct SOS-I constraints using the [`SOS1`](@ref) set:
```jldoctest con_sos; setup=:(model = Model())
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SOS1())
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([1.0, 2.0, 3.0])
```

Although not required for feasibility, solvers can benefit from an ordering of
the variables (e.g., the variables represent different factories to build, at
most one factory can be built, and the factories can be ordered according to
cost). To induce an ordering, weights can be provided; as such, they should be
unique values. The kth element in the ordered set corresponds to the kth weight
in weights when the weights are sorted.

For example, in the constraint:
```jldoctest con_sos
julia> @constraint(model, x in SOS1([3.1, 1.2, 2.3]))
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([3.1, 1.2, 2.3])
```
the variables `x` have precedence `x[2]`, `x[3]`, `x[1]`.

### Type 2

In a Special Ordered Set of Type 2 (SOS-II), at most two elements can be
non-zero, and if there are two non-zeros, they must be consecutive according to
the ordering induced by a weight vector.

Construct SOS-II constraints using the [`SOS2`](@ref) set.
```jldoctest con_sos
julia> @constraint(model, x in SOS2([3.0, 1.0, 2.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])
```
In the following constraint, the possible non-zero pairs are (`x[1]` and `x[3]`)
and (`x[2]` and `x[3]`):

If the weight vector is omitted, JuMP induces an ordering from `1:length(x)`:
```jldoctest con_sos
julia> @constraint(model, x in SOS2())
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([1.0, 2.0, 3.0])
```
