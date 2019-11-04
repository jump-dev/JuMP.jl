```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Constraints

This page explains how to write various types of constraints in JuMP. Before
reading further, please make sure you are familiar with JuMP models, and JuMP
[Variables](@ref). For nonlinear constraints, see [Nonlinear Modeling](@ref)
instead.

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

JuMP adopts the notion of [conic duality from MOI](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apimanual/#Duals-1).
For linear programs, a feasible dual on a `>=` constraint is nonnegative and a
feasible dual on a `<=` constraint is nonpositive. If the constraint is an
equality constraint, it depends on which direction is binding.

!!! note
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.**

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
    mock = backend(model).optimizer.model;
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
the returned constraint reference. Note that in linear programming,
the duals on variable bounds are also called the reduced costs
(although the sign might differ from the one you expect).

```@meta
DocTestSetup = quote
    using JuMP
end
```

## Constraint names

The name, i.e. the value of the `MOI.ConstraintName` attribute, of a constraint
can be obtained by [`name(::JuMP.ConstraintRef)`](@ref) and set by
[`set_name(::JuMP.ConstraintRef, ::String)`](@ref).
```@docs
name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex})
set_name(::JuMP.ConstraintRef{Model, <:JuMP.MOI.ConstraintIndex}, ::String)
```

The constraint can also be retrieved from its name using
[`constraint_by_name`](@ref).
```@docs
constraint_by_name
```

## Constraint containers

So far, we've added constraints one-by-one. However, just like
[Variable containers](@ref), JuMP provides a mechanism for building groups of
constraints compactly. References to these groups of constraints are returned in
*containers*. Three types of constraint containers are supported: `Array`s,
`DenseAxisArray`s, and `SparseAxisArray`s. We explain each of these in the
following.

### [Arrays](@id constraint_arrays)

One way of adding a group of constraints compactly is the following:
```jldoctest constraint_arrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:3], i * x <= i + 1)
3-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 con[1] : x <= 2.0
 con[2] : 2 x <= 3.0
 con[3] : 3 x <= 4.0
```
JuMP returns references to the three constraints in an `Array` that is bound to
the Julia variable `con`. This array can be accessed and sliced as you would
with any Julia array:
```jldoctest constraint_arrays
julia> con[1]
con[1] : x <= 2.0

julia> con[2:3]
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 con[2] : 2 x <= 3.0
 con[3] : 3 x <= 4.0
```

Anonymous containers can also be constructed by dropping the name (e.g. `con`)
before the square brackets:
```jldoctest constraint_arrays
julia> @constraint(model, [i = 1:2], i * x <= i + 1)
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 x <= 2.0
 2 x <= 3.0
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
2-dimensional DenseAxisArray{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 2:3
And data, a 2×2 Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},2}:
 con[1,2] : x <= 3.0    con[1,3] : x <= 4.0
 con[2,2] : 2 x <= 3.0  con[2,3] : 2 x <= 4.0
```

### SparseAxisArrays

The syntax for constructing a
[`SparseAxisArray`](@ref Containers.SparseAxisArray) of constraints is very
similar to the [syntax for constructing](@ref variable_sparseaxisarrays) a
`SparseAxisArray` of variables.

```jldoctest constraint_jumparrays; setup=:(model=Model(); @variable(model, x))
julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},2,Tuple{Int64,Int64}} with 2 entries:
  [1, 2]  =  con[1,2] : x <= 3.0
  [2, 1]  =  con[2,1] : 2 x <= 2.0
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
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}},ScalarShape},1}:
 x[1] + 2 x[2] == 5.0
 3 x[1] + 4 x[2] == 6.0
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
[MOI documentation](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apireference/#Sets-1)
for more information.

Note also that for the first time we have used an explicit *function-in-set*
description of the constraint. Read more about this representation for
constraints in the
[MOI documentation](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apimanual/#Constraints-by-function-set-pairs-1).

## Constraints on a single variable

In [Variables](@ref), we saw how to modify the variable bounds, as well as add
binary and integer restrictions to the domain of each variable. This can also be
achieved using the [`@constraint`](@ref) macro. For example, `MOI.ZeroOne()`
restricts the domain to ``\{0, 1\}``:
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

## Quadratic constraints

In addition to affine functions, JuMP also supports constraints with quadratic
terms. (For more general nonlinear functions, see [Nonlinear Modeling](@ref).)
For example:
```jldoctest con_quadratic; setup=:(model=Model())
julia> @variable(model, x[i=1:2])
2-element Array{VariableRef,1}:
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
3-element Array{VariableRef,1}:
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
and power cones. See the [MathOptInterface documentation](http://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apireference/#Sets-1)
for more information.

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

## Semidefinite constraints

JuMP provides a special syntax for constraining a matrix to be symmetric
positive semidefinite (PSD) with the [`@SDconstraint`](@ref) macro.
In the context of this macro, the inequality `A >= B` between two square
matrices `A` and `B` is understood as constraining `A - B` to be symmetric
positive semidefinite.
```jldoctest con_psd; setup=:(model = Model())
julia> @variable(model, x)
x

julia> @SDconstraint(model, [x 2x; 3x 4x] >= ones(2, 2))
[x - 1    2 x - 1;
 3 x - 1  4 x - 1] ∈ PSDCone()
```

Solvers supporting such constraints usually expect to be given a matrix that
is *symbolically* symmetric, that is, for which the expression in corresponding
off-diagonal entries are the same. In our example, the expressions of entries
`(1, 2)` and `(2, 1)` are respectively `2x - 1` and `3x - 1` which are
different. To bridge the gap between the constraint modeled and what the solver
expects, JuMP creates an equality constraint `3x - 1 == 2x - 1` and constrains
the symmetric matrix `[x - 1, 2 x - 1, 2 x - 1, 4 x - 1]` to be positive
semidefinite.

!!! note
    If the matrix provided is already symbolically symmetric, the equality
    constrains are equivalent to `0 = 0` and are not added. In practice, if
    all coefficients are smaller than `1e-10`, the constraint is ignored, if
    all coefficients are smaller than `1e-8` but some are larger than `1e-10`,
    it is ignored but a warning is displayed, otherwise if at least one
    coefficient is larger than `1e-8`, the constraint is added.

If the matrix is known to be symmetric, the PSD constraint can be added as
follows:
```jldoctest con_psd
julia> using LinearAlgebra

julia> @constraint(model, Symmetric([x 2x; 2x 4x] - ones(2, 2)) in PSDCone())
[x - 1    2 x - 1;
 2 x - 1  4 x - 1] ∈ PSDCone()
```

Note that the lower triangular entries are silently ignored even if they are
different so use it with caution:
```jldoctest con_psd
julia> cref = @constraint(model, Symmetric([x 2x; 3x 4x]) in PSDCone())
[x    2 x;
 2 x  4 x] ∈ PSDCone()

julia> jump_function(constraint_object(cref))
3-element Array{GenericAffExpr{Float64,VariableRef},1}:
 x
 2 x
 4 x

julia> moi_set(constraint_object(cref))
MathOptInterface.PositiveSemidefiniteConeTriangle(2)
```

## Constraint modifications

A common paradigm, especially in linear programming, is to repeatedly solve a
model with different coefficients.

### Modifying a constant term

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
    JuMP normalizes constraints into a standard form by moving all constant terms
    onto the right-hand side of the constraint.
    ```julia
    @constraint(model, 2x - 1 <= 2)
    ```
    will be normalized to
    ```julia
    @constraint(model, 2x <= 3)
    ```
    [`set_normalized_rhs`](@ref) sets the right-hand side term of the
    normalized constraint.

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

## Constraint deletion

Constraints can be deleted from a model using [`delete`](@ref). Just like
variable references, it is possible to check if a constraint reference is valid
using [`is_valid`](@ref). Here is an example of deleting a constraint:
```jldoctest; setup = :(model=Model(); @variable(model, x))
julia> @constraint(model, con, 2x <= 1)
con : 2 x <= 1.0

julia> is_valid(model, con)
true

julia> delete(model, con)

julia> is_valid(model, con)
false
```

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
3-element Array{Tuple{DataType,DataType},1}:
 (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64})
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (VariableRef, MathOptInterface.Integer)

julia> num_constraints(model, VariableRef, MOI.Integer)
2

julia> all_constraints(model, VariableRef, MOI.Integer)
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.Integer},ScalarShape},1}:
 x[1] integer
 x[2] integer

julia> num_constraints(model, VariableRef, MOI.GreaterThan{Float64})
2

julia> all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
2-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}},ScalarShape},1}:
 x[1] ≥ 1.0
 x[2] ≥ 2.0

julia> num_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
1

julia> less_than_constraints = all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64})
1-element Array{ConstraintRef{Model,MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}},ScalarShape},1}:
 x[1] + x[2] ≤ 1.0

julia> con = constraint_object(less_than_constraints[1])
ScalarConstraint{GenericAffExpr{Float64,VariableRef},MathOptInterface.LessThan{Float64}}(x[1] + x[2], MathOptInterface.LessThan{Float64}(1.0))

julia> con.func
x[1] + x[2]

julia> con.set
MathOptInterface.LessThan{Float64}(1.0)
```


## Reference

```@docs
@constraint
@SDconstraint
SecondOrderCone
RotatedSecondOrderCone
PSDCone
shadow_price
normalized_coefficient
set_normalized_coefficient
normalized_rhs
set_normalized_rhs
add_to_function_constant
is_valid
JuMP.delete
LowerBoundRef
UpperBoundRef
FixRef
ConstraintRef
list_of_constraint_types
all_constraints
num_constraints
constraint_object
AbstractConstraint
ScalarConstraint
VectorConstraint
index(::ConstraintRef)
optimizer_index(::ConstraintRef{Model})
```

## Constructing constraints without adding them to the model

For advanced use cases.

```@docs
@build_constraint
```
