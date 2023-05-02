```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
    import HiGHS
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Constraints](@id jump_constraints)

JuMP is based on the [MathOptInterface (MOI) API](@ref moi_documentation).
Because of this, JuMP uses the following standard form to represent problems:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x)
    \\
    & \;\;\text{s.t.} & f_i(x) & \in \mathcal{S}_i & i = 1 \ldots m
\end{align}
```
Each constraint, ``f_i(x) \in \mathcal{S}_i``, is composed of a function and a
set. For example, instead of calling ``a^\top x \le b`` a
*less-than-or-equal-to* constraint, we say that it is a
*scalar-affine-in-less-than* constraint, where the function ``a^\top x`` belongs
to the *less-than* set ``(-\infty, b]``. We use the shorthand
*function-in-set* to refer to constraints composed of different types of
functions and sets.

This page explains how to write various types of constraints in JuMP. For
nonlinear constraints, see [Nonlinear Modeling](@ref) instead.

## Add a constraint

Add a constraint to a JuMP model using the [`@constraint`](@ref) macro. The
syntax to use depends on the type of constraint you wish to add.

### Add a linear constraint

Create linear constraints using the [`@constraint`](@ref) macro:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @constraint(model, c1, sum(x) <= 1)
c1 : x[1] + x[2] + x[3] ≤ 1

julia> @constraint(model, c2, x[1] + 2 * x[3] >= 2)
c2 : x[1] + 2 x[3] ≥ 2

julia> @constraint(model, c3, sum(i * x[i] for i in 1:3) == 3)
c3 : x[1] + 2 x[2] + 3 x[3] = 3

julia> @constraint(model, c4, 4 <= 2 * x[2] <= 5)
c4 : 2 x[2] ∈ [4, 5]
```

### Normalization

JuMP normalizes constraints by moving all of the terms containing variables to
the left-hand side and all of the constant terms to the right-hand side. Thus,
we get:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, c, 2x + 1 <= 4x + 4)
c : -2 x ≤ 3
```

### [Add a quadratic constraint](@id quad_constraints)

In addition to affine functions, JuMP also supports constraints with quadratic
terms. For example:
```jldoctest con_quadratic
julia> model = Model();

julia> @variable(model, x[i=1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @variable(model, t >= 0)
t

julia> @constraint(model, my_q, x[1]^2 + x[2]^2 <= t^2)
my_q : x[1]² + x[2]² - t² ≤ 0
```

!!! tip
    Because solvers can take advantage of the knowledge that a constraint is
    quadratic, prefer adding quadratic constraints using [`@constraint`](@ref),
    rather than [`@NLconstraint`](@ref).

## Vectorized constraints

You can also add constraints to JuMP using vectorized linear algebra. For
example:
```jldoctest con_vector
julia> model = Model();

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

julia> @constraint(model, con_vector, A * x == b)
con_vector : [x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Zeros(2)

julia> @constraint(model, con_scalar, A * x .== b)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 con_scalar : x[1] + 2 x[2] = 5
 con_scalar : 3 x[1] + 4 x[2] = 6
```

The two constraints, `==` and `.==` are similar, but subtly different. The first
creates a single constraint that is a [`MOI.VectorAffineFunction`](@ref) in
[`MOI.Zeros`](@ref) constraint. The second creates a vector of
[`MOI.ScalarAffineFunction`](@ref) in [`MOI.EqualTo`](@ref) constraints.

Which formulation to choose depends on the solver, and what you want to do with
the constraint object `con_vector` or `con_scalar`.

 * If you are using a conic solver, expect the dual of `con_vector` to be a
   `Vector{Float64}`, and do not intend to delete a row in the constraint,
   choose the `==` formulation.
 * If you are using a solver that expects a list of scalar constraints, for
   example HiGHS, or you wish to delete part of the constraint or access a
   single row of the constraint, for example, `dual(con_scalar[2])`, then use
   the broadcast `.==`.

JuMP reformulates both constraints into the other form if needed by the solver,
but choosing the right format for a particular solver is more efficient.

You can also use `<=`, `.<=` , `>=`, and `.>=` as comparison operators in the
constraint.

```jldoctest con_vector
julia> @constraint(model, A * x <= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Nonpositives(2)

julia> @constraint(model, A * x .<= b)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x[1] + 2 x[2] ≤ 5
 3 x[1] + 4 x[2] ≤ 6

julia> @constraint(model, A * x >= b)
[x[1] + 2 x[2] - 5, 3 x[1] + 4 x[2] - 6] ∈ MathOptInterface.Nonnegatives(2)

julia> @constraint(model, A * x .>= b)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.GreaterThan{Float64}}, ScalarShape}}:
 x[1] + 2 x[2] ≥ 5
 3 x[1] + 4 x[2] ≥ 6
```

### Vectorized matrix constraints

In most cases, you cannot use the non-broadcasting syntax for general matrices.
For example:

```jldoctest
julia> model = Model();

julia> @variable(model, X[1:2, 1:2])
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[2,1]  X[2,2]

julia> @constraint(model, X >= 0)
ERROR: At none:1: `@constraint(model, X >= 0)`: Unsupported matrix in vector-valued set. Did you mean to use the broadcasting syntax `.>=` instead? Alternatively, perhaps you are missing a set argument like `@constraint(model, X >= 0, PSDCone())` or `@constraint(model, X >= 0, HermmitianPSDCone())`.
Stacktrace:
[...]
```

Instead, to represent matrix inequalities you must always use the element-wise
broadcasting `.==`, `.>=`, or `.<=`, or use the [Set inequality syntax](@ref).

There are two exceptions: if the result of the left-hand side minus the
right-hand side is a `LinearAlgebra.Symmetric` matrix or a `LinearAlgebra.Hermitian`
matrix, you may use the non-broadcasting equality syntax:

```jldoctest con_symmetric_zeros
julia> using LinearAlgebra

julia> model = Model();

julia> @variable(model, X[1:2, 1:2], Symmetric)
2×2 Symmetric{VariableRef, Matrix{VariableRef}}:
 X[1,1]  X[1,2]
 X[1,2]  X[2,2]

julia> @constraint(model, X == LinearAlgebra.I)
[X[1,1] - 1  X[1,2];
 X[1,2]      X[2,2] - 1] ∈ Zeros()
```

Despite the model showing the matrix in [`Zeros`](@ref), this will add only
three rows to the constraint matrix because the symmetric constraints are
redundant. In contrast, the broadcasting syntax adds four linear constraints:

```jldoctest con_symmetric_zeros
julia> @constraint(model, X .== LinearAlgebra.I)
2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.EqualTo{Float64}}, ScalarShape}}:
 X[1,1] = 1  X[1,2] = 0
 X[1,2] = 0  X[2,2] = 1
```

The same holds for `LinearAlgebra.Hermitian` matrices:

```jldoctest con_hermitian_zeros
julia> using LinearAlgebra

julia> model = Model();

julia> @variable(model, X[1:2, 1:2] in HermitianPSDCone())
2×2 Hermitian{GenericAffExpr{ComplexF64, VariableRef}, Matrix{GenericAffExpr{ComplexF64, VariableRef}}}:
 real(X[1,1])                    real(X[1,2]) + imag(X[1,2]) im
 real(X[1,2]) - imag(X[1,2]) im  real(X[2,2])

julia> @constraint(model, X == LinearAlgebra.I)
[real(X[1,1]) - 1                real(X[1,2]) + imag(X[1,2]) im;
 real(X[1,2]) - imag(X[1,2]) im  real(X[2,2]) - 1] ∈ Zeros()

julia> @constraint(model, X .== LinearAlgebra.I)
2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{ComplexF64}, MathOptInterface.EqualTo{ComplexF64}}, ScalarShape}}:
 real(X[1,1]) = 1                    real(X[1,2]) + imag(X[1,2]) im = 0
 real(X[1,2]) - imag(X[1,2]) im = 0  real(X[2,2]) = 1
```

## Containers of constraints

The [`@constraint`](@ref) macro supports creating collections of constraints.
We'll cover some brief syntax here; read the [Constraint containers](@ref)
section for more details:

Create arrays of constraints:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @constraint(model, c[i=1:3], x[i] <= i^2)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c[1] : x[1] ≤ 1
 c[2] : x[2] ≤ 4
 c[3] : x[3] ≤ 9

julia> c[2]
c[2] : x[2] ≤ 4
```

Sets can be any Julia type that supports iteration:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @constraint(model, c[i=2:3, ["red", "blue"]], x[i] <= i^2)
2-dimensional DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, ["red", "blue"]
And data, a 2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c[2,red] : x[2] ≤ 4  c[2,blue] : x[2] ≤ 4
 c[3,red] : x[3] ≤ 9  c[3,blue] : x[3] ≤ 9

julia> c[2, "red"]
c[2,red] : x[2] ≤ 4
```

Sets can depend upon previous indices:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @constraint(model, c[i=1:3, j=i:3], x[i] <= j)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  c[1,1] : x[1] ≤ 1
  [1, 2]  =  c[1,2] : x[1] ≤ 2
  [1, 3]  =  c[1,3] : x[1] ≤ 3
  [2, 2]  =  c[2,2] : x[2] ≤ 2
  [2, 3]  =  c[2,3] : x[2] ≤ 3
  [3, 3]  =  c[3,3] : x[3] ≤ 3
```
and you can filter elements in the sets using the `;` syntax:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:9]);

julia> @constraint(model, c[i=1:9; mod(i, 3) == 0], x[i] <= i)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 1, Tuple{Int64}} with 3 entries:
  [3]  =  c[3] : x[3] ≤ 3
  [6]  =  c[6] : x[6] ≤ 6
  [9]  =  c[9] : x[9] ≤ 9
```

## Registered constraints

When you create constraints, JuMP registers them inside the model using their
corresponding symbol. Get a registered name using `model[:key]`:
```jldoctest
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x)
x

julia> @constraint(model, my_c, 2x <= 1)
my_c : 2 x ≤ 1

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: my_c, x

julia> model[:my_c] === my_c
true
```

## Anonymous constraints

To reduce the likelihood of accidental bugs, and because JuMP registers
constraints inside a model, creating two constraints with the same name is an
error:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, c, 2x <= 1)
c : 2 x ≤ 1

julia> @constraint(model, c, 2x <= 1)
ERROR: An object of name c is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :c)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:c]`.
[...]
```

A common reason for encountering this error is adding constraints in a loop.

As a work-around, JuMP provides *anonymous* constraints. Create an anonymous
constraint by omitting the name argument:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> c = @constraint(model, 2x <= 1)
2 x ≤ 1
```

Create a container of anonymous constraints by dropping the name in front of
the `[`:
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> c = @constraint(model, [i = 1:3], x[i] <= i)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x[1] ≤ 1
 x[2] ≤ 2
 x[3] ≤ 3
```

## Constraint names

In addition to the symbol that constraints are registered with, constraints have
a `String` name that is used for printing and writing to file formats.

Get and set the name of a constraint using [`name(::JuMP.ConstraintRef)`](@ref)
and [`set_name(::JuMP.ConstraintRef, ::String)`](@ref):
```jldoctest
julia> model = Model(); @variable(model, x);

julia> @constraint(model, con, x <= 1)
con : x ≤ 1

julia> name(con)
"con"

julia> set_name(con, "my_con_name")

julia> con
my_con_name : x ≤ 1
```

Override the default choice of name using the `base_name` keyword:
```jldoctest constraint_name_vector
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1
 my_con[2] : x ≤ 2
```

Note that names apply to each element of the container, not to the container of
constraints:
```jldoctest constraint_name_vector
julia> name(con[1])
"my_con[1]"

julia> set_name(con[1], "c")

julia> con
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 c : x ≤ 1
 my_con[2] : x ≤ 2
```

!!! tip
    For some models, setting the string name of each constraint can take a
    non-trivial portion of the total time required to build the model. Turn off
    `String` names by passing `set_string_name = false` to [`@constraint`](@ref):
    ```jldoctest
    julia> model = Model();

    julia> @variable(model, x);

    julia> @constraint(model, con, x <= 2, set_string_name = false)
    x ≤ 2
    ```
    See [Disable string names](@ref) for more information.

### Retrieve a constraint by name

Retrieve a constraint from a model using [`constraint_by_name`](@ref):
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "c")
c : x ≤ 1
```

If the name is not present, `nothing` will be returned:
```jldoctest constraint_name_vector
julia> constraint_by_name(model, "bad_name")
```

You can only look up individual constraints using [`constraint_by_name`](@ref).
Something like this will not work:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> con = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1
 my_con[2] : x ≤ 2

julia> constraint_by_name(model, "my_con")
```

To look up a collection of constraints, do not use [`constraint_by_name`](@ref).
Instead, register them using the `model[:key] = value` syntax:
```jldoctest
julia> model = Model(); @variable(model, x);

julia> model[:con] = @constraint(model, [i=1:2], x <= i, base_name = "my_con")
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1
 my_con[2] : x ≤ 2

julia> model[:con]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 my_con[1] : x ≤ 1
 my_con[2] : x ≤ 2
```

## String names, symbolic names, and bindings

It's common for new users to experience confusion relating to constraints.
Part of the problem is the difference between the name that a constraint is
registered under and the `String` name used for printing.

Here's a summary of the differences:

 * Constraints are created using [`@constraint`](@ref).
 * Constraints can be named or anonymous.
 * Named constraints have the form `@constraint(model, c, expr)`. For named
   constraints:
   * The `String` name of the constraint is set to `"c"`.
   * A Julia variable `c` is created that binds `c` to  the JuMP constraint.
   * The name `:c` is registered as a key in the model with the value `c`.
 * Anonymous constraints have the form `c = @constraint(model, expr)`. For
   anonymous constraints:
   * The `String` name of the constraint is set to `""`.
   * You control the name of the Julia variable used as the binding.
   * No name is registered as a key in the model.
 * The `base_name` keyword can override the `String` name of the constraint.
 * You can manually register names in the model via `model[:key] = value`.

Here's an example of the differences:
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> c_binding = @constraint(model, 2x <= 1, base_name = "c")
c : 2 x ≤ 1

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x

julia> c
ERROR: UndefVarError: c not defined

julia> c_binding
c : 2 x ≤ 1

julia> name(c_binding)
"c"

julia> model[:c_register] = c_binding
c : 2 x ≤ 1

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
`AffExpr`-in-`MathOptInterface.LessThan{Float64}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: c_register, x

julia> model[:c_register]
c : 2 x ≤ 1

julia> model[:c_register] === c_binding
true

julia> c
ERROR: UndefVarError: c not defined
```

## The `@constraints` macro

If you have many [`@constraint`](@ref) calls, use the [`@constraints`](@ref)
macro to improve readability:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraints(model, begin
           2x <= 1
           c, x >= -1
       end)
(2 x ≤ 1, c : x ≥ -1)

julia> print(model)
Feasibility
Subject to
 c : x ≥ -1
 2 x ≤ 1
```
The [`@constraints`](@ref) macro returns a tuple of the constraints that were
defined.

## [Duality](@id constraint_duality)

JuMP adopts the notion of [conic duality from MathOptInterface](@ref Duality).
For linear programs, a feasible dual on a `>=` constraint is nonnegative and a
feasible dual on a `<=` constraint is nonpositive. If the constraint is an
equality constraint, it depends on which direction is binding.

!!! warning
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.** If you have a linear program, and you
    want the textbook definition, you probably want to use [`shadow_price`](@ref)
    and [`reduced_cost`](@ref) instead.

The dual value associated with a constraint in the most recent solution can be
accessed using the [`dual`](@ref) function. Use [`has_duals`](@ref) to check if
the model has a dual solution available to query. For example:
```jldoctest con_duality
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x)
x

julia> @constraint(model, con, x <= 1)
con : x ≤ 1

julia> @objective(model, Min, -2x)
-2 x

julia> has_duals(model)
false

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

To help users who may be less familiar with conic duality, JuMP provides
[`shadow_price`](@ref), which returns a value that can be interpreted as the
improvement in the objective in response to an infinitesimal relaxation (on the
scale of one unit) in the right-hand side of the constraint.
[`shadow_price`](@ref) can be used only on linear constraints with a `<=`, `>=`,
or `==` comparison operator.

In the example above, `dual(con)` returned `-2.0` regardless of the optimization
sense. However, in the second case when the optimization sense is `Max`,
[`shadow_price`](@ref) returns:
```jldoctest con_duality
julia> shadow_price(con)
2.0
```

### Duals of variable bounds

To query the dual variables associated with a variable bound, first obtain a
constraint reference using one of [`UpperBoundRef`](@ref),
[`LowerBoundRef`](@ref), or [`FixRef`](@ref), and then call [`dual`](@ref) on
the returned constraint reference. The [`reduced_cost`](@ref) function may
simplify this process as it returns the shadow price of an active bound of
a variable (or zero, if no active bound exists).
```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x <= 1)
x

julia> @objective(model, Min, -2x)
-2 x

julia> optimize!(model)

julia> dual(UpperBoundRef(x))
-2.0

julia> reduced_cost(x)
-2.0
```

## Modify a constant term

This section explains how to modify the constant term in a constraint. There are
multiple ways to achieve this goal; we explain three options.

### Option 1: change the right-hand side

Use [`set_normalized_rhs`](@ref) to modify the right-hand side (constant)
term of a linear or quadratic  constraint. Use [`normalized_rhs`](@ref) to query
the right-hand side term.
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 2x <= 1)
con : 2 x ≤ 1

julia> set_normalized_rhs(con, 3)

julia> con
con : 2 x ≤ 3

julia> normalized_rhs(con)
3.0
```

!!! warning
    [`set_normalized_rhs`](@ref) sets the right-hand side term of the
    normalized constraint. See [Normalization](@ref) for more details.

### Option 2: use fixed variables

If constraints are complicated, for example, they are composed of a number of
components, each of which has a constant term, then it may be difficult to
calculate what the right-hand side term is in the standard form.

For this situation, JuMP includes the ability to *fix* variables to a
value using the [`fix`](@ref) function. Fixing a variable sets its lower
and upper bound to the same value. Thus, changes in a constant term can be
simulated by adding a new variable and fixing it to different values. Here is
an example:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, const_term)
const_term

julia> @constraint(model, con, 2x <= const_term + 1)
con : 2 x - const_term ≤ 1

julia> fix(const_term, 1.0)
```
The constraint `con` is now equivalent to `2x <= 2`.

!!! warning
    Fixed variables are not replaced with constants when communicating the
    problem to a solver. Therefore, even though `const_term` is fixed, it is
    still a decision variable, and so `const_term * x` is bilinear.

### Option 3: modify the function's constant term

The third option is to use [`add_to_function_constant`](@ref). The constant
given is added to the function of a `func`-in-`set` constraint. In the following
example, adding `2` to the function has the effect of removing `2` to the
right-hand side:
```jldoctest con_add
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 2x <= 1)
con : 2 x ≤ 1

julia> add_to_function_constant(con, 2)

julia> con
con : 2 x ≤ -1

julia> normalized_rhs(con)
-1.0
```

In the case of interval constraints, the constant is removed from each bound:
```jldoctest con_add_interval
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 0 <= 2x + 1 <= 2)
con : 2 x ∈ [-1, 1]

julia> add_to_function_constant(con, 3)

julia> con
con : 2 x ∈ [-4, -2]
```

## Modify a variable coefficient

### Scalar constraints

To modify the coefficients for a linear term (modifying the coefficient of a
quadratic term is not supported) in a constraint, use
[`set_normalized_coefficient`](@ref). To query the current coefficient, use
[`normalized_coefficient`](@ref).
```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @constraint(model, con, 2x[1] + x[2] <= 1)
con : 2 x[1] + x[2] ≤ 1

julia> set_normalized_coefficient(con, x[2], 0)

julia> con
con : 2 x[1] ≤ 1

julia> normalized_coefficient(con, x[2])
0.0
```

!!! warning
    [`set_normalized_coefficient`](@ref) sets the coefficient of the normalized
    constraint. See [Normalization](@ref) for more details.

### Vector constraints

To modify the coefficients of a vector-valued constraint, use
[`set_normalized_coefficients`](@ref).
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, [2x + 3x, 4x] in MOI.Nonnegatives(2))
con : [5 x, 4 x] ∈ MathOptInterface.Nonnegatives(2)

julia> set_normalized_coefficients(con, x, [(1, 3.0)])

julia> con
con : [3 x, 4 x] ∈ MathOptInterface.Nonnegatives(2)

julia> set_normalized_coefficients(con, x, [(1, 2.0), (2, 5.0)])

julia> con
con : [2 x, 5 x] ∈ MathOptInterface.Nonnegatives(2)
```

## Delete a constraint

Use [`delete`](@ref) to delete a constraint from a model. Use [`is_valid`](@ref)
to check if a constraint belongs to a model and has not been deleted.

```jldoctest constraints_delete
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con, 2x <= 1)
con : 2 x ≤ 1

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
con : 2 x ≤ 1
```

!!! info
    [`delete`](@ref) does not automatically [`unregister`](@ref) because we do
    not distinguish between names that are automatically registered by JuMP
    macros, and names that are manually registered by the user by setting values
    in [`object_dictionary`](@ref). In addition, deleting a constraint and then
    adding a new constraint of the same name is an easy way to introduce bugs
    into your code.

## Start values

Provide a starting value (also called warmstart) for a constraint's primal and
dual solutions using [`set_start_value`](@ref) and
[`set_dual_start_value`](@ref).

Query the starting value for a constraint's primal and dual solution using
[`start_value`](@ref) and [`dual_start_value`](@ref). If no start value has been
set, the methods will return `nothing`.

```jldoctest constraint_dual_start
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, con, x >= 10)
con : x ≥ 10

julia> start_value(con)

julia> set_start_value(con, 10.0)

julia> start_value(con)
10.0

julia> dual_start_value(con)

julia> set_dual_start_value(con, 2)

julia> dual_start_value(con)
2.0
```

Vector-valued constraints require a vector:
```jldoctest constraint_dual_start_vector
julia> model = Model();

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

!!! tip
    To simplify setting start values for all variables and constraints in a
    model, see [`set_start_values`](@ref). The [Primal and dual warm-starts](@ref)
    tutorial also gives a detailed description of how to iterate over constraints
    in the model to set custom start values.

## Constraint containers

Like [Variable containers](@ref), JuMP provides a mechanism for building groups
of constraints compactly. References to these groups of constraints are returned
in *containers*. Three types of constraint containers are supported: `Array`s,
`DenseAxisArray`s, and `SparseAxisArray`s. We explain each of these in the
following.

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### [Arrays](@id constraint_arrays)

One way of adding a group of constraints compactly is the following:
```jldoctest constraint_arrays
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con[i = 1:3], i * x <= i + 1)
3-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1] : x ≤ 2
 con[2] : 2 x ≤ 3
 con[3] : 3 x ≤ 4
```
JuMP returns references to the three constraints in an `Array` that is bound to
the Julia variable `con`. This array can be accessed and sliced as you would
with any Julia array:
```jldoctest constraint_arrays
julia> con[1]
con[1] : x ≤ 2

julia> con[2:3]
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[2] : 2 x ≤ 3
 con[3] : 3 x ≤ 4
```

Anonymous containers can also be constructed by dropping the name (for example, `con`)
before the square brackets:
```jldoctest constraint_arrays
julia> con = @constraint(model, [i = 1:2], i * x <= i + 1)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 x ≤ 2
 2 x ≤ 3
```

Just like [`@variable`](@ref), JuMP will form an `Array` of constraints when it
can determine at parse time that the indices are one-based integer ranges.
Therefore `con[1:b]` will create an `Array`, but `con[a:b]` will not. A special
case is `con[Base.OneTo(n)]` which will produce an `Array`. If JuMP cannot
determine that the indices are one-based integer ranges (for example, in the case of
`con[a:b]`), JuMP will create a `DenseAxisArray` instead.

### DenseAxisArrays

The syntax for constructing a [`DenseAxisArray`](@ref Containers.DenseAxisArray)
of constraints is very similar to the
[syntax for constructing](@ref variable_jump_arrays) a `DenseAxisArray` of
variables.

```jldoctest constraint_jumparrays
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con[i = 1:2, j = 2:3], i * x <= j + 1)
2-dimensional DenseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape},2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, 2:3
And data, a 2×2 Matrix{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}:
 con[1,2] : x ≤ 3    con[1,3] : x ≤ 4
 con[2,2] : 2 x ≤ 3  con[2,3] : 2 x ≤ 4
```

### SparseAxisArrays

The syntax for constructing a
[`SparseAxisArray`](@ref Containers.SparseAxisArray) of constraints is very
similar to the [syntax for constructing](@ref variable_sparseaxisarrays) a
`SparseAxisArray` of variables.

```jldoctest constraint_jumparrays
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, con[i = 1:2, j = 1:2; i != j], i * x <= j + 1)
JuMP.Containers.SparseAxisArray{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}, 2, Tuple{Int64, Int64}} with 2 entries:
  [1, 2]  =  con[1,2] : x ≤ 3
  [2, 1]  =  con[2,1] : 2 x ≤ 2
```

!!! warning
    If you have many index dimensions and a large amount of sparsity, read
    [Performance considerations](@ref).

### Forcing the container type

When creating a container of constraints, JuMP will attempt to choose the
tightest container type that can store the constraints. However, because this
happens at parse time, it does not always make the best choice. Just like in
[`@variable`](@ref), you can force the type of container using the `container`
keyword. For syntax and the reason behind this, take a look at the
[variable docs](@ref variable_forcing).

### Constraints with similar indices

Containers are often used to create constraints over a set of indices. However,
you'll often have cases in which you are repeating the indices:
```jldoctest repeated_containers
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> @variable(model, y[1:2]);

julia> @constraints(model, begin
           [i=1:2, j=1:2, k=1:2], i * x[j] <= k
           [i=1:2, j=1:2, k=1:2], i * y[j] <= k
       end);
```
This is hard to read and leads to a lot of copy-paste. A more readable way is to
use a for-loop:
```jldoctest repeated_containers
julia> for i=1:2, j=1:2, k=1:2
           @constraints(model, begin
               i * x[j] <= k
               i * y[j] <= k
           end)
       end
```

## Accessing constraints from a model

Query the types of function-in-set constraints in a model using
[`list_of_constraint_types`](@ref):
```jldoctest con_access
julia> model = Model();

julia> @variable(model, x[i=1:2] >= i, Int);

julia> @constraint(model, x[1] + x[2] <= 1);

julia> list_of_constraint_types(model)
3-element Vector{Tuple{Type, Type}}:
 (AffExpr, MathOptInterface.LessThan{Float64})
 (VariableRef, MathOptInterface.GreaterThan{Float64})
 (VariableRef, MathOptInterface.Integer)
```

For a given combination of function and set type, use
[`num_constraints`](@ref) to access the number of constraints and
[`all_constraints`](@ref) to access a list of their references:
```jldoctest con_access
julia> num_constraints(model, VariableRef, MOI.Integer)
2

julia> cons = all_constraints(model, VariableRef, MOI.Integer)
2-element Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex, MathOptInterface.Integer}, ScalarShape}}:
 x[1] integer
 x[2] integer
```

You can also count the total number of constraints in the model, but you must
explicitly choose whether to count `VariableRef` constraints such as bound and
integrality constraints:
```jldoctest con_access
julia> num_constraints(model; count_variable_in_set_constraints = true)
5

julia> num_constraints(model; count_variable_in_set_constraints = false)
1
```
The same also applies for [`all_constraints`](@ref):
```jldoctest con_access
julia> all_constraints(model; include_variable_in_set_constraints = true)
5-element Vector{ConstraintRef}:
 x[1] + x[2] ≤ 1
 x[1] ≥ 1
 x[2] ≥ 2
 x[1] integer
 x[2] integer

julia> all_constraints(model; include_variable_in_set_constraints = false)
1-element Vector{ConstraintRef}:
 x[1] + x[2] ≤ 1
```

If you need finer-grained control on which constraints to include, use a variant
of:
```jldoctest con_access
julia> sum(
           num_constraints(model, F, S) for
           (F, S) in list_of_constraint_types(model) if F != VariableRef
       )
1
```

Use [`constraint_object`](@ref) to get an instance of an
[`AbstractConstraint`](@ref) object that stores the constraint data:
```jldoctest con_access
julia> con = constraint_object(cons[1])
ScalarConstraint{VariableRef, MathOptInterface.Integer}(x[1], MathOptInterface.Integer())

julia> con.func
x[1]

julia> con.set
MathOptInterface.Integer()
```

## MathOptInterface constraints

Because JuMP is based on MathOptInterface, you can add any constraints supported
by MathOptInterface using the function-in-set syntax. For a list of supported
functions and sets, read [Standard form problem](@ref).

!!! note
    We use `MOI` as an alias for the `MathOptInterface` module. This alias is
    defined by `using JuMP`. You may also define it in your code as follows:
    ```julia
    import MathOptInterface as MOI
    ```

For example, the following two constraints are equivalent:
```jldoctest moi
julia> model = Model();

julia> @variable(model, x[1:3]);

julia> @constraint(model, 2 * x[1] <= 1)
2 x[1] ≤ 1

julia> @constraint(model, 2 * x[1] in MOI.LessThan(1.0))
2 x[1] ≤ 1
```

You can also use any set defined by MathOptInterface:
```jldoctest moi
julia> @constraint(model, x - [1; 2; 3] in MOI.Nonnegatives(3))
[x[1] - 1, x[2] - 2, x[3] - 3] ∈ MathOptInterface.Nonnegatives(3)

julia> @constraint(model, x in MOI.ExponentialCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.ExponentialCone()
```

!!! info
    Similar to how JuMP defines the `<=` and `>=` syntax as a convenience way to
    specify [`MOI.LessThan`](@ref) and [`MOI.GreaterThan`](@ref) constraints,
    the remaining sections in this page describe functions and syntax that have
    been added for the convenience of common modeling situations.

## Set inequality syntax

For modeling convenience, the syntax `@constraint(model, x >= y, Set())` is
short-hand for `@constraint(model, x - y in Set())`. Therefore, the following
calls are equivalent:
```jldoctest set_inequality
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> y = [0.5, 0.75];

julia> @constraint(model, x >= y, MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)

julia> @constraint(model, y <= x, MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)

julia> @constraint(model, x - y in MOI.Nonnegatives(2))
[x[1] - 0.5, x[2] - 0.75] ∈ MathOptInterface.Nonnegatives(2)
```

Non-zero constants are not supported in this syntax:
```jldoctest set_inequality
julia> @constraint(model, x >= 1, MOI.Nonnegatives(2))
ERROR: Operation `sub_mul` between `Vector{VariableRef}` and `Int64` is not allowed. This most often happens when you write a constraint like `x >= y` where `x` is an array and `y` is a constant. Use the broadcast syntax `x .- y >= 0` instead.
Stacktrace:
[...]
```
Use instead:
```jldoctest set_inequality
julia> @constraint(model, x .- 1 >= 0, MOI.Nonnegatives(2))
[x[1] - 1, x[2] - 1] ∈ MathOptInterface.Nonnegatives(2)
```

## Second-order cone constraints

A [`SecondOrderCone`](@ref) constrains the variables `t` and `x` to the set:
```math
||x||_2 \le t,
```
and ``t \ge 0``. It can be added as follows:
```jldoctest
julia> model = Model();

julia> @variable(model, t)
t

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, [t; x] in SecondOrderCone())
[t, x[1], x[2]] ∈ MathOptInterface.SecondOrderCone(3)
```

## Rotated second-order cone constraints

A [`RotatedSecondOrderCone`](@ref) constrains the variables `t`, `u`, and `x`
to the set:
```math
||x||_2^2 \le 2 t \cdot u
```
and ``t, u \ge 0``. It can be added as follows:
```jldoctest
julia> model = Model();

julia> @variable(model, t)
t

julia> @variable(model, u)
u

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> @constraint(model, [t; u; x] in RotatedSecondOrderCone())
[t, u, x[1], x[2]] ∈ MathOptInterface.RotatedSecondOrderCone(4)
```

## Semi-integer and semi-continuous variables

Semi-continuous variables are constrained to the set
``x \in \{0\} \cup [l, u]``.

Create a semi-continuous variable using the [`Semicontinuous`](@ref) set:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, x in Semicontinuous(1.5, 3.5))
x in MathOptInterface.Semicontinuous{Float64}(1.5, 3.5)
```

Semi-integer variables  are constrained to the set
``x \in \{0\} \cup \{l, l+1, \dots, u\}``.

Create a semi-integer variable using the [`Semiinteger`](@ref) set:
```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @constraint(model, x in Semiinteger(1.0, 3.0))
x in MathOptInterface.Semiinteger{Float64}(1.0, 3.0)
```

## Special Ordered Sets of Type 1

In a Special Ordered Set of Type 1 (often denoted SOS-I or SOS1), at most one
element can take a non-zero value.

Construct SOS-I constraints using the [`SOS1`](@ref) set:
```jldoctest con_sos
julia> model = Model();

julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SOS1())
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([1.0, 2.0, 3.0])
```

Although not required for feasibility, solvers can benefit from an ordering of
the variables (for example, the variables represent different factories to build, at
most one factory can be built, and the factories can be ordered according to
cost). To induce an ordering, a vector of weights can be provided, and the
variables are ordered according to their corresponding weight.

For example, in the constraint:
```jldoctest con_sos
julia> @constraint(model, x in SOS1([3.1, 1.2, 2.3]))
[x[1], x[2], x[3]] in MathOptInterface.SOS1{Float64}([3.1, 1.2, 2.3])
```
the variables `x` have precedence `x[2]`, `x[3]`, `x[1]`.

## Special Ordered Sets of Type 2

In a Special Ordered Set of Type 2 (SOS-II), at most two elements can be
non-zero, and if there are two non-zeros, they must be consecutive according to
the ordering induced by a weight vector.

Construct SOS-II constraints using the [`SOS2`](@ref) set:
```jldoctest con_sos
julia> @constraint(model, x in SOS2([3.0, 1.0, 2.0]))
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([3.0, 1.0, 2.0])
```
The possible non-zero pairs are (`x[1]`, `x[3]`) and (`x[2]`, `x[3]`):

If the weight vector is omitted, JuMP induces an ordering from `1:length(x)`:
```jldoctest con_sos
julia> @constraint(model, x in SOS2())
[x[1], x[2], x[3]] in MathOptInterface.SOS2{Float64}([1.0, 2.0, 3.0])
```

## Indicator constraints

Indicator constraints consist of a binary variable and a linear constraint. The
constraint holds when the binary variable takes the value `1`. The constraint
may or may not hold when the binary variable takes the value `0`.

To enforce the constraint `x + y <= 1` when the binary variable `a` is `1`, use:
```jldoctest indicator
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, y)
y

julia> @variable(model, a, Bin)
a

julia> @constraint(model, a => {x + y <= 1})
a => {x + y ≤ 1}
```

If the constraint must hold when `a` is zero, add `!` or `¬` before the binary
variable;
```jldoctest indicator
julia> @constraint(model, !a => {x + y <= 1})
!a => {x + y ≤ 1}
```

## Semidefinite constraints

To constrain a matrix to be positive semidefinite (PSD), use [`PSDCone`](@ref):
```jldoctest con_psd
julia> model = Model();

julia> @variable(model, X[1:2, 1:2])
2×2 Matrix{VariableRef}:
 X[1,1]  X[1,2]
 X[2,1]  X[2,2]

julia> @constraint(model, X >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[2,1]  X[2,2]] ∈ PSDCone()
```

!!! tip
    Where possible, prefer constructing a matrix of
    [Semidefinite variables](@ref) using the [`@variable`](@ref) macro, rather
    than adding a constraint like `@constraint(model, X >= 0, PSDCone())`. In
    some solvers, adding the constraint via [`@constraint`](@ref) is less
    efficient, and can result in additional intermediate variables and
    constraints being added to the model.

The inequality `X >= Y` between two square matrices `X` and `Y` is understood as
constraining `X - Y` to be positive semidefinite.
```jldoctest con_psd
julia> Y = [1 2; 2 1]
2×2 Matrix{Int64}:
 1  2
 2  1

julia> @constraint(model, X >= Y, PSDCone())
[X[1,1] - 1  X[1,2] - 2;
 X[2,1] - 2  X[2,2] - 1] ∈ PSDCone()
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

Note that the lower triangular entries are ignored even if they are
different so use it with caution:
```jldoctest con_psd
julia> @constraint(model, LinearAlgebra.Symmetric(X) >= 0, PSDCone())
[X[1,1]  X[1,2];
 X[1,2]  X[2,2]] ∈ PSDCone()
```
(Note the `(2, 1)` element of the constraint is `X[1,2]`, not `X[2,1]`.)

## Complementarity constraints

A mixed complementarity constraint `F(x) ⟂ x` consists of finding `x` in the
interval `[lb, ub]`, such that the following holds:

- `F(x) == 0` if `lb < x < ub`
- `F(x) >= 0` if `lb == x`
- `F(x) <= 0` if `x == ub`

JuMP supports mixed complementarity constraints via `complements(F(x), x)` or
`F(x) ⟂ x` in the [`@constraint`](@ref) macro. The interval set `[lb, ub]` is
obtained from the variable bounds on `x`.

For example, to define the problem `2x - 1 ⟂ x` with `x ∈ [0, ∞)`, do:
```jldoctest complementarity
julia> model = Model();

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
