```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

# [Variables](@id jump_variables)

The term *variable* in mathematical optimization has many meanings. For example,
*optimization* variables (also called decision variables) are the unknowns ``x``
that we are solving for in the problem:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x)
    \\
    & \;\;\text{s.t.} & f_i(x) & \in \mathcal{S}_i & i = 1 \ldots m
\end{align}
```

To complicate things, Julia uses *variable* to mean a binding between a name and
a value. For example, in the statement:
```jldoctest
julia> x = 1
1
```
`x` is a variable that stores the value `1`.

JuMP uses *variable* in a third way, to mean an instance of the
[`VariableRef`](@ref) struct. JuMP variables are the link between Julia and the
optimization variables inside a JuMP model.

This page explains how to create and manage JuMP variables in a variety of
contexts.

## Create a variable

Create variables using the [`@variable`](@ref) macro. When creating a variable,
you can also specify variable bounds:
```jldoctest variables_2
model = Model()
@variable(model, x_free)
@variable(model, x_lower >= 0)
@variable(model, x_upper <= 1)
@variable(model, 2 <= x_interval <= 3)
@variable(model, x_fixed == 4)
print(model)

# output

Feasibility
Subject to
 x_fixed = 4.0
 x_lower ≥ 0.0
 x_interval ≥ 2.0
 x_upper ≤ 1.0
 x_interval ≤ 3.0
```

!!! warning
    When creating a variable with a single lower- or upper-bound, and the
    value of the bound is not a numeric literal (e.g., `1` or `1.0`), the name
    of the variable _must_ appear on the left-hand side. Putting the name on the
    right-hand side is an error. For example to create a variable `x`:
    ```julia
    a = 1
    @variable(model, x >= 1)      # ✓ Okay
    @variable(model, 1.0 <= x)    # ✓ Okay
    @variable(model, x >= a)      # ✓ Okay
    @variable(model, a <= x)      # × Not okay
    @variable(model, x >= 1 / 2)  # ✓ Okay
    @variable(model, 1 / 2 <= x)  # × Not okay
    ```

### Containers of variables

The [`@variable`](@ref) macro also supports creating collections of JuMP
variables. We'll cover some brief syntax here; read the [Variable containers](@ref)
section for more details.

You can create arrays of JuMP variables:
```jldoctest; setup=:(model = Model())
julia> @variable(model, x[1:2, 1:2])
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]

julia> x[1, 2]
x[1,2]
```

Index sets can be named, and bounds can depend on those names:
```jldoctest; setup=:(model = Model())
julia> @variable(model, sqrt(i) <= x[i = 1:3] <= i^2)
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> x[2]
x[2]
```

Sets can be any Julia type that supports iteration:
```jldoctest; setup=:(model = Model())
julia> @variable(model, x[i = 2:3, j = 1:2:3, ["red", "blue"]] >= 0)
3-dimensional DenseAxisArray{VariableRef,3,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
    Dimension 3, ["red", "blue"]
And data, a 2×2×2 Array{VariableRef, 3}:
[:, :, "red"] =
 x[2,1,red]  x[2,3,red]
 x[3,1,red]  x[3,3,red]

[:, :, "blue"] =
 x[2,1,blue]  x[2,3,blue]
 x[3,1,blue]  x[3,3,blue]

julia> x[2, 1, "red"]
x[2,1,red]
```

Sets can depend upon previous indices:
```jldoctest; setup=:(model = Model())
julia> @variable(model, u[i = 1:2, j = i:3])
JuMP.Containers.SparseAxisArray{VariableRef, 2, Tuple{Int64, Int64}} with 5 entries:
  [1, 1]  =  u[1,1]
  [1, 2]  =  u[1,2]
  [1, 3]  =  u[1,3]
  [2, 2]  =  u[2,2]
  [2, 3]  =  u[2,3]
```
and we can filter elements in the sets using the `;` syntax:
```jldoctest; setup=:(model = Model())
julia> @variable(model, v[i = 1:9; mod(i, 3) == 0])
JuMP.Containers.SparseAxisArray{VariableRef, 1, Tuple{Int64}} with 3 entries:
  [3]  =  v[3]
  [6]  =  v[6]
  [9]  =  v[9]
```

## Registered variables

When you create variables, JuMP registers them inside the model using their
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

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x

julia> model[:x] === x
true
```

Registered names are most useful when you start to write larger models and
want to break up the model construction into functions:
```jldoctest
julia> function set_objective(model::Model)
           @objective(model, Min, 2 * model[:my_x] + 1)
           return
       end
set_objective (generic function with 1 method)

julia> model = Model();

julia> @variable(model, my_x);

julia> set_objective(model)

julia> print(model)
Min 2 my_x + 1
Subject to
```

## [Anonymous variables](@id anonymous_variables)

To reduce the likelihood of accidental bugs, and because JuMP registers
variables inside a model, creating two variables with the same name is an error:
```julia
julia> model = Model();

julia> @variable(model, x)
x

julia> @variable(model, x)
ERROR: An object of name x is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :x)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:x]`.
[...]
```

A common reason for encountering this error is adding variables in a loop.

As a work-around, JuMP provides *anonymous* variables. Create a scalar valued
anonymous variable by ommitting the name argument:
```jldoctest; setup=:(model=Model())
julia> x = @variable(model)
noname
```
Create a container of anonymous JuMP variables by dropping the name in front of
the `[`:
```jldoctest; setup=:(model=Model())
julia> y = @variable(model, [1:2])
2-element Vector{VariableRef}:
 noname
 noname
```

The `<=` and `>=` short-hand cannot be used to set bounds on scalar-valued
anonymous JuMP variables. Instead, use the `lower_bound` and `upper_bound`
keywords:
```jldoctest; setup=:(model=Model())
julia> x_lower = @variable(model, lower_bound = 1.0)
noname
julia> x_upper = @variable(model, upper_bound = 2.0)
noname

julia> x_interval = @variable(model, lower_bound = 3.0, upper_bound = 4.0)
noname
```

## Variable names

In addition to the symbol that variables are registered with, JuMP variables
have a `String` name that is used for printing and writing to file formats.

Get and set the name of a variable using [`name`](@ref) and [`set_name`](@ref):
```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> name(x)
"x"

julia> set_name(x, "my_x_name")

julia> x
my_x_name
```

Override the default choice of name using the `base_name` keyword:
```jldoctest variable_name_vector
julia> model = Model();

julia> @variable(model, x[i=1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]
```

Note that names apply to each element of the container, not to the container of
variables:
```jldoctest variable_name_vector
julia> name(x[1])
"my_var[1]"

julia> set_name(x[1], "my_x")

julia> x
2-element Vector{VariableRef}:
 my_x
 my_var[2]
```

### Retrieve a variable by name

Retrieve a variable from a model using [`variable_by_name`](@ref):
```jldoctest variable_name_vector
julia> variable_by_name(model, "my_x")
my_x
```
If the name is not present, `nothing` will be returned:
```jldoctest variable_name_vector
julia> variable_by_name(model, "bad_name")
```

You can only look up invididual variables using [`variable_by_name`](@ref).
Something like this will not work:
```jldoctest
julia> model = Model();

julia> @variable(model, [i = 1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]

julia> variable_by_name(model, "my_var")
```

To look up a collection of variables, do not use [`variable_by_name`](@ref).
Instead, register them using the `model[:key] = value` syntax:
```jldoctest
julia> model = Model();

julia> model[:x] = @variable(model, [i = 1:2], base_name = "my_var")
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]

julia> model[:x]
2-element Vector{VariableRef}:
 my_var[1]
 my_var[2]
```

## String names, symbolic names, and bindings

It's common for new users to experience confusion relating to JuMP variables.
Part of the problem the overloaded use of "variable" in mathematical
optimization, along with the difference between the name that a variable is
registered under and the `String` name used for printing.

Here's a summary of the differences:

 * JuMP variables are created using [`@variable`](@ref).
 * JuMP variables can be named or anonymous.
 * Named JuMP variables have the form `@variable(model, x)`. For named
   variables:
   * The `String` name of the variable is set to `"x"`.
   * A Julia variable `x` is created that binds `x` to  the JuMP variable.
   * The name `:x` is registered as a key in the model with the value `x`.
 * Anonymous JuMP variables have the form `x = @variable(model)`. For anonymous
   variables:
   * The `String` name of the variable is set to `""`. When printed, this is
     replaced with `"noname"`.
   * You control the name of the Julia variable used as the binding.
   * No name is registered as a key in the model.
 * The `base_name` keyword can override the `String` name of the variable.
 * You can manually register names in the model via `model[:key] = value`

Here's an example that should make things clearer:
```jldoctest
julia> model = Model();

julia> x_binding = @variable(model, base_name = "x")
x

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> x
ERROR: UndefVarError: x not defined

julia> x_binding
x

julia> name(x_binding)
"x"

julia> model[:x_register] = x_binding
x

julia> model
A JuMP Model
Feasibility problem with:
Variable: 1
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: x_register

julia> model[:x_register]
x

julia> model[:x_register] === x_binding
true

julia> x
ERROR: UndefVarError: x not defined
```

## Create, delete, and modify variable bounds

Query whether a variable has a bound using [`has_lower_bound`](@ref),
[`has_upper_bound`](@ref), and [`is_fixed`](@ref):
```jldoctest variables_2
julia> has_lower_bound(x_free)
false

julia> has_upper_bound(x_upper)
true

julia> is_fixed(x_fixed)
true
```

If a variable has a particular bound, query the value of it using
[`lower_bound`](@ref), [`upper_bound`](@ref), and [`fix_value`](@ref):
```jldoctest variables_2
julia> lower_bound(x_interval)
2.0

julia> upper_bound(x_interval)
3.0

julia> fix_value(x_fixed)
4.0
```

Querying the value of a bound that does not exist will result in an error.

Delete variable bounds using [`delete_lower_bound`](@ref),
[`delete_upper_bound`](@ref), and [`unfix`](@ref):
```jldoctest variables_2
julia> delete_lower_bound(x_lower)

julia> has_lower_bound(x_lower)
false

julia> delete_upper_bound(x_upper)

julia> has_upper_bound(x_upper)
false

julia> unfix(x_fixed)

julia> is_fixed(x_fixed)
false
```

Set or update variable bounds using [`set_lower_bound`](@ref),
[`set_upper_bound`](@ref), and [`fix`](@ref):
```jldoctest variables_2
julia> set_lower_bound(x_lower, 1.1)

julia> set_upper_bound(x_upper, 2.1)

julia> fix(x_fixed, 4.1)
```

Fixing a variable with existing bounds will throw an error. To delete the bounds
prior to fixing, use `fix(variable, value; force = true)`.

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 1)
x

julia> fix(x, 2)
ERROR: Unable to fix x to 2 because it has existing variable bounds. Consider calling `JuMP.fix(variable, value; force=true)` which will delete existing bounds before fixing the variable.

julia> fix(x, 2; force = true)


julia> fix_value(x)
2.0
```

!!! tip
    Use [`fix`](@ref) instead of `@constraint(model, x == 2)`. The former
    modifies variable bounds, while the latter adds a new linear constraint to
    the problem.

## Binary variables

Binary variables are constrained to the set ``x \in \{0, 1\}``.

Create a binary variable by passing `Bin` as an optional positional argument:
```jldoctest variables_binary; setup=:(model=Model())
julia> @variable(model, x, Bin)
x
```

Check if a variable is binary using [`is_binary`](@ref):
```jldoctest variables_binary
julia> is_binary(x)
true
```
Delete a binary constraint using [`unset_binary`](@ref):
```jldoctest variables_binary
julia> unset_binary(x)

julia> is_binary(x)
false
```

Binary variables can also be created by setting the `binary` keyword to `true`:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, binary=true)
x
```
or by using [`set_binary`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x)
x

julia> set_binary(x)
```

## Integer variables

Integer variables are constrained to the set ``x \in \mathbb{Z}``.

Create an integer variable by passing `Int` as an optional positional argument:
```jldoctest variables_integer; setup=:(model=Model())
julia> @variable(model, x, Int)
x
```

Check if a variable is integer using [`is_integer`](@ref):
```jldoctest variables_integer
julia> is_integer(x)
true
```
Delete an integer constraint using [`unset_integer`](@ref).
```jldoctest variables_integer
julia> unset_integer(x)

julia> is_integer(x)
false
```

Integer variables can also be created by setting the `integer` keyword to
`true`:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, integer=true)
x
```
or by using [`set_integer`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x)
x

julia> set_integer(x)
```

!!! tip
    The [`relax_integrality`](@ref) function relaxes all integrality constraints
    in the model, returning a function that can be called to undo the operation
    later on.

## Start values

There are two ways to provide a primal starting solution (also called MIP-start
or a warmstart) for each variable:

 - using the `start` keyword in the [`@variable`](@ref) macro
 - using [`set_start_value`](@ref)

The starting value of a variable can be queried using [`start_value`](@ref). If
no start value has been set, [`start_value`](@ref) will return `nothing`.

```jldoctest variables_start; setup=:(model=Model())
julia> @variable(model, x)
x

julia> start_value(x)

julia> @variable(model, y, start = 1)
y

julia> start_value(y)
1.0

julia> set_start_value(y, 2)

julia> start_value(y)
2.0
```

!!! warning
    Some solvers do not support start values. If a solver does not support start
    values, an `MathOptInterface.UnsupportedAttribute{MathOptInterface.VariablePrimalStart}`
    error will be thrown.

!!! tip
    To set the optimal solution from a previous solve as a new starting value,
    use [`all_variables`](@ref) to get a vector of all the variables in the
    model, then run:
    ```julia
    x = all_variables(model)
    x_solution = value.(x)
    set_start_value.(x, x_solution)
    ```

## [Delete a variable](@id delete_a_variable)

Use [`delete`](@ref) to delete a variable from a model. Use [`is_valid`](@ref)
to check if a variable belongs to a model and has not been deleted.

```jldoctest variables_delete; setup=:(model=Model())
julia> @variable(model, x)
x

julia> is_valid(model, x)
true

julia> delete(model, x)

julia> is_valid(model, x)
false
```

Deleting a variable does not unregister the corresponding name from the model.
Therefore, creating a new variable of the same name will throw an error:
```jldoctest variables_delete
julia> @variable(model, x)
ERROR: An object of name x is already attached to this model. If this
    is intended, consider using the anonymous construction syntax, e.g.,
    `x = @variable(model, [1:N], ...)` where the name of the object does
    not appear inside the macro.

    Alternatively, use `unregister(model, :x)` to first unregister
    the existing name from the model. Note that this will not delete the
    object; it will just remove the reference at `model[:x]`.
[...]
```

After calling [`delete`](@ref), call [`unregister`](@ref) to remove the symbolic
reference:
```jldoctest variables_delete
julia> unregister(model, :x)

julia> @variable(model, x)
x
```

!!! info
    [`delete`](@ref) does not automatically [`unregister`](@ref) because we do
    not distinguish between names that are automatically registered by JuMP
    macros and names that are manually registered by the user by setting values
    in [`object_dictionary`](@ref). In addition, deleting a variable and then
    adding a new variable of the same name is an easy way to introduce bugs into
    your code.

## Variable containers

JuMP provides a mechanism for creating collections of variables in three types
of data structures, which we refer to as *containers*.

The three types are `Array`s, `DenseAxisArray`s, and `SparseAxisArray`s. We
explain each of these in the following.

!!! tip
    You can read more about containers in the [Containers](@ref) section.

### Arrays

We have already seen the creation of an array of JuMP variables with the
`x[1:2]` syntax. This can naturally be extended to create multi-dimensional
arrays of JuMP variables. For example:
```jldoctest variables_arrays; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2])
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]
```

Arrays of JuMP variables can be indexed and sliced as follows:
```jldoctest variables_arrays
julia> x[1, 2]
x[1,2]

julia> x[2, :]
2-element Vector{VariableRef}:
 x[2,1]
 x[2,2]
```

Variable bounds can depend upon the indices:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:2, j=1:2] >= 2i + j)
2×2 Matrix{VariableRef}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]

julia> lower_bound.(x)
2×2 Matrix{Float64}:
 3.0  4.0
 5.0  6.0
```

JuMP will form an `Array` of JuMP variables when it can determine at compile
time that the indices are one-based integer ranges. Therefore `x[1:b]` will
create an `Array` of JuMP variables, but `x[a:b]` will not. If JuMP cannot
determine that the indices are one-based integer ranges (e.g., in the case of
`x[a:b]`), JuMP will create a `DenseAxisArray` instead.

### [DenseAxisArrays](@id variable_jump_arrays)

We often want to create arrays where the indices are not one-based integer
ranges. For example, we may want to create a variable indexed by the name of a
product or a location. The syntax is the same as that above, except with an
arbitrary vector as an index as opposed to a one-based range. The biggest
difference is that instead of returning an `Array` of JuMP variables, JuMP will
return a `DenseAxisArray`. For example:
```jldoctest variables_jump_arrays; setup=:(model=Model())
julia> @variable(model, x[1:2, [:A,:B]])
2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:
    Dimension 1, Base.OneTo(2)
    Dimension 2, [:A, :B]
And data, a 2×2 Matrix{VariableRef}:
 x[1,A]  x[1,B]
 x[2,A]  x[2,B]
```

DenseAxisArrays can be indexed and sliced as follows:
```jldoctest variables_jump_arrays
julia> x[1, :A]
x[1,A]

julia> x[2, :]
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, [:A, :B]
And data, a 2-element Vector{VariableRef}:
 x[2,A]
 x[2,B]
```

Bounds can depend upon indices:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=2:3, j=1:2:3] >= 0.5i + j)
2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Matrix{VariableRef}:
 x[2,1]  x[2,3]
 x[3,1]  x[3,3]

julia> lower_bound.(x)
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Matrix{Float64}:
 2.0  4.0
 2.5  4.5
```

### [SparseAxisArrays](@id variable_sparseaxisarrays)

The third container type that JuMP natively supports is `SparseAxisArray`.
These arrays are created when the indices do not form a rectangular set.
For example, this applies when indices have a dependence upon previous
indices (called *triangular indexing*). JuMP supports this as follows:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:2, j=i:2])
JuMP.Containers.SparseAxisArray{VariableRef, 2, Tuple{Int64, Int64}} with 3 entries:
  [1, 1]  =  x[1,1]
  [1, 2]  =  x[1,2]
  [2, 2]  =  x[2,2]
```

We can also conditionally create variables via a JuMP-specific syntax. This
syntax appends a comparison check that depends upon the named indices and is
separated from the indices by a semi-colon (`;`). For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:4; mod(i, 2)==0])
JuMP.Containers.SparseAxisArray{VariableRef, 1, Tuple{Int64}} with 2 entries:
  [2]  =  x[2]
  [4]  =  x[4]
```

#### Performance considerations

When using the semi-colon as a filter, JuMP iterates over *all* indices and
evaluates the conditional for each combination. If there are many index
dimensions and a large amount of sparsity, this can be inefficient.

For example:
```jldoctest; setup=:(model=Model()), filter=r"[0-9\.]+ seconds.+?time"
julia> N = 10
10

julia> S = [(1, 1, 1), (N, N, N)]
2-element Vector{Tuple{Int64, Int64, Int64}}:
 (1, 1, 1)
 (10, 10, 10)

julia> @time @variable(model, x1[i=1:N, j=1:N, k=1:N; (i, j, k) in S])
  0.203861 seconds (392.22 k allocations: 23.977 MiB, 99.10% compilation time)
JuMP.Containers.SparseAxisArray{VariableRef, 3, Tuple{Int64, Int64, Int64}} with 2 entries:
  [1, 1, 1   ]  =  x1[1,1,1]
  [10, 10, 10]  =  x1[10,10,10]

julia> @time @variable(model, x2[S])
  0.045407 seconds (65.24 k allocations: 3.771 MiB, 99.15% compilation time)
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, [(1, 1, 1), (10, 10, 10)]
And data, a 2-element Vector{VariableRef}:
 x2[(1, 1, 1)]
 x2[(10, 10, 10)]
```

The first option is slower because it is equivalent to:
```julia
x1 = Dict()
for i in 1:N
    for j in 1:N
        for k in 1:N
            if (i, j, k) in S
                x1[i, j, k] = @variable(model)
            end
        end
    end
end
```
If performance is a concern, explicitly construct the set of indices instead of
using the filtering syntax.

### [Forcing the container type](@id variable_forcing)

When creating a container of JuMP variables, JuMP will attempt to choose the
tightest container type that can store the JuMP variables. Thus, it will prefer
to create an Array before a DenseAxisArray and a DenseAxisArray before a
SparseAxisArray. However, because this happens at compile time, JuMP does not
always make the best choice. To illustrate this, consider the following example:
```jldoctest variable_force_container; setup=:(model=Model())
julia> A = 1:2
1:2

julia> @variable(model, x[A])
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Vector{VariableRef}:
 x[1]
 x[2]
```
Since the value (and type) of `A` is unknown at parsing time, JuMP is unable to
infer that `A` is a one-based integer range. Therefore, JuMP creates a
`DenseAxisArray`, even though it could store these two variables in a standard
one-dimensional `Array`.

We can share our knowledge that it is possible to store these JuMP variables as
an array by setting the `container` keyword:
```jldoctest variable_force_container
julia> @variable(model, y[A], container=Array)
2-element Vector{VariableRef}:
 y[1]
 y[2]
```
JuMP now creates a vector of JuMP variables instead of a DenseAxisArray.
Choosing an invalid container type will throw an error.

### User-defined containers

In addition to the built-in container types, you can create your own collections
of JuMP variables.

!!! tip
    This is a point that users often overlook: you are not restricted to the
    built-in container types in JuMP.

For example, the following code creates a dictionary with symmetric matrices as
the values:
```jldoctest; setup=:(model=Model())
julia> variables = Dict{Symbol,Array{VariableRef,2}}(
           key => @variable(model, [1:2, 1:2], Symmetric, base_name = "$(key)")
           for key in [:A, :B]
       )
Dict{Symbol, Matrix{VariableRef}} with 2 entries:
  :A => [A[1,1] A[1,2]; A[1,2] A[2,2]]
  :B => [B[1,1] B[1,2]; B[1,2] B[2,2]]
```

Another common scenario is a request to add variables to existing containers,
for example:
```julia
using JuMP
model = Model()
@variable(model, x[1:2] >= 0)
# Later I want to add
@variable(model, x[3:4] >= 0)
```
This is not possible with the built-in JuMP container types. However, you can
use regular Julia types instead:
```jldoctest
model = Model()
x = model[:x] = @variable(model, [1:2], lower_bound = 0, base_name = "x")
append!(x, @variable(model, [1:2], lower_bound = 0, base_name = "y"))
model[:x]

# output

4-element Vector{VariableRef}:
 x[1]
 x[2]
 y[1]
 y[2]
```

## Semidefinite variables

A square symmetric matrix ``X`` is positive semidefinite if all eigenvalues are
nonnegative.

Declare a matrix of JuMP variables to be positive semidefinite by passing `PSD`
as an optional positional argument:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], PSD)
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

!!! note
    `x` must be a square 2-dimensional `Array` of JuMP variables; it cannot be a
    DenseAxisArray or a SparseAxisArray.

## Symmetric variables

Declare a square matrix of JuMP variables to be symmetric by passing
`Symmetric`  as an optional positional argument:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], Symmetric)
2×2 LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

## Skew-symmetric variables

Declare a matrix of JuMP variables to be skew-symmetric using
[`SkewSymmetricMatrixSpace`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2] in SkewSymmetricMatrixSpace())
2×2 Matrix{AffExpr}:
 0        x[1,2]
 -x[1,2]  0
```
Note that even though `x` is 2 by 2, only one decision variable is added to
`model`; the remaining elements in `x` are linear transformations of the single
variable.

## [The `@variables` macro](@id variables)

If you have many [`@variable`](@ref) calls, JuMP provides the macro
[`@variables`](@ref) that can improve readability:
```jldoctest; setup=:(model=Model())
julia> @variables(model, begin
           x
           y[i=1:2] >= i, (start = i, base_name = "Y_$i")
           z, Bin
       end)

julia> print(model)
Feasibility
Subject to
 Y_1[1] ≥ 1.0
 Y_2[2] ≥ 2.0
 z binary
```

!!! note
    Keyword arguments must be contained within parentheses.

## [Variables constrained on creation](@id jump_variables_on_creation)

All uses of the [`@variable`](@ref) macro documented so far translate into
separate calls for variable creation and the adding of any bound or integrality
constraints.

For example, `@variable(model, x >= 0, Int)`, is equivalent to:
```julia
@variable(model, x)
set_lower_bound(x, 0.0)
set_integer(x)
```
Importantly, the bound and integrality constraints are added _after_ the
variable has been created.

However, some solvers require a constraining set _at creation time_.
We say that these variables are _constrained on creation_.

Use `in` within `@variable` to access the special syntax for constraining
variables on creation. For example, the following creates a vector of variables
constrained on creation to belong to the [`SecondOrderCone`](@ref):
```jldoctest constrained_variables; setup=:(model=Model())
julia> @variable(model, y[1:3] in SecondOrderCone())
3-element Vector{VariableRef}:
 y[1]
 y[2]
 y[3]
```

For contrast, the standard syntax is as follows:
```jldoctest constrained_variables
julia> @variable(model, x[1:3])
3-element Vector{VariableRef}:
 x[1]
 x[2]
 x[3]

julia> @constraint(model, x in SecondOrderCone())
[x[1], x[2], x[3]] ∈ MathOptInterface.SecondOrderCone(3)
```

An alternate syntax to `x in Set` is to use the `set` keyword of
[`@variable`](@ref). This is most useful when creating anonymous variables:
```julia
x = @variable(model, [1:2, 1:2], set = PSDCone())
```

!!! note
    You cannot delete the constraint associated with a variable constrained on
    creation.

### Why use variables constrained on creation?

For most users, it does not matter if you use the constrained on creation
syntax. Therefore, use whatever syntax you find most convenient.

However, if you use [`direct_model`](@ref), you may be forced to use the
constrained on creation syntax.

The technical difference between variables constrained on creation and the
standard JuMP syntax is that variables constrained on creation calls
`MOI.add_constrained_variables`, while the standard JuMP syntax calls
`MOI.add_variables` and then `MOI.add_constraint`.

Consult the implementation of solver package you are using to see if your solver
requires `MOI.add_constrained_variables`.
