```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

# Variables

## What is a JuMP variable?

The term *variable* in mathematical optimization has many meanings. Here, we
distinguish between the following three types of variables:
1. *optimization* variables, which are the mathematical ``x`` in the problem
   ``\max\{f_0(x) | f_i(x) \in S_i\}``.
2. *Julia* variables, which are bindings between a name and a value, for example
   `x = 1`. (See [here](https://docs.julialang.org/en/v1.0.0/manual/variables/)
   for the Julia docs.)
3. *JuMP* variables, which are instances of the `VariableRef` struct
   defined by JuMP that contains a reference to an optimization variable in a
   model. (Extra for experts: the `VariableRef` struct is a thin wrapper around
   a `MOI.VariableIndex`, and also contains a reference to the JuMP model.)

To illustrate these three types of variables, consider the following JuMP code
(the full syntax is explained below):
```jldoctest variables
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, x[1:2])
2-element Array{VariableRef,1}:
 x[1]
 x[2]
```
This code does three things:
1. It adds two *optimization* variables to `model`.
2. It creates two *JuMP* variables that act as references to those optimization
   variables.
3. It binds those JuMP variables as a vector with two elements to the *Julia*
   variable `x`.

To reduce confusion, we will attempt, where possible, to always refer to
variables with their corresponding prefix.

JuMP variables can have attributes, such as names or an initial primal start
value. We illustrate the name attribute in the following example:
```jldoctest variables
julia> @variable(model, y, base_name="decision variable")
decision variable
```
This code does four things:
1. It adds one *optimization* variable to `model`.
2. It creates one *JuMP* variable that acts as a reference to that optimization
   variable.
3. It binds the JuMP variable to the Julia variable `y`.
4. It tells JuMP that the *name* attribute of this JuMP variable is "decision
   variable". JuMP uses the value of `base_name` when it has to print the variable
   as a string.

For example, when we print `y` at the REPL we get:
```jldoctest variables
julia> y
decision variable
```

Because `y` is a Julia variable, we can bind it to a different value. For
example, if we write:
```jldoctest variables
julia> y = 1
1
```
`y` is no longer a binding to a JuMP variable. This does not mean that the JuMP
variable has been destroyed. It still exists and is still a reference to the
same optimization variable. The binding can be reset by querying the model for
the symbol *as it was written in the `@variable` macro*. For example:
```jldoctest variables
julia> model[:y]
decision variable
```

This act of looking up the JuMP variable by using the symbol is most useful when
composing JuMP models across multiple functions, as illustrated by the following
example:
```julia
function add_component_to_model(model::JuMP.Model)
    x = model[:x]
    # ... code that uses x
end
function build_model()
    model = Model()
    @variable(model, x)
    add_component_to_model(model)
end
```

Now that we understand the difference between *optimization*, *JuMP*, and *Julia*
variables, we can introduce more of the functionality of the [`@variable`](@ref)
macro.

## Variable bounds

We have already seen the basic usage of the [`@variable`](@ref) macro. The next
extension is to add lower- and upper-bounds to each optimization variable. This
can be done as follows:
```jldoctest variables_2; setup=:(model=Model())
julia> @variable(model, x_free)
x_free

julia> @variable(model, x_lower >= 0)
x_lower

julia> @variable(model, x_upper <= 1)
x_upper

julia> @variable(model, 2 <= x_interval <= 3)
x_interval

julia> @variable(model, x_fixed == 4)
x_fixed
```
In the above examples, `x_free` represents an unbounded optimization variable,
`x_lower` represents an optimization variable with a lower bound and so forth.

!!! note
    When creating a variable with only a lower-bound or an upper-bound, and the
    value of the bound is not a numeric literal, the name of the variable must appear on the
    left-hand side. Putting the name on the right-hand side will result in an
    error. For example:
    ```julia
    @variable(model, 1 <= x)  # works
    a = 1
    @variable(model, a <= x)  # errors
    @variable(model, x >= a)  # works
    ```

We can query whether an optimization variable has a lower- or upper-bound via
the [`has_lower_bound`](@ref) and [`has_upper_bound`](@ref) functions. For
example:
```jldoctest variables_2
julia> has_lower_bound(x_free)
false

julia> has_upper_bound(x_upper)
true
```

If a variable has a lower or upper bound, we can query the value of it via the
[`lower_bound`](@ref) and [`upper_bound`](@ref) functions. For example:
```jldoctest variables_2
julia> lower_bound(x_interval)
2.0

julia> upper_bound(x_interval)
3.0
```
Querying the value of a bound that does not exist will result in an error.

Instead of using the `<=` and `>=` syntax, we can also use the `lower_bound` and
`upper_bound` keyword arguments. For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, lower_bound=1, upper_bound=2)
x

julia> lower_bound(x)
1.0
```

Another option is to use the [`set_lower_bound`](@ref) and
[`set_upper_bound`](@ref) functions. These can also be used to modify an
existing variable bound. For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x >= 1)
x

julia> lower_bound(x)
1.0

julia> set_lower_bound(x, 2)

julia> lower_bound(x)
2.0
```

We can delete variable bounds using [`delete_lower_bound`](@ref) and
[`delete_upper_bound`](@ref):
```jldoctest; setup=:(model=Model())
julia> @variable(model, 1 <= x <= 2)
x

julia> lower_bound(x)
1.0

julia> delete_lower_bound(x)

julia> has_lower_bound(x)
false

julia> upper_bound(x)
2.0

julia> delete_upper_bound(x)

julia> has_upper_bound(x)
false
```

In addition to upper and lower bounds, JuMP variables can also be fixed to a
value using [`fix`](@ref). See also [`is_fixed`](@ref), [`fix_value`](@ref), and
[`unfix`](@ref).
```jldoctest; setup=:(model=Model())
julia> @variable(model, x == 1)
x

julia> is_fixed(x)
true

julia> fix_value(x)
1.0

julia> unfix(x)

julia> is_fixed(x)
false
```
Fixing a variable with existing bounds will throw an error. To delete the bounds
prior to fixing, use `fix(variable, value; force = true)`.

```jldoctest; setup=:(model=Model())
julia> @variable(model, x >= 1)
x

julia> fix(x, 2)
ERROR: Unable to fix x to 2 because it has existing variable bounds. Consider calling `JuMP.fix(variable, value; force=true)` which will delete existing bounds before fixing the variable.

julia> fix(x, 2; force = true)


julia> fix_value(x)
2.0
```

## Variable names

The name, i.e. the value of the `MOI.VariableName` attribute, of a variable can
be obtained by [`JuMP.name(::JuMP.VariableRef)`](@ref) and set by
[`JuMP.set_name(::JuMP.VariableRef, ::String)`](@ref).
```@docs
name(::JuMP.VariableRef)
set_name(::JuMP.VariableRef, ::String)
```

The variable can also be retrieved from its name using
[`JuMP.variable_by_name`](@ref).
```@docs
variable_by_name
```

## Variable containers

In the examples above, we have mostly created scalar variables. By scalar, we
mean that the Julia variable is bound to exactly one JuMP variable. However, it
is often useful to create collections of JuMP variables inside more complicated
data structures.

JuMP provides a mechanism for creating three types of these data structures,
which we refer to as *containers*. The three types are `Array`s, `DenseAxisArray`s,
and `SparseAxisArray`s. We explain each of these in the following.

### Arrays

We have already seen the creation of an array of JuMP variables with the
`x[1:2]` syntax. This can naturally be extended to create multi-dimensional
arrays of JuMP variables. For example:
```jldoctest variables_arrays; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2])
2×2 Array{VariableRef,2}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]
```

Arrays of JuMP variables can be indexed and sliced as follows:
```jldoctest variables_arrays
julia> x[1, 2]
x[1,2]

julia> x[2, :]
2-element Array{VariableRef,1}:
 x[2,1]
 x[2,2]
```

Variable bounds can depend upon the indices:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:2, j=1:2] >= 2i + j)
2×2 Array{VariableRef,2}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]

julia> lower_bound.(x)
2×2 Array{Float64,2}:
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
    Dimension 2, Symbol[:A, :B]
And data, a 2×2 Array{VariableRef,2}:
 x[1,A]  x[1,B]
 x[2,A]  x[2,B]
```

DenseAxisArrays can be indexed and sliced as follows:
```jldoctest variables_jump_arrays
julia> x[1, :A]
x[1,A]

julia> x[2, :]
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, Symbol[:A, :B]
And data, a 2-element Array{VariableRef,1}:
 x[2,A]
 x[2,B]
```

Similarly to the `Array` case, bounds can depend upon indices. For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=2:3, j=1:2:3] >= 0.5i + j)
2-dimensional DenseAxisArray{VariableRef,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Array{VariableRef,2}:
 x[2,1]  x[2,3]
 x[3,1]  x[3,3]

julia> lower_bound.(x)
2-dimensional DenseAxisArray{Float64,2,...} with index sets:
    Dimension 1, 2:3
    Dimension 2, 1:2:3
And data, a 2×2 Array{Float64,2}:
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
JuMP.Containers.SparseAxisArray{VariableRef,2,Tuple{Int64,Int64}} with 3 entries:
  [1, 2]  =  x[1,2]
  [2, 2]  =  x[2,2]
  [1, 1]  =  x[1,1]
```

We can also conditionally create variables via a JuMP-specific syntax. This
syntax appends a comparison check that depends upon the named indices and is
separated from the indices by a semi-colon (`;`). For example:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[i=1:4; mod(i, 2)==0])
JuMP.Containers.SparseAxisArray{VariableRef,1,Tuple{Int64}} with 2 entries:
  [4]  =  x[4]
  [2]  =  x[2]
```

### [Forcing the container type](@id variable_forcing)

When creating a container of JuMP variables, JuMP will attempt to choose the
tightest container type that can store the JuMP variables. Thus, it will prefer
to create an Array before a DenseAxisArray and a DenseAxisArray before a
SparseAxisArray. However, because this happens at compile time, it does not
always make the best choice. To illustrate this, consider the following example:
```jldoctest variable_force_container; setup=:(model=Model())
julia> A = 1:2
1:2

julia> @variable(model, x[A])
1-dimensional DenseAxisArray{VariableRef,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Array{VariableRef,1}:
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
2-element Array{VariableRef,1}:
 y[1]
 y[2]
```
JuMP now creates a vector of JuMP variables instead of a DenseAxisArray. Note
that choosing an invalid container type will throw an error.

## Integrality shortcuts

Adding integrality constraints to a model such as `@constraint(model, x in MOI.ZeroOne())`
and `@constraint(model, x in MOI.Integer())` is a common operation. Therefore,
JuMP supports two shortcuts for adding such constraints.

#### Binary (ZeroOne) constraints

Binary optimization variables are constrained to the set ``x \in \{0, 1\}``. (The
`MOI.ZeroOne` set in MathOptInterface.) Binary optimization variables can be
created in JuMP by passing `Bin` as an optional positional argument:
```jldoctest variables_binary; setup=:(model=Model())
julia> @variable(model, x, Bin)
x
```
We can check if an optimization variable is binary by calling
[`is_binary`](@ref) on the JuMP variable, and binary constraints can be removed
with [`unset_binary`](@ref).
```jldoctest variables_binary
julia> is_binary(x)
true

julia> unset_binary(x)

julia> is_binary(x)
false
```

Binary optimization variables can also be created by setting the `binary`
keyword to `true`.
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, binary=true)
x
```
#### Integer constraints

Integer optimization variables are constrained to the set ``x \in \mathbb{Z}``.
(The `MOI.Integer` set in MathOptInterface.) Integer optimization variables can
be created in JuMP by passing `Int` as an optional positional argument:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x, Int)
x
```
Integer optimization variables can also be created by setting the `integer`
keyword to `true`.
```jldoctest variables_integer; setup=:(model=Model())
julia> @variable(model, x, integer=true)
x
```
We can check if an optimization variable is integer by calling
[`is_integer`](@ref) on the JuMP variable, and integer constraints can be
removed with [`unset_integer`](@ref).
```jldoctest variables_integer
julia> is_integer(x)
true

julia> unset_integer(x)

julia> is_integer(x)
false
```

## Semidefinite variables

JuMP also supports modeling with semidefinite variables. A square symmetric
matrix ``X`` is positive semidefinite if all eigenvalues are nonnegative. We can
declare a matrix of JuMP variables to be positive semidefinite as follows:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], PSD)
2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```

Note that `x` must be a square 2-dimensional `Array` of JuMP variables; it
cannot be a DenseAxisArray or a SparseAxisArray.
(See [Variable containers](@ref), above, for more on this.)

You can also impose a weaker constraint that the square matrix is only symmetric
(instead of positive semidefinite) as follows:
```jldoctest; setup=:(model=Model())
julia> @variable(model, x[1:2, 1:2], Symmetric)
2×2 LinearAlgebra.Symmetric{VariableRef,Array{VariableRef,2}}:
 x[1,1]  x[1,2]
 x[1,2]  x[2,2]
```
## Anonymous JuMP variables

In all of the above examples, we have created *named* JuMP variables. However,
it is also possible to create so called *anonymous* JuMP variables. To create an
anonymous JuMP variable, we drop the name of the variable from the macro call.
This means dropping the second positional argument if the JuMP variable is a
scalar, or dropping the name before the square bracket (`[`) if a container is
being created. For example:
```jldoctest; setup=:(model=Model())
julia> x = @variable(model)
noname
```
This shows how `(model, x)` is really short for:
```jldoctest anon_variables; setup=:(model=Model())
julia> x = model[:x] = @variable(model, base_name="x")
x
```
An `Array` of anonymous JuMP variables can be created as follows:
```jldoctest; setup=:(model=Model())
julia> y = @variable(model, [i=1:2])
2-element Array{VariableRef,1}:
 noname
 noname
```
If necessary, you can store `x` in `model` as follows:
```
julia> model[:x] = x
```
The `<=` and `>=` short-hand cannot be used to set bounds on anonymous JuMP
variables. Instead, you should use the `lower_bound` and `upper_bound` keywords.

Passing the `Bin` and `Int` variable types are also invalid. Instead, you should
use the `binary` and `integer` keywords.

Thus, the anonymous variant of `@variable(model, x[i=1:2] >= i, Int)` is:
```jldoctest; setup=:(model=Model())
julia> x = @variable(model, [i=1:2], base_name="x", lower_bound=i, integer=true)
2-element Array{VariableRef,1}:
 x[1]
 x[2]
```

!!! warn
    Creating two named JuMP variables with the same name results in an error at
    runtime. Use anonymous variables as an alternative.

## User-defined containers

In the section [Variable containers](@ref), we explained how JuMP supports the
efficient creation of collections of JuMP variables in three types of
containers. However, users are also free to create collections of JuMP variables
in their own datastructures. For example, the following code creates a
dictionary with symmetric matrices as the values:
```jldoctest; setup=:(model=Model())
julia> variables = Dict{Symbol, Array{VariableRef,2}}()
Dict{Symbol,Array{VariableRef,2}} with 0 entries

julia> for key in [:A, :B]
           global variables[key] = @variable(model, [1:2, 1:2])
       end

julia> variables
Dict{Symbol,Array{VariableRef,2}} with 2 entries:
  :A => VariableRef[noname noname; noname noname]
  :B => VariableRef[noname noname; noname noname]
```

## Deleting variables

JuMP supports the deletion of optimization variables.  To delete variables, we
can use the [`delete`](@ref) method. We can also check whether `x` is a valid
JuMP variable in `model` using the [`is_valid`](@ref) method:
```jldoctest variables_delete; setup=:(model=Model())
julia> @variable(model, x)
x

julia> is_valid(model, x)
true

julia> delete(model, x)

julia> is_valid(model, x)
false
```

## Listing all variables

Use [`JuMP.all_variables`](@ref) to obtain a list of all variables present
in the model. This is useful for performing operations like:

- relaxing all integrality constraints in the model
- setting the starting values for variables to the result of the last solve

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

!!! note
    Prior to JuMP 0.19, the previous solution to a solve was automatically set
    as the new starting value. JuMP 0.19 no longer does this automatically. To
    reproduce the functionality, use:
    ```julia
    set_start_value.(all_variables(model), value.(all_variables(model)))
    ```

## [The `@variables` macro](@id variables)

If you have many [`@variable`](@ref) calls, JuMP provides the macro `@variables`
that can improve readability:

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
    Keyword arguments must be contained within parentheses. (See the example
    above.)

## Reference

```@docs
@variable
owner_model
VariableRef
all_variables
num_variables

has_lower_bound
lower_bound
set_lower_bound
delete_lower_bound

has_upper_bound
upper_bound
set_upper_bound
delete_upper_bound

is_fixed
fix_value
fix
unfix

is_integer
set_integer
unset_integer
IntegerRef

is_binary
set_binary
unset_binary
BinaryRef

index(::VariableRef)
optimizer_index(::VariableRef)

set_start_value
start_value
```
