```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
```

Variables
=========

What is a JuMP variable?
---------------------------

The term *variable* in computational optimization has many meanings. Here, we
distinguish between the following three types of variables:
1. *optimization* variables, which are the mathematical ``x`` in the problem
   ``\\max\{c^\\top x | Ax = b\}``.
2. *Julia* variables, which are bindings between a name and a value, for example
   `x = 1`. (See [here](https://docs.julialang.org/en/stable/manual/variables/)
   for the Julia docs.)
3. *JuMP* variables, which are instances of the `JuMP.VariableRef` struct
   defined by JuMP that contains a reference to an optimization variable in a
   model. (Extra for experts: the `VariableRef` struct is a thin wrapper around
   a `MOI.VariableIndex`, and also contains a reference to the JuMP model.)

To illustrate these three types of variables, consider the following JuMP code
(the full syntax is explained below):
```jldoctest variables
julia> model = Model()
A JuMP Model

julia> @variable(model, x[1:2])
2-element Array{JuMP.VariableRef,1}:
 x[1]
 x[2]
```
This code does three things:
1. it adds two *optimization* variables to `model`
2. it creates two *JuMP* variables that act as references to those optimization
   variables
3. it binds those JuMP variables as a vector with two elements to the *Julia*
   variable `x`.

To reduce confusion, we will attempt, where possible, to always refer to
variables with their corresponding prefix.

JuMP variables can have attributes, such as names or an initial primal start
value. We illustrate the name attribute in the following example:
```jldoctest variables
julia> @variable(model, y, basename="decision variable")
decision variable
```
This code does four things:
1. it adds one *optimization* variable to `model`
2. it creates one *JuMP* variable that acts as a reference to that optimization
   variable
3. it binds the JuMP variable to the Julia variable `y`
4. it tells JuMP that the *name* attribute of this JuMP variable is "decision
variable". JuMP uses the value of `basename` when it has to print the variable
as a string.

For example, when we print `y` at the REPL we get:
```jldoctest variables
julia> y
decision variable
```

Because `y` is a JuMP variable, I can bind it to a different value. For example,
if I go:
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
# TODO(@odow): add a section on looking up by string
```

Now that we understand the different between *optimization*, *JuMP*, and *Julia*
variables, we can properly introduce the `@variable` macro.

The `@variable` macro
---------------------

We have already seen the basic usage of the `@variable` macro. The next
extension is to add lower and upper bounds to each optimization variable. This
can be done as follows:

```jldoctest variables
julia> @variable(model, a)
a

julia> @variable(model, b >= 0)
b

julia> @variable(model, c <= 1)
c

julia> @variable(model, 2 <= d <= 3)
d
```

We can query whether an optimization variable has a lower or upper bound via the `JuMP.haslowerbound` and `JuMP.hasupperbound` functions. For example:
```jldoctest variables
julia> JuMP.haslowerbound(a)
false

julia> JuMP.hasupperbound(c)
true
```

If a variable has a lower or upper bound, we can query the value of it via the `JuMP.lowerbound` and `JuMP.upperbound` functions. For example:
```jldoctest variables
julia> JuMP.lowerbound(d)
2.0

julia> JuMP.upperbound(d)
3.0
```
Querying the value of a bound that does not exist will result in an error.

Instead of using the `<=` and `>=` syntax, we can also use the `lowerbound` and `upperbound` keyword arguments. For example:
```jldoctest variables
julia> @variable(model, kw_var, lowerbound=1, upperbound=2)
kw_var

julia> JuMP.lowerbound(kw_var)
1.0
```

## Containers


### Arrays

```julia
julia> @variable(model, x[1:2, 1:2])
2x2 Array{JuMP.VariableRef,2}:
 x[1,1]  x[1,2]
 x[2,1]  x[2,2]
```

### JuMPArrays

```julia
julia> @variable(model, x[1:2, [:A,:B]])
2-dimensional JuMPArray{JuMP.Variableref,2,...} with index sets:
    Dimension 1, 1:2
    Dimension 2, Symbol[:A, :B]
And data, a 2x2 Array{JuMP.Variableref,2}:
 x[1,A]  x[1,B]
 x[2,A]  x[2,B]
```

### Dictionaries

```julia
julia> @variable(model, x[i=1:4; mod(i, 2)==0])
Dict{Any,JuMP.VariableRef} with 2 entries:
 4 => x[4]
 2 => x[2]
```

### Forcing the container type

```julia
julia> A = 1:2
1:2

julia> @variable(model, x[A])
1-dimensional JuMPArray{JuMP.Variableref,1,...} with index sets:
    Dimension 1, 1:2
And data, a 2-element Array{JuMP.VariableRef,1}:
 x[1]
 x[2]

julia> @variable(model, x[A], container=Array)
2-element Array{JuMP.VariableRef,1}:
 x[1]
 x[2]
```

## Variable types

`Bin`, `Int`, `Symmetric`, `PSD`

## Anonymous variables


```@docs
@variable
```

How to delete variables.
