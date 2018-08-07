Variables
=========

What is a JuMP variable?
---------------------------

The term _variable_ in computational optimization has many meanings. Here, we
distinguish between the following three types of variables:
1. _optimization_ variables, which are the mathematical ``x`` in the problem
``\\max{c^\top x | Ax = b\}``.
2. _Julia_ variables, which are bindings between a name and a value, for example
`x = 1`. (See
[here](https://docs.julialang.org/en/stable/manual/variables/) for the Julia
docs.)
3. _JuMP_ variables, which are instances of the `JuMP.VariableRef` struct
defined by JuMP that contains a reference to an optimization variable in a
model. (Extra for experts: the `VariableRef` struct is a thin wrapper around a
`MOI.VariableIndex`, and also contains a reference to the JuMP model.)

To illustrate these three types of variables, consider the following JuMP code
(the full syntax is explained below):
```julia
model = Model()
@variable(model, x[1:2])
```
This code does three things:
1. it tells `model` to add two _optimization_ variables
2. it creates two _JuMP_ variables that act as references to those
optimization variables
3. it binds those JuMP variables as a vector with two elements to the name `x`.

Therefore, we have:
```julia
julia> x
2-element Array{JuMP.VariableRef, 1}
 x[1]
 x[2]
```

To reduce confusion, we will attempt, where possible, to always refer to
variables with their corresponding prefix.

JuMP variables can have attributes, such as names or an initial primal start
value. We illustrate the name attribute in the following example:
```julia
model = Model()
@variable(model, x, basename="decision variable")
```
This code does four things:
1. it tells `model` to add one _optimization_ variable
2. it creates one _JuMP_ variable that acts as a reference to that
optimization variables
3. it binds the JuMP variable to the name `x`
4. it tells JuMP that the _name_ attribute of this JuMP variable is "decision
variable". JuMP uses the value of `basename` when it has to print the variable
as a string.

For example:
```julia
julia> x
decision variable

julia> typeof(x)
JuMP.VariableRef
```

Because `x` is a JuMP variable, I can bind it to a different value. For example,
if I go:
```julia
julia> x = 1
1
```
`x` is no longer a binding to a JuMP variable. This does not mean that the JuMP
variable has been destroyed. It still exists and is still a reference to the
same optimization variable. The binding can be reset by querying the model for
the symbol _as it was written in the `@variable` macro_. For example:
```julia
julia> model[:x]
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

The `@variable` macro
---------------------

DRAFT: Describe the complete syntax of the `@variable` macro. Anonymous versus
named variables. Describe the three possible container types returned and how
to use them (`Array`, `JuMPArray`, and `Dict`).

```@docs
@variable
```

How to delete variables.
