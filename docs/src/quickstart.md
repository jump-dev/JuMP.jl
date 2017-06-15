Quick Start Guide
=================

This quick start guide will introduce the main concepts of JuMP. If you are familiar with another modeling language embedded in a high-level language such as PuLP (Python) or a solver-specific interface you will find most of this familiar, with the exception of *macros*. A deep understanding of macros is not essential, but if you would like to know more please see the [Julia documentation](http://docs.julialang.org/en/latest/manual/metaprogramming/). If you are coming from an AMPL or similar background, you may find some of the concepts novel but the general appearance will still be familiar.

Creating a Model
----------------

**Models** are Julia objects. They are created by calling the constructor:

```julia
m = Model()
```

All variables and constraints are associated with a `Model` object. Usually, you'll also want to provide a solver object here by using the `solver=` keyword argument; see the simple example below. For a list of all functions related to `Model`, see [Models](@ref).

Defining Variables
------------------

**Variables** are also Julia objects, and are defined using the `@variable` macro. The first argument will always be the `Model` to associate this variable with. In the examples below we assume `m` is already defined. The second argument is an expression that declares the variable name and optionally allows specification of lower and upper bounds. For example:

```julia
@variable(m, x )              # No bounds
@variable(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
@variable(m, x <= ub )        # Upper bound only
@variable(m, lb <= x <= ub )  # Lower and upper bounds
```

All these variations introduce a new variable `x` in the local scope. The names of your variables must be valid Julia variable names. For information about common operations on variables, e.g. changing their bounds, see the [Variables](@ref) section.

**Integer** and **binary** restrictions can optionally be specified with a third argument, `Int` or `Bin`.

To create arrays of variables we append brackets to the variable name. For example:

```julia
@variable(m, x[1:M,1:N] >= 0 )
```

will create an `M` by `N` array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form `a:b` where `a` is an explicit integer, not a variable. Using ranges will generally be faster than using arbitrary symbols. You can mix both ranges and lists of symbols, as in the following example:

```julia
s = ["Green", "Blue"]
@variable(m, x[-10:10,s], Int )
# e.g. x[-4, "Green"]
```

Finally, bounds can depend on variable indices:

```julia
@variable(m, x[i=1:10] >= i )
```

Objective and Constraints
-------------------------

JuMP allows users to use a natural notation to describe linear expressions. To add constraints, use the `@constraint()` and `@objective()` macros, e.g.:

```julia
@constraint(m, x[i] - s[i] <= 0)  # Other options: == and >=
@constraint(m, sum(x[i] for i=1:numLocation) == 1)
@objective(m, Max, 5x + 22y + (x+y)/2) # or Min
```

!!! note
    The `sense` passed to `@objective` must be a [symbol](http://docs.julialang.org/en/latest/manual/metaprogramming/#symbols) type: `:Min` or `:Max`, although the macro accepts `:Min` and `:Max`, as well as `Min` and `Max` (without the colon) directly.

The `sum()` syntax directly follows Julia's own generator expression syntax. You may use conditions within sums, e.g.:

```julia
sum(expression for i = I1, j = I2 if cond)
```

which is equivalent to:

```julia
a = zero(AffExpr)
for i = I1
    for j = I2
        ...
        if cond
            a += expression
        end
        ...
    end
end
```

!!! note
    JuMP previously used a special curly brace syntax for `sum{}`, `prod{}`, and `norm2{}`. This has been entirely replaced by `sum()`, `prod()`, and `norm()` since Julia 0.5. The curly brace syntax is deprecated and will be removed in a future release.


Simple Example
==============

In this section we will construct a simple model and explain every step along the way. The are more complex examples in the `JuMP/examples/` [folder](https://github.com/JuliaOpt/JuMP.jl/tree/master/examples). Here is the code we will walk through:

```julia
using JuMP
using Clp

m = Model(solver = ClpSolver())
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )

@objective(m, Max, 5x + 3*y )
@constraint(m, 1x + 5y <= 3.0 )

print(m)

status = solve(m)

println("Objective value: ", getobjectivevalue(m))
println("x = ", getvalue(x))
println("y = ", getvalue(y))
```

Once JuMP is installed [Installation Guide](@ref), to use JuMP in your programs, you just need to say:

```julia
using JuMP
```

We also need to include a Julia package which provides an appropriate solver. In this case, we'll use Clp:

```julia
using Clp
```

Models are created with the `Model()` function. The `solver=` keyword argument is used to specify the solver to be used:

```julia
m = Model(solver = ClpSolver())
```

!!! note
    Your model doesn't have to be called `m` - it's just a name.

There are a few options for defining a variable, depending on whether you want to have lower bounds, upper bounds, both bounds, or even no bounds. The following commands will create two variables, `x` and `y`, with both lower and upper bounds. Note the first argument is our model variable `m`. These variables are associated with this model and cannot be used in another model.:

```julia
@variable(m, 0 <= x <= 2 )
@variable(m, 0 <= y <= 30 )
```

Next we'll set our objective. Note again the `m`, so we know which model's objective we are setting! The objective sense, `Max` or `Min`, should be provided as the second argument. Note also that we don't have a multiplication `*` symbol between 5 and our variable `x` - Julia is smart enough to not need it! Feel free to stick with `*` if it makes you feel more comfortable, as we have done with `3*y`:

```julia
@objective(m, Max, 5x + 3*y )
```

Adding constraints is a lot like setting the objective. Here we create a less-than-or-equal-to constraint using `<=`, but we can also create equality constraints using `==` and greater-than-or-equal-to constraints with `>=`:

```julia
@constraint(m, 1x + 5y <= 3.0 )
```

If you want to see what your model looks like in a human-readable format, the `print` function is defined for models.

```julia
print(m)
```

Models are solved with the `solve()` function. This function will not raise an error if your model is infeasible - instead it will return a flag. In this case, the model is feasible so the value of `status` will be `:Optimal`, where `:` again denotes a symbol. The possible values of `status` are described in [Solve Status](@ref).

```julia
status = solve(m)
```

Finally, we can access the results of our optimization. Getting the objective value is simple:

```julia
println("Objective value: ", getobjectivevalue(m))
```

To get the value from a variable, we call the `getvalue()` function. If `x` is not a single variable, but instead a range of variables, `getvalue()` will return a list. In this case, however, it will just return a single value.

```julia
println("x = ", getvalue(x))
println("y = ", getvalue(y))
```
