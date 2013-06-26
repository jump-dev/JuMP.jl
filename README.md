MathProg.jl
===========

MathProg.jl is a linear/integer/quadratic programming modelling language
embedded in the Julia language. It is currently connected to the COIN CLP and
CBC solvers, GNU GLPK, and the Gurobi solver. One the best features of MathProg
is its speed - benchmarking has shown that it can create problems at similar
speeds to special-purpose modelling languages such as AMPL.

This package is not related to GNU MathProg.

[![Build Status](https://travis-ci.org/IainNZ/MathProg.jl.png)](https://travis-ci.org/IainNZ/MathProg.jl)

# Installation

You can install MathProg through the Julia package manager (version 0.2 prerelease required):

```jl
Pkg.add("MathProg")
```

This will install MathProg's dependencies **[MathProgBase]** (which provides
a generic solver interface) and the **[Clp]** and **[CoinMP]** packages which
link with powerful open-source LP and MILP solvers, which use by default.

For more installation information, see the corresponding READMEs for more instructions.

[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
[Clp]: https://github.com/mlubin/Clp.jl
[CoinMP]: https://github.com/mlubin/CoinMP.jl

# Simple Example

```jl
using MathProg

m = Model(:Max)
@defVar(m, 0 <= x <= 2 )
@defVar(m, 0 <= y <= 30 )

@setObjective(m, 5x + 3y )
@addConstraint(m, 1x + 5y <= 3.0 )
    
print(m)
    
status = solve(m)
    
println("Objective value: ", m.objVal)
println("x = ", getValue(x))
println("y = ", getValue(y))
```

See the [SimpleExample] wiki page for a detailed explanation, and for 
more examples, see the [examples/] folder.

[SimpleExample]: https://github.com/IainNZ/MathProg.jl/wiki/Simple-Example
[examples/]: https://github.com/IainNZ/MathProg.jl/tree/master/examples

# Quick Start Guide

For a full user guide, please see the [wiki](https://github.com/IainNZ/MathProg.jl/wiki).

### Defining Variables

Variables are defined using the ``@defVar`` macro, where the first argument
will always be the ``Model``. In the examples below we assume ``m`` is already
defined. The second argument is an expression that declares the variable name
and optionally allows specification of lower and upper bounds. For example:

```jl
@defVar(m, x )              # No bounds
@defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
@defVar(m, x <= ub )        # Upper bound only
@defVar(m, lb <= x <= ub )  # Lower and upper bounds
```

All these variations create a new local variable ``x``. For information about
common operations on variables, e.g. changing their bounds, see the
[wiki](https://github.com/IainNZ/MathProg.jl/wiki/Variables).

Integer and binary restrictions can optionally be specified with a third 
argument, ``Int`` or ``Bin``.

To create an arrays of variables we appending brackets to the variable name.
For example

```jl
@defVar(m, x[1:M,1:N] >= 0 )
```

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary
iterable sets are supported as index sets. Currently we only support ranges
of the form ``a:b`` where ``a`` is an explicit integer, not a variable. 
Using ranges will generally be faster than using arbitrary symbols. You can
mix both ranges and lists of symbols, as in the following example:

```jl
s = ["Green","Blue"]
@defVar(m, x[-10:10,s] , Int)
x[-4,"Green"]
```

Finally, bounds can depend on variable indices:
```jl
@defVar(m, x[i=1:10] >= i )
```

### Objective and Constraints

MathProg allows users to use a natural notation to describe linear expressions.
There are two ways to do so. The first is very similar to other modelling
languages and has no restrictions. The second utilizes Julia's powerful
metaprogramming features to get excellent performance even for large problems,
but has some restrictions on how they can be used.

To add constraints in the first way, use the ``addConstraint()`` and ``setObjective()``
functions, e.g.

```jl
setObjective(m, 5x + 22y + (x+y)/2)
addConstraint(m, y + z == 4)  # Other options: <= and >=
```

The second way is very similar, and uses the ``@addConstraint`` and ``@setObjective``
macros, e.g.

```jl
@addConstraint(m, x[i] - s[i] <= 0)  
@setObjective(m, sum{x[i], i=1:numLocation} )
```
    
There are some restrictions on what can go inside the expression:
 * If there is a product between coefficients and variables, the variables
   must appear last. That is, Coefficient times Variable is good, but 
   Variable times Coefficient is bad.
   * However, division by constants is supported.

You may have noticed a special ``sum{}`` operator above. The syntax is of the 
form
```jl	
	sum{expression, i = I1, j = I2, ...}
```
which is equivalent to
```jl
a = AffExpr()  # Create a new empty affine expression
for i = I1
    for j = I2
        ...
        a += expression
        ...
    end
end
```

You can also put a condition in:
```jl
sum{expression, i = I1, j = I2, ...; cond} 
```
which is equivalent to

```jl
a = AffExpr()
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

### Quadratic Objectives ###

There is preliminary support for convex quadratic objectives. Currently the
only supported solver is ``Gurobi``; it must be set as the ``lpsolver`` or 
``mipsolver`` when solving QPs or mixed-integer QPs, respectively. The 
``@setObjective`` macro does not yet support quadratic terms, but you may
use instead the (slower) operator overloading functionality and the 
``setObjective`` function:

```jl
MathProgBase.setlpsolver(:Gurobi)
m = Model(:Min)
@defVar(m, 0 <= x <= 2 )
@defVar(m, 0 <= y <= 30 )

setObjective(m, x*x+ 2x*y + y*y )
@addConstraint(m, x + y >= 1 )
  
print(m)

status = solve(m)
```
