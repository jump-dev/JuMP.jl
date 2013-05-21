MathProg.jl
===========

Linear Programming, Quadratic Programming and Integer Programming 
modeling with Julia. The goal is to have the speed of AMPL embedded in
a fully functional language. Compare with PuLP/Pyomo/YALMIP/...

This package is not related to GNU MathProg.

# Installation

You can install MathProg through the Julia package manager (version 0.2 prerelease required):

    Pkg.add("MathProg")
    
MathProg is built on **[MathProgBase]** which provides solver-independent
LP and MILP functionality. See that package for how to select solvers.
By default we install and use the **[Clp]** and **[CoinMP]** packages which 
link with powerful open-source LP and MILP solvers. Check the corresponding
READMEs for platform-specific installation instructions.

[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
[Clp]: https://github.com/mlubin/Clp.jl
[CoinMP]: https://github.com/mlubin/CoinMP.jl

# Simple Example

    using MathProg

    m = Model("max")
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    @setObjective(m, 5x + 3y )
    @addConstraint(m, 1x + 5y <= 3.0 )
    
    print(m)
    
    status = solve(m)
    
    println("Objective value: ", m.objVal)
    println("x = ", getValue(x))
    println("y = ", getValue(y))

# Defining variables

Variables are defined using the ``@defVar`` macro. The first argument must be a ``Model``.
In the examples below we assume ``m`` is already defined.
The second is an expression that declares the variable name and optionally allows specification
of lower and upper bounds. The possible combinations are

No bounds:

	@defVar(m, x )

Lower bound only:

	@defVar(m, x >= lb )

Upper bound only:

	@defVar(m, x <= ub )

Lower and upper bounds:

	@defVar(m, lb <= x <= ub )

Where ``x`` must be a valid symbol which will be assigned to in the local context.

Integer and binary restrictions can optionally be specified with a third argument, ``Int`` or ``Bin``.

Arrays of variable objects can created by appending brackets to the variable name.
For example,

	@defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form ``a:b`` where ``a`` is an explicit integer, not a variable. Using ranges will generally be faster, as this avoids dictionary lookups. The following code

	s = ["Green","Blue"]
	@defVar(m, x[-10:10,s] , Int)
	x[-4,"Green"]

works fine.

Bounds can depend on variable indices:

	@defVar(m, x[i=1:10] >= i )

The macro is just a wrapper for the ``Variable`` constuctor, which can be called directly by advanced users.

    
# Defining linear expressions

Using Julia's powerful metaprogramming features, we can turn easy-to-read
statements into a sparse internal representation very quickly. To 
invoke this, use the ``@addConstraint`` (and ``@setObjective``) macros. Here are some examples:

    @addConstraint(m, x[i] - s[i] <= 0)
    
	@setObjective(m, sum{x[i], i=1:numLocation} )
    
There are some restrictions on what can go inside the macros
 * If there is a product between coefficients and variables, the variables
   must appear last. That is, Coefficient times Variable is good, but 
   Variable times Coefficient is bad.
 * Addition, subtraction, distributive rule should all mostly work, 
   but the parser has not been well tested. Please report any issues.

The ``sum`` "operator" behaves as follows.
	
	sum{expression, i = I1, j = I2, ...}

is equivalent to

	x = AffExpr()
	for i = I1
		for j = I2
			...
				x += expression
			...
		end
	end


We also allow conditions:

	sum{expression, i = I1, j = I2, ...; cond} 

is equivalent to

	x = AffExpr()
	for i = I1
		for j = I2
			...
				if cond
					x += expression
				end
			...
		end
	end


# Function listing

`Model(sense)` 
 * Construct a Model with the objective sense provided. Use either "max" or "min"
 

`print(model)`
 * Displays the Model to the screen in a user-readable way


`Variable(model, lower, upper, category[, name])`
 * Create a new variable in the model. `lower` and `upper` are the bounds on,
   variable, and `category` should be one of CONTINUOUS, INTEGER,
   or BINARY. `name` is an optional string argument. The ``@defVar`` macro should
   typically be used instead of this.


`setName(v,name)`,`getName(v)`

`setLower(v,lower)`,`getLower(v)`

`setUpper(v,upper)`,`getUpper(v)`

 * Update variable options.

`getValue(v)`

 * Retrieve the optimal value of a variable after a solution. May also be
   called with variable collections, in which cases it returns a collection of
   solution values.

`setObjective(model, affexpr)`
 * Sets the objective to `affexpr`, where `affexpr` is any valid combination of 
   variables and constants ("affine expression"). This is the slow approach based
   on operator overloading. Use ``@setObjective`` instead.

`addConstraint(model, constraint)`
 * Adds a constraint, where `constraint` is of the form `affexpr <= number`, 
   `affexpr == number`, or `affexpr >= number`. This is the slow approach based on
   operator overloading. Use ``@addConstraint`` instead.

`writeLP(model, filename)`
 * Writes the model out to LP format, which should be readable by most solvers.

`writeMPS(model, filename)`
 * Writes the model out to MPS format, which should be readable by most solvers.

`solve(model)`
 * Solves the model. If any integer variables a present, it will use CoinMP, 
   otherwise Clp is used. It returns a flag representing the solver status.
   It also sets the objVal and colVal fields of the model.

# Quadratic objective support

There is partial support for quadratic objectives. You can formulate problems with them, 
and they can be written out to MPS files. More support to follow.
