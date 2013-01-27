Jump
====

Linear Programming, Quadratic Programming and Integer Programming 
modelling with Julia. The goal is to have the speed of AMPL embedded in
a fully functional language. Compare with PuLP/Pyomo/YALMIP/...

    Julia + Mathematical Programming = JuMP

# Installation

You can install Jump through the Julia package manager:
    Pkg.add("Jump")
    
Jump is currently dependent on the CLP package which provides a powerful
open-source LP solver. As the infrastructure and interfaces for solvers
are developed in Julia we will extend the functionality. Check the Clp.jl
README to see if your platform is supported.

# Simple Example

    using Jump

    m = Model(JUMP_MAX)

    x = addVar(m, 0,  2, JUMP_CONTINUOUS)
    y = addVar(m, 0, 30, JUMP_CONTINUOUS, "y")

    setObjective(m, 5x + 3y)
    addConstraint(m, 1x + 5y <= 3.0)
    
    print(m)
    
    status = solveClp(m)
    
    println("Objective value: ", m.objVal)
    println("x = ", getValue(x))
    println("y = ", getValue(y))
    
# Getting Speed

Using Julia's powerful metaprogramming features, we can turn easy-to-read
statements into the a sparse internal representation very quickly. To 
invoke this, use the @addConstraint (and @setObjective) macros. Here are some examples:

    @addConstraint(m, x[i] - s[i] <= 0)
    
	@setObjective(m, sum{x[i], i=1:numLocation} )
    
There are some restrictions on what can go inside the macros
 * If there is a product between coefficients and variables, the variables
   must appear last. That is, Coefficient times Variable is good, but 
   Variable times Coefficient is bad.
 * Addition, subtraction, distributive rule should all mostly work, 
   but the parser has not been well tested. Please report any issues.

The ``sum`` expression behaves as follows.
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


# Full function listing

`Model(sense)` 
 * Construct a Model with the objective sense provided. There are two
   constants exported by Jump, JUMP_MAXIMIZE and JUMP_MINIMIZE, that
   you should use.
 

`print(model)`
 * Displays the Model to the screen in a user-readable way


`addVar(model, lower, upper, category[, name])`
 * Add a variable to the model. `lower` and `upper` are the bounds on,
   variable, and `category` should be one of JUMP_CONTINUOUS, JUMP_INTEGER,
   or JUMP_BINARY. `name` is an optional string argument.


`addVars(model, lower, upper, category, dims[, name])`
 * Add multiple variables, returns them as a list of variables. Same as
   above, except for dims. If dims is a single integer, it indexes the variables
   from 1 to dims, e.g. `x[3]`. If dims is a tuple of two integers, e.g. (10,5) 
   it indexes them along the two dimensions, e.g. `x[1,1]` to `x[10,5]`

`setName(v,name)`,`getName(v)`

`setLower(v,lower)`,`getLower(v)`

`setUpper(v,upper)`,`getUpper(v)`

`getValue(v)`
 * Update variable options, and retrieve the value of the variable after solution.

`setObjective(model, affexpr)`
 * Sets the objective to `affexpr`, where `affexpr` is any valid combination of 
   variables and constants ("affine expression")

`addConstraint(model, constraint)`
 * Adds a constraint, where `constraint` is of the form `affexpr <= number`, 
   `affexpr == number`, or `affexpr >= number`. This is very limiting form
   but will be extended soon.

`writeLP(model, filename)`
 * Writes the model out to LP format, which should be readable by most solvers.

`writeMPS(model, filename)`
 * Writes the model out to MPS format, which should be readable by most solvers.

# Quadratic objective support

There is partial support for quadratic objectives. You can formulate problems with them, 
and they can be written out to MPS files. More support to follow.
