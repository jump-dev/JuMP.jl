:mod:`JuMP` --- Julia for Mathematical Programming (JuMP)
=================================================================

.. module:: JuMP
   :synopsis: Julia for Mathematical Programming

JuMP is a domain-specific modeling language
for `mathematical programming <http://en.wikipedia.org/wiki/Mathematical_optimization>`_ 
embedded in `Julia <http://julialang.org/>`_. 
It it currently supports a number of open-source and commerical 
solvers (COIN Clp, COIN Cbc, GNU GLPK, and Gurobi) via a generic solver-independent
interface provided by the `MathProgBase <https://github.com/mlubin/MathProgBase.jl>`_ package. 
One the best features of JuMP is its speed -
benchmarking has shown that it can create problems at similar
speeds to special-purpose modeling languages such as AMPL
while maintaing the expressiveness of a generic high-level programming language.

.. _jump-installation:

------------
Installation
------------

JuMP can be installed through the Julia package manager (version 0.2 required)::

    julia> Pkg.add("JuMP")

This will install JuMP's default solvers: `Clp <https://github.com/mlubin/Clp.jl>`_ and `Cbc <https://github.com/mlubin/Cbc.jl>`_. All Julia-supported platforms (including Linux, OS X, and Windows) are supported by JuMP. JuMP itself is written in pure Julia; however, external solvers introduce binary dependencies. See the README for `Cbc <https://github.com/mlubin/Cbc.jl>`_ for platform-specifc installation details.  


-----------------
Quick Start Guide
-----------------

Defining Variables
^^^^^^^^^^^^^^^^^^

Variables, which are Julia objects,
are defined using the ``@defVar`` macro, where the first argument
will always be the ``Model``. In the examples below we assume ``m`` is already
defined. The second argument is an expression that declares the variable name
and optionally allows specification of lower and upper bounds. For example::

    @defVar(m, x )              # No bounds
    @defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @defVar(m, x <= ub )        # Upper bound only
    @defVar(m, lb <= x <= ub )  # Lower and upper bounds

All these variations introduce a new variable ``x`` in the local scope. For information about
common operations on variables, e.g. changing their bounds, see the
:ref:`variables <jump-variables>` section.

Integer and binary restrictions can optionally be specified with a third 
argument, ``Int`` or ``Bin``.

To create arrays of variables we append brackets to the variable name.
For example::

    @defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary
iterable sets are supported as index sets. Currently we only support ranges
of the form ``a:b`` where ``a`` is an explicit integer, not a variable. 
Using ranges will generally be faster than using arbitrary symbols. You can
mix both ranges and lists of symbols, as in the following example::

    s = ["Green","Blue"]
    @defVar(m, x[-10:10,s] , Int)
    x[-4,"Green"]

Finally, bounds can depend on variable indices::

@defVar(m, x[i=1:10] >= i )


Objective and Constraints
^^^^^^^^^^^^^^^^^^^^^^^^^

JuMP allows users to use a natural notation to describe linear expressions.
There are two ways to do so. The first is very similar to other modeling
languages and has no restrictions. The second utilizes Julia's powerful
metaprogramming features to get excellent performance even for large problems,
but has some restrictions on how they can be used.

To add constraints in the first way, use the ``addConstraint()`` and ``setObjective()``
functions, e.g.::

    setObjective(m, 5x + 22y + (x+y)/2)
    addConstraint(m, y + z == 4)  # Other options: <= and >=

The second way is very similar, and uses the ``@addConstraint`` and ``@setObjective``
macros, e.g.::

    @addConstraint(m, x[i] - s[i] <= 0)  
    @setObjective(m, sum{x[i], i=1:numLocation} )
    
There are some restrictions on what can go inside the expression:
 * If there is a product between coefficients and variables, the variables
   must appear last. That is, Coefficient times Variable is good, but 
   Variable times Coefficient is bad.
 * However, division by constants is supported.

You may have noticed a special ``sum{}`` operator above. The syntax is of the 
form::

	sum{expression, i = I1, j = I2, ...}

which is equivalent to::

    a = AffExpr()  # Create a new empty affine expression
    for i = I1
        for j = I2
            ...
            a += expression
            ...
        end
    end


You can also put a condition in::

    sum{expression, i = I1, j = I2, ...; cond} 

which is equivalent to::

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

Quadratic Objectives
^^^^^^^^^^^^^^^^^^^^

There is preliminary support for convex quadratic objectives. Currently the
only supported solver is ``Gurobi``; it must be set as the ``lpsolver`` or 
``mipsolver`` when solving QPs or mixed-integer QPs, respectively. The 
``@setObjective`` macro does not yet support quadratic terms, but you may
use instead the (slower) operator overloading functionality and the 
``setObjective`` function::

    MathProgBase.setlpsolver(:Gurobi)
    m = Model(:Min)
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    setObjective(m, x*x+ 2x*y + y*y )
    @addConstraint(m, x + y >= 1 )
      
    print(m)

    status = solve(m)

Quadratic Constraints
^^^^^^^^^^^^^^^^^^^^^

There is preliminary support for convex quadratic constraints. Currently the 
only supported solver is ``Gurobi``; it must be set as the ``lpsolver`` or
``mipsolver`` when solving QC programs. The ``@addConstraint`` macro does not 
yet support quadratic expressions, but you may instead use the (slower) 
operator overloading functionality via the ``addConstraint`` function::

    MathProgBase.setlpsolver(:Gurobi)
    m = Model(:Min)
    @defVar(m, -1 <= x <= 1)
    @defVar(m, -1 <= y <= 1)

    @setObjective(m, x + y)
    addConstraint(m, x*x + y*y <= 1)

    print(m)

    status = solve(m)


--------------
Simple Example
--------------

In this section we will construct a simple model and explain every step along the way. If you have used an algebraic modelling language before like AMPL, PULP or YALMIP then the following will most likely be familiar to you.

Heres the full piece of code, from the README file::

    using JuMP

    m = Model(:Max)
    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

    @setObjective(m, 5x + 3*y )
    @addConstraint(m, 1x + 5y <= 3.0 )
        
    print(m)
        
    status = solve(m)
        
    println("Objective value: ", getObjectiveValue(m))
    println("x = ", getValue(x))
    println("y = ", getValue(y))

Explanation
^^^^^^^^^^^

Once JuMP is :ref:`installed <jump-installation>`, to use JuMP in your programs, you just need to say::

    using JuMP

Models are created with the ``Model()`` function. This function takes one argument, the model sense. The two options are ``:Max`` and ``:Min``. Note: your model doesn't have to be called m - its just a variable name! Also, in case you were wondering, the colon ``:`` operator defines a `symbol <http://docs.julialang.org/en/latest/manual/metaprogramming/#symbols>`_::

    m = Model(:Max)

Defining variables is also easy. There are many options, depending on whether you want to have lower, upper, both, or no bounds. The names of your variables must be valid variable names. The following commands will create two variables, ``x`` and ``y``, with both lower and upper bounds. Note the first argument is our model variable 'm'. These variables are associated with this model, they cannot be used in another model.

::

    @defVar(m, 0 <= x <= 2 )
    @defVar(m, 0 <= y <= 30 )

Here we set our objective. Note again the m, so we know which model's objective we are setting! Note also that we don't have a multiplication ``*`` symbol between 5 and our variable name x. Julia is smart enough to not need it! Feel free to stick with ``*`` if it makes you feel more comfortable, as we have done with 3*y.

::

    @setObjective(m, 5x + 3*y )

Adding constraints is much like setting the objective. Here we create a less-than-or-equal-to constraint using <=, but we can also create equality constraints using == and greater-than-or-equal-to constraints with >=.

::

    @addConstraint(m, 1x + 5y <= 3.0 )

If you want to see what your model looks like in a human-readable format, the print function is defined for models. As you get to large models, this is probably not going to fit on your screen very well!

::

    print(m)

Models are solved with the solve() function. This function will never cause an error if your model is infeasible - instead it will return a flag. In this case, the model is feasible so the value of status will be ``:Optimal``, where ``:`` again denotes a symbol.

::

    status = solve(m)

Finally, we can access the results of our optimization. Getting the objective value is simple - its stored in the ``objVal`` field of the ``Model`` object.

::
    
    println("Objective value: ", getObjectiveValue(m))

To get the value from a variable, we call the ``getValue()`` function If ``x`` is not a single variable, but instead a range of variables, ``getValue()`` will return a list. In this case, however, it will just return a single value.

::
    
    println("x = ", getValue(x))
    println("y = ", getValue(y))

----------------------------------------------------
Expressions, constraints, and the objective function
----------------------------------------------------

Macros vs operator overloading @addConstraint duals @setObjective


.. _jump-variables:

---------
Variables
---------


Variables, also known as columns or decision variables, are the results of the optimization.

Creation
^^^^^^^^

The primary way to create variables is with the ``@defVar`` macro. The first argument will always be a ``Model``. In the examples below we assume ``m`` is already defined. The second argument is an expression that declares the variable name and optionally allows specification of lower and upper bounds. For example::

    @defVar(m, x )              # No bounds
    @defVar(m, x >= lb )        # Lower bound only (note: 'lb <= x' is not valid)
    @defVar(m, x <= ub )        # Upper bound only
    @defVar(m, lb <= x <= ub )  # Lower and upper bounds

All these variations create a new local variable, in this case ``x``. Integer and binary restrictions can optionally be specified with a third argument, ``Int`` or ``Bin``.

Arrays of variables
^^^^^^^^^^^^^^^^^^^

To create arrays of variables we append brackets to the variable name. For example::

    @defVar(m, x[1:M,1:N] >= 0 )

will create an ``M`` by ``N`` array of variables. Both ranges and arbitrary iterable sets are supported as index sets. Currently we only support ranges of the form ``a:b`` where ``a`` is an explicit integer, not a variable. Using ranges will generally be faster than using arbitrary symbols. You can mix both ranges and lists of symbols, as in the following example::

    s = ["Green","Blue"]
    @defVar(m, x[-10:10,s] , Int)
    x[-4,"Green"]

An interesting feature is that bounds can depend on variable indices::

    @defVar(m, x[i=1:10] >= i )

Variables can be constructed manually, one-by-one, using

::

    x = Variable(model, lower, upper, category)

but this is not considered idiomatic.

Modification
^^^^^^^^^^^^

Bounds
++++++
* ``setLower(x, lower)``, ``getLower(x)`` - Set/get the lower bound of a variable.
* ``setUpper(x, upper)``, ``getUpper(x)`` - Set/get the upper bound of a variable.

Values at solution
++++++++++++++++++
* ``getValue(x)`` - Get the value of this variable in the solution.
* ``getDual(x)`` - Get the reduced cost of this variable in the solution.

Names
+++++
Variables can have internal names that can be used for writing models to file. This is currently not done for performance reasons, but may be added if there is demand.
* ``setName(x, newName)``, ``getName(x)`` - Set/get the variable's internal name.



.. _jump-function-list:

-------------------------------
List of functions in JuMP module
-------------------------------

.. function:: setName(x, newName)

    Assigns a name to the the variable ``x``.


