.. _extensions:

--------------
Extending JuMP
--------------

.. note::
    This section documents internal JuMP APIs which may change between versions without warning.

In addition to the problem classes directly supported by JuMP, it is possible write
your own extensions to JuMP to handle problem structures which don't fall nicely into
the standard linear-programming-type modeling framework. Extensions are particularly
powerful for making advanced algorithmic techniques for specialized problem classes
accessible to users who just want to write down and solve their problem.
If only more papers came with a JuMP extension!

[List existing extensions?]

Transformations
^^^^^^^^^^^^^^^

The simplest way to write a JuMP extension is by providing a Julia function that
performs an immediate transformation. For example, the function::

    function abs(x::JuMP.Variable)
        model = x.m # access internal field of a JuMP Variable
        @defVar(model, z >= 0)
        @addConstraint(model, z >= x)
        @addConstraint(model, z >= -x)
        return z
    end

returns a variable ``z`` which satisfies :math:`z \geq |x|`, a
standard transformation in linear programming.
Users could then use ``abs`` within their JuMP model, e.g.::

    @setObjective(m, Min, sum{ abs(x[i]), i in 1:N })

to represent minimization of an :math:`\ell_1` norm as a linear objective.

For a slightly more complex example, consider the constraint that a vector of variables
``x`` belongs to a polyhedron given by the the convex hull of a fixed set of vectors.
That is, :math:`x \in \operatorname{conv}\{z_1,z_2,\ldots,z_k\}`. This holds iff
:math:`\exists\, \lambda \geq 0 \in \mathbb{R}^k \text{ with } \sum_{i=1}^k \lambda_i = 1  \text{ such that } x = \sum_{i=1}^k \lambda_i z_i`.
This constraint can be implemented as::

    function inconvexhull(m::Model, x, z::Vector{Vector{Float64}})
        k = length(z)
        n = length(x)
        @defVar(m, λ[1:length(z)] >= 0)
        @addConstraint(m, sum(λ) == 1)
        @addConstraint(m, cvxhull[i=1:k,j=1:n], x[j] == sum{ z[i][j]*λ[i], i = 1:k } )
    end


A library of such transformations could well be considered a JuMP extension.
These types of extensions require very little access to JuMP internals,
although they allow little flexibility in syntax; the user interacts with
the extension by calling a plain Julia function with specified arguments.
At the same time, one could write macros that
invoke these functions with special syntax.
For example in the above case, a macro could
translate a statement like ``@inconv m x {z1,z2,z3}`` to
``inconvexhull(m, x, [z1,z2,z3])``, again without needing to access JuMP internals.


Specialized solution methods for structured problems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, transformations to a standard form recognized by
JuMP are either impossible or not computationally advantageous.
In this case, it's important to store some extra structure
of the model as it's being constructed and then apply
a specialized algorithm to solve the optimization problem.
These types of extensions require deeper knowledge of JuMP's
internals and macro metaprogramming.

[ Walk through developing an extension for the convex hull constraint above? ]


