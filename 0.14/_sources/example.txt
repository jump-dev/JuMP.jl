.. _simple-example:

Simple Example
^^^^^^^^^^^^^^

In this section we will construct a simple model and explain every step along the way.
The are more complex examples in the ``JuMP/examples/`` `folder <https://github.com/JuliaOpt/JuMP.jl/tree/master/examples>`_. Here is the code we will walk through::

    using JuMP

    m = Model()
    @variable(m, 0 <= x <= 2 )
    @variable(m, 0 <= y <= 30 )

    @objective(m, Max, 5x + 3*y )
    @constraint(m, 1x + 5y <= 3.0 )

    print(m)

    status = solve(m)

    println("Objective value: ", getobjectivevalue(m))
    println("x = ", getvalue(x))
    println("y = ", getvalue(y))

Once JuMP is :ref:`installed <jump-installation>`, to use JuMP in your
programs, you just need to say::

    using JuMP

Models are created with the ``Model()`` function::

    m = Model()

.. note::
   Your model doesn't have to be called m - it's just a name.

There are a few options for defining a variable, depending on whether you want
to have lower bounds, upper bounds, both bounds, or even no bounds. The following
commands will create two variables, ``x`` and ``y``, with both lower and upper bounds.
Note the first argument is our model variable ``m``. These variables are associated
with this model and cannot be used in another model.::

    @variable(m, 0 <= x <= 2 )
    @variable(m, 0 <= y <= 30 )

Next we'll set our objective. Note again the ``m``, so we know which model's
objective we are setting! The objective sense, ``Max`` or ``Min``, should
be provided as the second argument. Note also that we don't have a multiplication ``*``
symbol between 5 and our variable ``x`` - Julia is smart enough to not need it!
Feel free to stick with ``*`` if it makes you feel more comfortable, as we have
done with ``3*y``::

    @objective(m, Max, 5x + 3*y )

Adding constraints is a lot like setting the objective. Here we create a
less-than-or-equal-to constraint using ``<=``, but we can also create equality
constraints using ``==`` and greater-than-or-equal-to constraints with ``>=``::

    @constraint(m, 1x + 5y <= 3.0 )

If you want to see what your model looks like in a human-readable format,
the ``print`` function is defined for models.

::

    print(m)

Models are solved with the ``solve()`` function. This function will not raise
an error if your model is infeasible - instead it will return a flag. In this
case, the model is feasible so the value of ``status`` will be ``:Optimal``,
where ``:`` again denotes a symbol. The possible values of ``status``
are described :ref:`here <solvestatus>`.

::

    status = solve(m)

Finally, we can access the results of our optimization. Getting the objective
value is simple::

    println("Objective value: ", getobjectivevalue(m))

To get the value from a variable, we call the ``getvalue()`` function. If ``x``
is not a single variable, but instead a range of variables, ``getvalue()`` will
return a list. In this case, however, it will just return a single value.

::

    println("x = ", getvalue(x))
    println("y = ", getvalue(y))
