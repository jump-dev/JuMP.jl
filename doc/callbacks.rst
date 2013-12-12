.. _callbacks:

----------------
Solver Callbacks
----------------

Many solvers offer the ability to modify the solve process. Examples include
changing branching decisions in branch-and-bound, adding custom cuts, providing
solvers with integer solutions, or adding new constraints only when they are
violated by the current solution (lazy constraints).

Solver-independent modelling languages do not, in general, provide a way to
provide callbacks that will work with any solver. However, JuMP does provide
limited support for this functionality. Currently we have cross-solver support
for adding "lazy constraints" for the Gurobi and CPLEX solvers.

Lazy Constraints
^^^^^^^^^^^^^^^^

Lazy constraints are useful when the full set of constraints is too large to
explicitly include in the initial formulation. When a MIP solver reaches a new
solution, for example with a heuristic or by solving a problem at a node in the
branch-and-bound tree, it will give the user the chance to provide constraint(s)
that would make the current solution infeasible. For some more information about
lazy constraints, see this blog post by `Paul Rubin <http://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html>`_.

There are three important steps to providing a lazy constraint callback. First we
must write a function that will analyze the current solution that takes a 
single argument, e.g. ``function myLazyCutGenerator(cb)``, where cb is a reference
to the callback management code inside JuMP. Next you will do whatever
analysis of the solution you need to inside your function to generate the new
constraint before adding it to the model with the JuMP function
``addLazyConstraint(cb, myconstraint)`` or the macro version
``@addLazyConstraint(cb, myconstraint)`` (same limitations as addConstraint).
Finally we notify JuMP that this function should be used for lazy constraint
generation using the ``setlazycallback(m, myLazyCutGenerator)`` function 
before we call ``solve(m)``.

The following is a simple example to make this more clear. In this two-dimensional
problem we have a set of box constraints explicitly provided and a set of two
lazy constraints we can add on the fly. The solution without the lazy constraints
will be either (0,2) or (2,2), and the final solution will be (1,2)::

    using JuMP
    using Gurobi

    # We will use Gurobi, which requires that we manually set the attribute
    # LazyConstraints to 1 if we use lazy constraint generation
    m = Model(solver=GurobiSolver(LazyConstraints=1))

    # Define our variables to be inside a box, and integer
    @defVar(m, 0 <= x <= 2, Int)
    @defVar(m, 0 <= y <= 2, Int)

    @setObjective(m, Max, y)

    # We now define our callback function that takes one argument,
    # the callback handle. Note that we can access m, x, and y because
    # this function is defined inside the same scope
    function corners(cb)
        x_val = getValue(x)
        y_val = getValue(y)
        println("In callback function, x=$x_val, y=$y_val")

        # We have two constraints, one cutting off the top
        # left corner and one cutting off the top right corner, e.g.
        # (0,2) +---+---+ (2,2)
        #       |xx/ \xx|
        #       |x/   \x|
        #       |/     \|
        #       +       +
        #       |       |
        #       |       |
        #       |       |
        # (0,0) +---+---+ (2,0)

        # Allow for some impreciseness in the solution
        TOL = 1e-6

        # Check top left, allowing some tolerance
        if y_val - x_val > 1 + TOL
            # Cut off this solution
            println("Solution was in top left, cut it off")
            # Use the original variables, but not m - cb instead
            @addLazyConstraint(cb, y - x <= 1)
        # Check top right
        elseif y_val + x_val > 3 + TOL
            # Cut off this solution
            println("Solution was in top right, cut it off")
            # Use the original variables, but not m - cb instead
            @addLazyConstraint(cb, y + x <= 3)
        end
    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    setlazycallback(m, corners)

    # Solve the problem
    solve(m)

    # Print our final solution
    println("Final solution: [ $(getValue(x)), $(getValue(y)) ]")

The code should print something like (amongst the output from Gurobi)::
    
    In callback function, x=2.0, y=2.0
    Solution was in top right, cut it off
    In callback function, x=0.0, y=2.0
    Solution was in top left, cut it off
    In callback function, x=1.0, y=2.0
    Final solution: [ 1.0, 2.0 ]

This code can also be found in ``/JuMP/examples/simplelazy.jl``.

Code Design Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above example the callback function is defined in the same scope as
the model and variable definitions, allowing us to access them. If we defined
the function in some other scope, or even file, we would not be able to access them directly.
The proposed solution to this design problem is to seperate the logic of analyzing the
current solution values from the callback itself. This has many benefits,
including writing unit tests for the callback function to check its
correctness. The callback function pased to JuMP is then simply a stub
that extracts the current solution and any other relevant information
and passes that to the constraint generation logic. To apply this to our
previous example, consider the following code::

    using JuMP
    using Gurobi
    using Base.Test

    function cornerChecker(x_val, y_val)
        # This function does not depend on the model, and could
        # be written anywhere. Instead, it returns a tuple of
        # values (newcut, x_coeff, y_coeff, rhs) where newcut is a
        # boolean if a cut is needed, x_coeff is the coefficient
        # on the x variable, y_coeff is the coefficient on
        # the y variable, and rhs is the right hand side
        TOL = 1e-6
        if y_val - x_val > 1 + TOL
            return (true, -1.0, 1.0, 1.0)  # Top left
        elseif y_val + x_val > 3 + TOL
            return (true,  1.0, 1.0, 3.0)  # Top right
        else
            return (false, 0.0, 0.0, 0.0)  # No cut
        end
    end

    # A unit test for the cornerChecker function
    function test_cornerChecker()
        # Test the four corners - only two should produce cuts

        newcut, x_coeff, y_coeff, rhs = cornerChecker(0, 0)
        @test !newcut

        newcut, x_coeff, y_coeff, rhs = cornerChecker(2, 0)
        @test !newcut

        newcut, x_coeff, y_coeff, rhs = cornerChecker(0, 2)
        @test newcut
        @test x_coeff == -1.0
        @test y_coeff ==  1.0
        @test rhs == 1.0

        newcut, x_coeff, y_coeff, rhs = cornerChecker(2, 2)
        @test newcut
        @test x_coeff ==  1.0
        @test y_coeff ==  1.0
        @test rhs == 3.0
    end

    function solveProblem()
        m = Model(solver=GurobiSolver(LazyConstraints=1))

        @defVar(m, 0 <= x <= 2, Int)
        @defVar(m, 0 <= y <= 2, Int)
        @setObjective(m, Max, y)

        # Note that the callback is now a stub that passes off
        # the work to the "algorithm"
        function corners(cb)
            x_val = getValue(x)
            y_val = getValue(y)
            println("In callback function, x=$x_val, y=$y_val")
            
            newcut, x_coeff, y_coeff, rhs = cornerChecker(x_val, y_val)

            if newcut
                @addLazyConstraint(cb, x_coeff*x + y_coeff*y <= rhs)
            end
        end  # End of callback function

        setlazycallback(m, corners)
        solve(m)
        println("Final solution: [ $(getValue(x)), $(getValue(y)) ]")
    end

    # Run tests
    test_cornerChecker()

    # Solve it
    solveProblem()

This code can also be found in ``/JuMP/examples/simplelazy2.jl``.