.. _callbacks:

.. include:: warn.rst

----------------
Solver Callbacks
----------------

Many mixed-integer programming solvers offer the ability to modify the solve process.
Examples include changing branching decisions in branch-and-bound, adding custom cutting planes, providing custom heuristics to find feasible solutions, or implementing on-demand separators to add new constraints only when they are violated by the current solution (also known as lazy constraints).

While historically this functionality has been limited to solver-specific interfaces,
JuMP provides *solver-independent* support for a number of commonly used solver callbacks. Currently, we support lazy constraints, user-provided cuts, and user-provided
heuristics for the Gurobi, CPLEX, and GLPK solvers. We do not yet support any
other class of callbacks, but they may be accessible by using the solver's
low-level interface.

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
single argument, e.g. ``function myLazyConGenerator(cb)``, where cb is a reference
to the callback management code inside JuMP. Next you will do whatever
analysis of the solution you need to inside your function to generate the new
constraint before adding it to the model with
``@lazyconstraint(cb, myconstraint)``.
Finally we notify JuMP that this function should be used for lazy constraint
generation using the ``addlazycallback(m, myLazyConGenerator)`` function
before we call ``solve(m)``.

The following is a simple example to make this more clear. In this two-dimensional
problem we have a set of box constraints explicitly provided and a set of two
lazy constraints we can add on the fly. The solution without the lazy constraints
will be either (0,2) or (2,2), and the final solution will be (1,2)::

    using JuMP
    using Gurobi

    # We will use Gurobi
    m = Model(solver=GurobiSolver())

    # Define our variables to be inside a box, and integer
    @variable(m, 0 <= x <= 2, Int)
    @variable(m, 0 <= y <= 2, Int)

    @objective(m, Max, y)

    # We now define our callback function that takes one argument,
    # the callback handle. Note that we can access m, x, and y because
    # this function is defined inside the same scope
    function corners(cb)
        x_val = getvalue(x)
        y_val = getvalue(y)
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
            @lazyconstraint(cb, y - x <= 1)
        # Check top right
        elseif y_val + x_val > 3 + TOL
            # Cut off this solution
            println("Solution was in top right, cut it off")
            # Use the original variables, but not m - cb instead
            @lazyconstraint(cb, y + x <= 3)
        end
    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    addlazycallback(m, corners)

    # Solve the problem
    solve(m)

    # Print our final solution
    println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

The code should print something like (amongst the output from Gurobi)::

    In callback function, x=2.0, y=2.0
    Solution was in top right, cut it off
    In callback function, x=0.0, y=2.0
    Solution was in top left, cut it off
    In callback function, x=1.0, y=2.0
    Final solution: [ 1.0, 2.0 ]

This code can also be found in ``/JuMP/examples/simplelazy.jl``.

There is an optional ``fractional`` keyword option to ``addlazycallback`` which
indicates that the callback may be called at solutions that do not satisfy
integrality constraints. For example, ``addlazycallback(m, myLazyConGenerator,
fractional=true)``. Depending on the solver, this may invoke the callback
after solving each LP relaxation in the Branch and Bound tree. By default, ``fractional`` is set to ``false``.


User Cuts
^^^^^^^^^

User cuts, or simply cuts, provide a way for the user to tighten the LP relaxation using problem-specific knowledge that the solver cannot or is unable to infer from the model. Just like with lazy constraints, when a MIP solver reaches a new node in the branch-and-bound tree, it will give the user the chance to provide cuts to make the current relaxed (fractional) solution infeasible in the hopes of obtaining an integer solution. For more details about the difference between user cuts and lazy constraints see the aforementioned `blog post <http://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html>`_.

Your user cuts should not change the set of integer feasible solutions. Equivalently, your cuts can only remove fractional solutions - that is, "tighten" the LP relaxation of the MILP. If you add a cut that removes an integer solution, the solver may return an incorrect solution.

Adding a user cut callback is similar to adding a lazy constraint callback. First we
must write a function that will analyze the current solution that takes a
single argument, e.g. ``function myUserCutGenerator(cb)``, where cb is a reference
to the callback management code inside JuMP. Next you will do whatever
analysis of the solution you need to inside your function to generate the new
constraint before adding it to the model with the JuMP macro
``@addusercut(cb, myconstraint)`` (same limitations as addConstraint).
Finally we notify JuMP that this function should be used for lazy constraint
generation using the ``addcutcallback(m, myUserCutGenerator)`` function
before we call ``solve(m)``.

Consider the following example which is related to the lazy constraint example. The problem is two-dimensional, and the objective sense prefers solution in the top-right of a 2-by-2 square. There is a single constraint that cuts off the top-right corner to make the LP relaxation solution fractional. We will exploit our knowledge of the problem structure to add a user cut that will make the LP relaxation integer, and thus solve the problem at the root node::

    using JuMP
    using Gurobi

    # We will use Gurobi, which requires that we manually set the attribute
    # PreCrush to 1 if we have user cuts. We will also disable PreSolve, Cuts,
    # and Heuristics so only our cut will be used
    m = Model(solver=GurobiSolver(PreCrush=1, Cuts=0, Presolve=0, Heuristics=0.0))

    # Define our variables to be inside a box, and integer
    @variable(m, 0 <= x <= 2, Int)
    @variable(m, 0 <= y <= 2, Int)

    # Optimal solution is trying to go towards top-right corner (2.0, 2.0)
    @objective(m, Max, x + 2y)

    # We have one constraint that cuts off the top right corner
    @constraint(m, y + x <= 3.5)

    # Optimal solution of relaxed problem will be (1.5, 2.0)
    # We can add a user cut that will cut of this fractional solution.

    # We now define our callback function that takes one argument,
    # the callback handle. Note that we can access m, x, and y because
    # this function is defined inside the same scope
    function mycutgenerator(cb)
        x_val = getvalue(x)
        y_val = getvalue(y)
        println("In callback function, x=$x_val, y=$y_val")

        # Allow for some impreciseness in the solution
        TOL = 1e-6

        # Check top right
        if y_val + x_val > 3 + TOL
            # Cut off this solution
            println("Fractional solution was in top right, cut it off")
            # Use the original variables
            @addusercut(cb, y + x <= 3)
        end
    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    addcutcallback(m, mycutgenerator)

    # Solve the problem
    solve(m)

    # Print our final solution
    println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

The code should print something like (amongst the output from Gurobi)::

    In callback function, x=1.5, y=2.0
    Fractional solution was in top right, cut it off
    In callback function, x=1.0, y=2.0
    Final solution: [ 1.0, 2.0 ]

This code can also be found in ``/JuMP/examples/simpleusercut.jl``.


User Heuristics
^^^^^^^^^^^^^^^

Integer programming solvers frequently include heuristics that run at the nodes of the branch-and-bound tree. They aim to find integer solutions quicker than plain branch-and-bound would to tighten the bound, allowing us to fathom nodes quicker and to tighten the integrality gap. Some heuristics take integer solutions and explore their "local neighborhood" (e.g. flipping binary variables, fix some variables and solve a smaller MILP, ...) and others take fractional solutions and attempt to round them in an intelligent way. You may want to add a heuristic of your own if you have some special insight into the problem structure that the solver is not aware of, e.g. you can consistently take fractional solutions and intelligently guess integer solutions from them.

The user heuristic callback is somewhat different from the previous two heuristics. The general concept is that we can create multiple partial solutions and submit them back to the solver - each solution must be submitted before a new solution is constructed. As before we provide a function that analyzes the current solution and takes a single argument, e.g. ``function myHeuristic(cb)``, where cb is a reference to the callback management code inside JuMP. You can build your solutions using ``setsolutionvalue(cb, x, value)`` and submit them with ``addsolution(cb)``. Note that ``addsolution`` will "wipe" the previous (partial) solution. Notify JuMP that this function should be used as a heuristic using the ``addheuristiccallback(m, myHeuristic)`` function before calling ``solve(m)``.

There is some unavoidable (for performance reasons) solver-dependent behavior - you should check your solver documentation for details. For example: GLPK will not check the feasibility of your heuristic solution. If you need to submit many heuristic solutions in one callback, there may be performance impacts from the "wiping" behavior of ``addsolution`` - please file an issue and we can address this issue.

Consider the following example, which is the same problem as seen in the user cuts section. The heuristic simply rounds the fractional variable to generate integer solutions.::

    using JuMP
    using Gurobi

    # We will use Gurobi and disable PreSolve, Cuts, and (in-built) Heuristics so
    # only our heuristic will be used
    m = Model(solver=GurobiSolver(Cuts=0, Presolve=0, Heuristics=0.0))

    # Define our variables to be inside a box, and integer
    @variable(m, 0 <= x <= 2, Int)
    @variable(m, 0 <= y <= 2, Int)

    # Optimal solution is trying to go towards top-right corner (2.0, 2.0)
    @objective(m, Max, x + 2y)

    # We have one constraint that cuts off the top right corner
    @constraint(m, y + x <= 3.5)

    # Optimal solution of relaxed problem will be (1.5, 2.0)

    # We now define our callback function that takes one argument,
    # the callback handle. Note that we can access m, x, and y because
    # this function is defined inside the same scope
    function myheuristic(cb)
        x_val = getvalue(x)
        y_val = getvalue(y)
        println("In callback function, x=$x_val, y=$y_val")

        setsolutionvalue(cb, x, floor(x_val))
        # Leave y undefined - solver should handle as it sees fit. In the case
        # of Gurobi it will try to figure out what it should be.
        addsolution(cb)

        # Submit a second solution
        setsolutionvalue(cb, x, ceil(x_val))
        addsolution(cb)
    end  # End of callback function

    # Tell JuMP/Gurobi to use our callback function
    addheuristiccallback(m, myheuristic)

    # Solve the problem
    solve(m)

    # Print our final solution
    println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")

The code should print something like::

    In callback function, x=1.5, y=2.0
         0     0    5.50000    0    1          -    5.50000     -      -    0s
    H    1     0                       5.0000000    5.50000  10.0%   0.0    0s

where the ``H`` denotes a solution found with a heuristic - our heuristic in this case. This code can also be found in ``/JuMP/examples/simpleheur.jl``.



Querying Solver Progress
^^^^^^^^^^^^^^^^^^^^^^^^

All JuMP callback methods must take a single argument, called ``cb`` by convention.
``cb`` is a handle to the internal callback system used by the underlying solver, and
allows the user to query solver state. There are a variety of methods available which
are listed in the `MathProgBase documentation <http://mathprogbasejl.readthedocs.org/en/latest/lpqcqp.html#mip-callbacks>`_
including::

    cbgetobj(cb)
    cbgetbestbound(cb)
    cbgetexplorednodes(cb)
    cbgetstate(cb)


Informational Callbacks
^^^^^^^^^^^^^^^^^^^^^^^

Sometimes it can be useful to track solver progress without actually changing the algorithm by adding cuts or heuristic solutions. In these cases, informational callbacks can be added, wherein statistics can be tracked via the ``cbget`` functions discussed in the previous section. Informational callbacks are added to a JuMP model with the ``addinfocallback(m::Model, f::Function)`` function.

For a simple example, we can add a function that tracks the best bound and incumbent objective value as the solver progresses through the branch-and-bound tree::

    type NodeData
        time::Float64  # in seconds since the epoch
        node::Int
        obj::Float64
        bestbound::Float64
    end

    # build model ``m`` up here

    bbdata = NodeData[]

    function infocallback(cb)
        node      = MathProgBase.cbgetexplorednodes(cb)
        obj       = MathProgBase.cbgetobj(cb)
        bestbound = MathProgBase.cbgetbestbound(cb)
        push!(bbdata, NodeData(time(),node,obj,bestbound))
    end
    addinfocallback(m, infocallback)

    solve(m)

    # Save results to file for analysis later
    open("bbtrack.csv","w") do fp
        println(fp, "time,node,obj,bestbound")
        for bb in bbdata
            println(fp, bb.time, ",", bb.node, ",",
                        bb.obj, ",", bb.bestbound)
        end
    end


Code Design Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above examples the callback function is defined in the same scope as the model and variable definitions, allowing us to access them. If we defined the function in some other scope, or even file, we would not be able to access them directly. The proposed solution to this design problem is to separate the logic of analyzing the current solution values from the callback itself. This has many benefits, including writing unit tests for the callback function to check its correctness. The callback function passed to JuMP is then simply a stub that extracts the current solution and any other relevant information and passes that to the constraint generation logic. To apply this to our previous lazy constraint example, consider the following code::

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
        m = Model(solver=GurobiSolver())

        @variable(m, 0 <= x <= 2, Int)
        @variable(m, 0 <= y <= 2, Int)
        @objective(m, Max, y)

        # Note that the callback is now a stub that passes off
        # the work to the "algorithm"
        function corners(cb)
            x_val = getvalue(x)
            y_val = getvalue(y)
            println("In callback function, x=$x_val, y=$y_val")

            newcut, x_coeff, y_coeff, rhs = cornerChecker(x_val, y_val)

            if newcut
                @lazyconstraint(cb, x_coeff*x + y_coeff*y <= rhs)
            end
        end  # End of callback function

        addlazycallback(m, corners)
        solve(m)
        println("Final solution: [ $(getvalue(x)), $(getvalue(y)) ]")
    end

    # Run tests
    test_cornerChecker()

    # Solve it
    solveProblem()

This code can also be found in ``/JuMP/examples/simplelazy2.jl``.

Exiting a callback early
^^^^^^^^^^^^^^^^^^^^^^^^

If you need to exit the optimization process earlier than a solver otherwise would, it is possible to throw a ``CallbackAbort`` exception in callback code::

    throw(CallbackAbort())

This will trigger the solver to exit immediately and return a ``:UserLimit`` status.
