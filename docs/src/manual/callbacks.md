```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Solver-independent Callbacks](@id callbacks_manual)

Some mixed-integer (linear, conic, and nonlinear) programming solvers offer
the ability to modify the solve process. Examples include changing branching
decisions in branch-and-bound, adding custom cutting planes, providing custom
heuristics to find feasible solutions, or implementing on-demand separators to
add new constraints only when they are violated by the current solution (also
known as lazy constraints).

While historically this functionality has been limited to solver-specific
interfaces, JuMP provides solver-independent support for three types of
callbacks:

 1. lazy constraints
 2. user-cuts
 3. heuristic solutions

Even though JuMP provides a solver-independent way of writing these three
callbacks, you should not assume that you will see identical behavior when
running the same code on different solvers. For example, some solvers may ignore
user-cuts for various reasons, while other solvers may add every user-cut. Read
the underlying solver's callback documentation to understand details specific to
each solver.

## Available solvers

Solver-independent callback support is limited to a few solvers, including:

 * [CPLEX](https://github.com/jump-dev/CPLEX.jl)
 * [GLPK](https://github.com/jump-dev/GLPK.jl)
 * [Gurobi](https://github.com/jump-dev/Gurobi.jl)
 * [SCIP](https://github.com/scipopt/SCIP.jl) (SCIP does not support lazy
   constraints).
 * [Xpress](https://github.com/jump-dev/Xpress.jl)

Each solver listed above also provides a solver-_dependent_ callback to provide
access to the full range of solver-specific features. Consult the solver's
README for an example of how to use the solver-dependent callback. This will
require you to understand the C interface of the solver.

## Things you can and cannot do during solver-independent callbacks

There is a limited range of things you can do during a callback. Only use the
functions and macros explicitly stated in this page of the documentation, or in
the [Callbacks tutorial](@ref callbacks_tutorial).

Using any other part of the JuMP API (for example, adding a constraint with [`@constraint`](@ref)
or modifying a variable bound with [`set_lower_bound`](@ref)) is undefined
behavior, and your solver may throw an error, return an incorrect solution, or
result in a segfault that aborts Julia.

In each of the three solver-independent callbacks, there are two things you may
query:
 - [`callback_node_status`](@ref) returns an [`MOI.CallbackNodeStatusCode`](@ref)
   enum indicating if the current primal solution is integer feasible.
 - [`callback_value`](@ref) returns the current primal solution of a variable.

If you need to query any other information, use a solver-dependent callback
instead. Each solver supporting a solver-dependent callback has information on
how to use it in the README of their GitHub repository.

You can only set each callback once. Calling [`set_attribute`](@ref) twice will
over-write the earlier callback. In addition, if you use a solver-independent
callback, you cannot set a solver-dependent callback.

## Lazy constraints

Lazy constraints are useful when the full set of constraints is too large to
explicitly include in the initial formulation. When a MIP solver reaches a new
solution, for example with a heuristic or by solving a problem at a node in
the branch-and-bound tree, it will give the user the chance to provide
constraints that would make the current solution infeasible. For some more
information about lazy constraints, see this [blog post by Paul Rubin](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).

A lazy constraint callback can be set using the following syntax:

```jldoctest
julia> model = Model();

julia> @variable(model, x <= 10, Int)
x

julia> @objective(model, Max, x)
x

julia> function my_callback_function(cb_data)
           status = callback_node_status(cb_data, model)
           if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
               # `callback_value(cb_data, x)` is not integer (to some tolerance).
               # If, for example, your lazy constraint generator requires an
               # integer-feasible primal solution, you can add a `return` here.
               return
           elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
               # `callback_value(cb_data, x)` is integer (to some tolerance).
           else
               @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
               # `callback_value(cb_data, x)` might be fractional or integer.
           end
           x_val = callback_value(cb_data, x)
           if x_val > 2 + 1e-6
               con = @build_constraint(x <= 2)
               MOI.submit(model, MOI.LazyConstraint(cb_data), con)
           end
           return
       end
my_callback_function (generic function with 1 method)

julia> set_attribute(model, MOI.LazyConstraintCallback(), my_callback_function)
```

### Notes

 * The lazy constraint callback _may_ be called at fractional or integer nodes
   in the branch-and-bound tree. Use [`callback_node_status`](@ref) to
   distinguish.

 * There is no guarantee that the callback is called at _every_ primal solution.

 * The solver may visit a point that was cut off by a previous lazy constraint,
   for example, because the earlier lazy constraint was removed during presolve.
   If this happens, you must re-add the lazy constraint.

 * Only add a lazy constraint if the primal solution violates the constraint.
   Adding the lazy constraint irrespective of feasibility may result in the
   solver returning an incorrect solution, or lead to many constraints being
   added, slowing down the solution process. For example, instead of:
   ```julia
   model = Model()
   @variable(model, x <= 10, Int)
   @objective(model, Max, x)
   function bad_callback_function(cb_data)
       con = @build_constraint(x <= 2)
       MOI.submit(model, MOI.LazyConstraint(cb_data), con)
       return
   end
   ```
   do
   ```julia
   function good_callback_function(cb_data)
       if callback_value(cb_data, x) > 2
           con = @build_constraint(x <= 2)
           MOI.submit(model, MOI.LazyConstraint(cb_data), con)
       end
       return
   end
   set_attribute(model, MOI.LazyConstraintCallback(), good_callback_function)
   ```

## User cuts

User cuts, or simply cuts, provide a way for the user to tighten the LP
relaxation using problem-specific knowledge that the solver cannot or is
unable to infer from the model. Just like with lazy constraints, when a MIP
solver reaches a new node in the branch-and-bound tree, it will give the user
the chance to provide cuts to make the current relaxed (fractional) solution
infeasible in the hopes of obtaining an integer solution. For more details
about the difference between user cuts and lazy constraints see the
aforementioned [blog post](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).

A user-cut callback can be added using the following syntax:

```jldoctest
julia> model = Model();

julia> @variable(model, x <= 10.5, Int)
x

julia> @objective(model, Max, x)
x

julia> function my_callback_function(cb_data)
           x_val = callback_value(cb_data, x)
           con = @build_constraint(x <= floor(x_val))
           MOI.submit(model, MOI.UserCut(cb_data), con)
           return
       end
my_callback_function (generic function with 1 method)

julia> set_attribute(model, MOI.UserCutCallback(), my_callback_function)
```

### Notes

 * User cuts **must not** change the set of integer feasible solutions.
   Equivalently, user cuts can only remove fractional solutions. If you add a
   cut that removes an integer solution (even one that is not optimal), the
   solver may return an incorrect solution.

 * The user-cut callback _may_ be called at fractional nodes in the
   branch-and-bound tree. There is no guarantee that the callback is called at
   _every_ fractional primal solution.

## Heuristic solutions

Integer programming solvers frequently include heuristics that run at the
nodes of the branch-and-bound tree. They aim to find integer solutions quicker
than plain branch-and-bound would to tighten the bound, allowing us to fathom
nodes quicker and to tighten the integrality gap.

Some heuristics take integer solutions and explore their "local neighborhood"
(for example, flipping binary variables, fix some variables and solve a smaller
MILP) and others take fractional solutions and attempt to round them in an
intelligent way.

You may want to add a heuristic of your own if you have some special insight
into the problem structure that the solver is not aware of, for example, you can
consistently take fractional solutions and intelligently guess integer
solutions from them.

A heuristic solution callback can be added using the following syntax:

```jldoctest
julia> model = Model();

julia> @variable(model, x <= 10.5, Int);

julia> @objective(model, Max, x);

julia> function my_callback_function(cb_data)
           status = MOI.submit(
               model,
               MOI.HeuristicSolution(cb_data),
               [x],
               # Heuristic solution by rounding down to nearest integer
               Float64[floor(Int, callback_value(cb_data, x))],
           )
           println(
               "I submitted a heuristic solution, and the status was: ",
               status,
           )
           return
       end
my_callback_function (generic function with 1 method)

julia> set_attribute(model, MOI.HeuristicCallback(), my_callback_function)
```

[`MOI.submit`](@ref) returns a [`MOI.HeuristicSolutionStatus`](@ref) enum that
indicates whether the solver accepted the solution. The possible return codes
are:

 - [`MOI.HEURISTIC_SOLUTION_ACCEPTED`](@ref)
 - [`MOI.HEURISTIC_SOLUTION_REJECTED`](@ref)
 - [`MOI.HEURISTIC_SOLUTION_UNKNOWN`](@ref)

### Notes

 * Some solvers may accept partial solutions. Others require a feasible integer
   solution for every variable. If in doubt, provide a complete solution.

 * The heuristic solution callback _may_ be called at fractional nodes in the
   branch-and-bound tree. There is no guarantee that the callback is called at
   _every_ fractional primal solution.
