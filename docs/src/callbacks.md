```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Callbacks

Many mixed-integer (linear, conic, and nonlinear) programming solvers offer
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

## Available solvers

Callback support is limited to a few solvers. This includes
[CPLEX](https://github.com/JuliaOpt/CPLEX.jl),
[GLPK](https://github.com/JuliaOpt/GLPK.jl), and
[Gurobi](https://github.com/JuliaOpt/Gurobi.jl).

!!! warning
    While JuMP provides a solver-independent way of accessing callbacks, you
    should not assume that you will see identical behavior when running the same
    code on different solvers. For example, some solvers may ignore user-cuts
    for various reasons, while other solvers may add every user-cut. Read the
    underlying solver's callback documentation to understand details specific to
    each solver.

## Things you can and cannot do during callbacks

There is a very limited range of things you can do during a callback. Only use 
the functions and macros explicitly stated in this page of the documentation, or
in the [Callbacks example](/examples/callbacks).

Using any other part of the JuMP API (e.g., adding a constraint with [`@constraint`](@ref)
or modifying a variable bound with [`set_lower_bound`](@ref)) is undefined
behavior, and your solver may throw an error, return an incorrect solution, or 
result in a segfault that aborts Julia.

In each of the three solver-independent callbacks, the only thing you may query is 
the primal value of the variables using [`callback_value`](@ref).

If you need to query any other information, use a solver-dependent callback 
instead. Each solver supporting a solver-dependent callback has information on 
how to use it in the README of their Github repository.

If you want to modify the problem in a callback, you _must_ use a lazy 
constraint.

## Lazy constraints

Lazy constraints are useful when the full set of constraints is too large to
explicitly include in the initial formulation. When a MIP solver reaches a new
solution, for example with a heuristic or by solving a problem at a node in
the branch-and-bound tree, it will give the user the chance to provide
constraint(s) that would make the current solution infeasible. For some more
information about lazy constraints, see this [blog post by Paul Rubin](https://orinanobworld.blogspot.com/2012/08/user-cuts-versus-lazy-constraints.html).

A lazy constraint callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    x_val = callback_value(cb_data, x)
    if x_val > 2 + 1e-6
        con = @build_constraint(x <= 2)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
end
MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function)
```

!!! info
    The lazy constraint callback _may_ be called at fractional or integer
    nodes in the branch-and-bound tree. There is no guarantee that the
    callback is called at _every_ primal solution.
    
!!! warn
    Only add a lazy constraint if your primal solution violates the constraint.
    Adding the lazy constraint irrespective of feasibility may result in the 
    solver returning an incorrect solution, or lead to a large number of 
    constraints being added, slowing down the solution process.
    ```julia
    model = Model(GLPK.Optimizer)
    @variable(model, x <= 10, Int)
    @objective(model, Max, x)
    function bad_callback_function(cb_data)
        # Don't do this!
        con = @build_constraint(x <= 2)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
    function good_callback_function(cb_data)
        if callback_value(x) > 2
            con = @build_constraint(x <= 2)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(model, MOI.LazyConstraintCallback(), good_callback_function)
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

A user-cut callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10.5, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    x_val = callback_value(cb_data, x)
    con = @build_constraint(x <= floor(x_val))
    MOI.submit(model, MOI.UserCut(cb_data), con)
end
MOI.set(model, MOI.UserCutCallback(), my_callback_function)
```

!!! warning
    Your user cuts should not change the set of integer feasible solutions.
    Equivalently, your cuts can only remove fractional solutions. If you add a
    cut that removes an integer solution, the solver may return an incorrect
    solution.

!!! info
    The user-cut callback _may_ be called at fractional nodes in the
    branch-and-bound tree. There is no guarantee that the callback is called
    at _every_ fractional primal solution.

## Heuristic solutions

Integer programming solvers frequently include heuristics that run at the
nodes of the branch-and-bound tree. They aim to find integer solutions quicker
than plain branch-and-bound would to tighten the bound, allowing us to fathom
nodes quicker and to tighten the integrality gap.

Some heuristics take integer solutions and explore their "local neighborhood"
(e.g., flipping binary variables, fix some variables and solve a smaller MILP)
and others take fractional solutions and attempt to round them in an
intelligent way.

You may want to add a heuristic of your own if you have some special insight
into the problem structure that the solver is not aware of, e.g. you can
consistently take fractional solutions and intelligently guess integer
solutions from them.

A heuristic solution callback can be set using the following syntax:

```julia
model = Model(GLPK.Optimizer)
@variable(model, x <= 10.5, Int)
@objective(model, Max, x)
function my_callback_function(cb_data)
    x_val = callback_value(cb_data, x)
    status = MOI.submit(
        model, MOI.HeuristicSolution(cb_data), [x], [floor(Int, x_val)]
    )
    println("I submitted a heuristic solution, and the status was: ", status)
end
MOI.set(model, MOI.HeuristicCallback(), my_callback_function)
```

The third argument to `submit` should be a vector of JuMP variables, and the
fourth argument should be a vector of values corresponding to each variable.

`MOI.submit` returns an enum that depends on whether the solver accepted the
solution. The possible return codes are:

 - `MOI.HEURISTIC_SOLUTION_ACCEPTED`
 - `MOI.HEURISTIC_SOLUTION_REJECTED`
 - `MOI.HEURISTIC_SOLUTION_UNKNOWN`

!!! warning
    Some solvers may accept partial solutions. Others require a feasible integer
    solution for every variable. If in doubt, provide a complete solution.

!!! info
    The heuristic solution callback _may_ be called at fractional nodes in the
    branch-and-bound tree. There is no guarantee that the callback is called
    at _every_ fractional primal solution.

## Reference

```@docs
callback_value
```
