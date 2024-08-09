```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, HiGHS
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Solutions](@id jump_solutions)

This section of the manual describes how to access a solved solution to a
problem. It uses the following model as an example:
```jldoctest solutions
julia> begin
           model = Model(HiGHS.Optimizer)
           set_silent(model)
           @variable(model, x >= 0)
           @variable(model, y[[:a, :b]] <= 1)
           @objective(model, Max, -12x - 20y[:a])
           @expression(model, my_expr, 6x + 8y[:a])
           @constraint(model, my_expr >= 100)
           @constraint(model, c1, 7x + 12y[:a] >= 120)
           optimize!(model)
           print(model)
       end
Max -12 x - 20 y[a]
Subject to
 6 x + 8 y[a] ≥ 100
 c1 : 7 x + 12 y[a] ≥ 120
 x ≥ 0
 y[a] ≤ 1
 y[b] ≤ 1
```

## Check if an optimal solution exists

Use [`is_solved_and_feasible`](@ref) to check if the solver found an optimal
solution:
```jldoctest solutions
julia> is_solved_and_feasible(model)
true
```

By default, [`is_solved_and_feasible`](@ref) returns `true` for both global and
local optima. Pass `allow_local = false` to check if the solver found a globally
optimal solution:
```jldoctest solutions
julia> is_solved_and_feasible(model; allow_local = false)
true
```

Pass `dual = true` to check if the solver found an optimal dual solution in
addition to an optimal primal solution:
```jldoctest solutions
julia> is_solved_and_feasible(model; dual = true)
true
```

If this function returns `false`, use the functions mentioned below like
[`solution_summary`](@ref), [`termination_status`](@ref), [`primal_status`](@ref),
and [`dual_status`](@ref) to understand what solution (if any) the solver found.

## Solutions summary

[`solution_summary`](@ref) can be used for checking the summary of the
optimization solutions.

```jldoctest solutions; filter=r"[0-9]+\.[0-9]+e[\+\-][0-9]+"
julia> solution_summary(model)
* Solver : HiGHS

* Status
  Result count       : 1
  Termination status : OPTIMAL
  Message from the solver:
  "kHighsModelStatusOptimal"

* Candidate solution (result #1)
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Objective value    : -2.05143e+02
  Objective bound    : -0.00000e+00
  Relative gap       : Inf
  Dual objective value : -2.05143e+02

* Work counters
  Solve time (sec)   : 6.70987e-04
  Simplex iterations : 2
  Barrier iterations : 0
  Node count         : -1

julia> solution_summary(model; verbose = true)
* Solver : HiGHS

* Status
  Result count       : 1
  Termination status : OPTIMAL
  Message from the solver:
  "kHighsModelStatusOptimal"

* Candidate solution (result #1)
  Primal status      : FEASIBLE_POINT
  Dual status        : FEASIBLE_POINT
  Objective value    : -2.05143e+02
  Objective bound    : -0.00000e+00
  Relative gap       : Inf
  Dual objective value : -2.05143e+02
  Primal solution :
    x : 1.54286e+01
    y[a] : 1.00000e+00
    y[b] : 1.00000e+00
  Dual solution :
    c1 : 1.71429e+00

* Work counters
  Solve time (sec)   : 6.70987e-04
  Simplex iterations : 2
  Barrier iterations : 0
  Node count         : -1
```

## Why did the solver stop?

Use[`termination_status`](@ref) to understand why the solver stopped.

```jldoctest solutions
julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1
```

The [`MOI.TerminationStatusCode`](@ref) enum describes the full list of statuses
that could be returned.

Common return values include [`OPTIMAL`](@ref), [`LOCALLY_SOLVED`](@ref),
[`INFEASIBLE`](@ref), [`DUAL_INFEASIBLE`](@ref), and [`TIME_LIMIT`](@ref).

!!! info
    A return status of [`OPTIMAL`](@ref) means the solver found (and proved) a
    globally optimal solution. A return status of [`LOCALLY_SOLVED`](@ref) means
    the solver found a locally optimal solution (which may also be globally
    optimal, but it could not prove so).

!!! warning
    A return status of [`DUAL_INFEASIBLE`](@ref) does not guarantee that the
    primal is unbounded. When the dual is infeasible, the primal is unbounded if
    there exists a feasible primal solution.

Use [`raw_status`](@ref) to get a solver-specific string explaining why the
optimization stopped:
```jldoctest solutions
julia> raw_status(model)
"kHighsModelStatusOptimal"
```

## Primal solutions

### Primal solution status

Use [`primal_status`](@ref) to return an [`MOI.ResultStatusCode`](@ref) enum
describing the status of the primal solution.
```jldoctest solutions
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Other common returns are [`NO_SOLUTION`](@ref), and [`INFEASIBILITY_CERTIFICATE`](@ref).
The first means that the solver doesn't have a solution to return, and the
second means that the primal solution is a certificate of dual infeasibility (a
primal unbounded ray).

You can also use [`has_values`](@ref), which returns `true` if there is a
solution that can be queried, and `false` otherwise.
```jldoctest solutions
julia> has_values(model)
true
```

### Objective values

The objective value of a solved problem can be obtained via
[`objective_value`](@ref). The best known bound on the optimal objective
value can be obtained via [`objective_bound`](@ref). If the solver supports it,
the value of the dual objective can be obtained via
[`dual_objective_value`](@ref).

```jldoctest solutions
julia> objective_value(model)
-205.14285714285714

julia> objective_bound(model)  # HiGHS only implements objective bound for MIPs
-0.0

julia> dual_objective_value(model)
-205.1428571428571
```

### Primal solution values

If the solver has a primal solution to return, use [`value`](@ref) to access it:
```jldoctest solutions
julia> value(x)
15.428571428571429
```

Broadcast [`value`](@ref) over containers:
```jldoctest solutions
julia> value.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, [:a, :b]
And data, a 2-element Vector{Float64}:
 1.0
 1.0
```

[`value`](@ref) also works on expressions:
```jldoctest solutions
julia> value(my_expr)
100.57142857142857
```
and constraints:
```jldoctest solutions
julia> value(c1)
120.0
```
!!! info
    Calling [`value`](@ref) on a constraint returns the constraint function
    evaluated at the solution.

## Dual solutions

### Dual solution status

Use [`dual_status`](@ref) to return an [`MOI.ResultStatusCode`](@ref) enum
describing the status of the dual solution.
```jldoctest solutions
julia> dual_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Other common returns are [`NO_SOLUTION`](@ref), and [`INFEASIBILITY_CERTIFICATE`](@ref).
The first means that the solver doesn't have a solution to return, and the
second means that the dual solution is a certificate of primal infeasibility (a
dual unbounded ray).

You can also use [`has_duals`](@ref), which returns `true` if there is a
solution that can be queried, and `false` otherwise.
```jldoctest solutions
julia> has_duals(model)
true
```

### Dual solution values

If the solver has a dual solution to return, use [`dual`](@ref) to access it:
```jldoctest solutions
julia> dual(c1)
1.7142857142857142
```

Query the duals of variable bounds using [`LowerBoundRef`](@ref),
[`UpperBoundRef`](@ref), and [`FixRef`](@ref):
```jldoctest solutions
julia> dual(LowerBoundRef(x))
0.0

julia> dual.(UpperBoundRef.(y))
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, [:a, :b]
And data, a 2-element Vector{Float64}:
 -0.5714285714285694
  0.0
```

!!! warning
    JuMP's definition of duality is independent of the objective sense. That is,
    the sign of feasible duals associated with a constraint depends on the
    direction of the constraint and not whether the problem is maximization or
    minimization. **This is a different convention from linear programming
    duality in some common textbooks.** If you have a linear program, and you
    want the textbook definition, you probably want to use [`shadow_price`](@ref)
    and [`reduced_cost`](@ref) instead.

```jldoctest solutions
julia> shadow_price(c1)
1.7142857142857142

julia> reduced_cost(x)
-0.0

julia> reduced_cost.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, [:a, :b]
And data, a 2-element Vector{Float64}:
  0.5714285714285694
 -0.0
```

## Recommended workflow

You should always check whether the solver found a solution before calling
solution functions like [`value`](@ref) or [`objective_value`](@ref).

A simple approach for small scripts and notebooks is to use
[`is_solved_and_feasible`](@ref):

```jldoctest solutions
julia> function solve_and_print_solution(model)
           optimize!(model)
           if !is_solved_and_feasible(model; dual = true)
               error(
                   """
                   The model was not solved correctly:
                   termination_status : $(termination_status(model))
                   primal_status      : $(primal_status(model))
                   dual_status        : $(dual_status(model))
                   raw_status         : $(raw_status(model))
                   """,
               )
           end
           println("Solution is optimal")
           println("  objective value = ", objective_value(model))
           println("  primal solution: x = ", value(x))
           println("  dual solution: c1 = ", dual(c1))
           return
       end
solve_and_print_solution (generic function with 1 method)

julia> solve_and_print_solution(model)
Solution is optimal
  objective value = -205.14285714285714
  primal solution: x = 15.428571428571429
  dual solution: c1 = 1.7142857142857142
```

For code like libraries that should be more robust to the range of possible
termination and result statuses, do some variation of the following:
```jldoctest solutions
julia> function solve_and_print_solution(model)
           status = termination_status(model)
           if status in (OPTIMAL, LOCALLY_SOLVED)
               println("Solution is optimal")
           elseif status in (ALMOST_OPTIMAL, ALMOST_LOCALLY_SOLVED)
               println("Solution is optimal to a relaxed tolerance")
           elseif status == TIME_LIMIT
               println(
                   "Solver stopped due to a time limit. If a solution is available, " *
                   "it may be suboptimal."
               )
           elseif status in (
               ITERATION_LIMIT, NODE_LIMIT, SOLUTION_LIMIT, MEMORY_LIMIT,
               OBJECTIVE_LIMIT, NORM_LIMIT, OTHER_LIMIT,
           )
               println(
                   "Solver stopped due to a limit. If a solution is available, it " *
                   "may be suboptimal."
               )
           elseif status in (INFEASIBLE, LOCALLY_INFEASIBLE)
               println("The problem is primal infeasible")
           elseif status == DUAL_INFEASIBLE
               println(
                   "The problem is dual infeasible. If a primal feasible solution " *
                   "exists, the problem is unbounded. To check, set the objective " *
                   "to `@objective(model, Min, 0)` and re-solve. If the problem is " *
                   "feasible, the primal is unbounded. If the problem is " *
                   "infeasible, both the primal and dual are infeasible.",
               )
           elseif status == INFEASIBLE_OR_UNBOUNDED
               println(
                   "The model is either infeasible or unbounded. Set the objective " *
                   "to `@objective(model, Min, 0)` and re-solve to disambiguate. If " *
                   "the problem was infeasible, it will still be infeasible. If the " *
                   "problem was unbounded, it will now have a finite optimal solution.",
               )
           else
               println(
                   "The model was not solved correctly. The termination status is $status",
               )
           end
           if primal_status(model) in (FEASIBLE_POINT, NEARLY_FEASIBLE_POINT)
               println("  objective value = ", objective_value(model))
               println("  primal solution: x = ", value(x))
           elseif primal_status(model) == INFEASIBILITY_CERTIFICATE
               println("  primal certificate: x = ", value(x))
           end
           if dual_status(model) in (FEASIBLE_POINT, NEARLY_FEASIBLE_POINT)
               println("  dual solution: c1 = ", dual(c1))
           elseif dual_status(model) == INFEASIBILITY_CERTIFICATE
               println("  dual certificate: c1 = ", dual(c1))
           end
           return
       end
solve_and_print_solution (generic function with 1 method)

julia> solve_and_print_solution(model)
Solution is optimal
  objective value = -205.14285714285714
  primal solution: x = 15.428571428571429
  dual solution: c1 = 1.7142857142857142
```

## OptimizeNotCalled errors

Due to differences in how solvers cache solutions internally, modifying a model
after calling [`optimize!`](@ref) will reset the model into the
[`OPTIMIZE_NOT_CALLED`](@ref) state. If you then attempt to query solution
information, an [`OptimizeNotCalled`](@ref) error will be thrown.

If you are iteratively querying solution information and modifying a model,
query all the results first, then modify the problem.

For example, instead of:
```jldoctest; filter = r"\@ JuMP.+/src/optimizer_interface.jl:[0-9]+"
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0);

julia> optimize!(model)

julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1

julia> set_upper_bound(x, 1)

julia> x_val = value(x)
┌ Warning: The model has been modified since the last call to `optimize!` (or `optimize!` has not been called yet). If you are iteratively querying solution information and modifying a model, query all the results first, then modify the model.
└ @ JuMP ~/work/JuMP.jl/JuMP.jl/src/optimizer_interface.jl:712
ERROR: OptimizeNotCalled()
Stacktrace:
[...]

julia> termination_status(model)
OPTIMIZE_NOT_CALLED::TerminationStatusCode = 0
```
do
```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0);

julia> optimize!(model);

julia> x_val = value(x)
0.0

julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1

julia> set_upper_bound(x, 1)

julia> set_lower_bound(x, x_val)

julia> termination_status(model)
OPTIMIZE_NOT_CALLED::TerminationStatusCode = 0
```

If you know that your particular solver supports querying solution information
after modifications, you can use [`direct_model`](@ref) to bypass the
[`OPTIMIZE_NOT_CALLED`](@ref) state:
```jldoctest
julia> model = direct_model(HiGHS.Optimizer());

julia> set_silent(model)

julia> @variable(model, x >= 0);

julia> optimize!(model)

julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1

julia> set_upper_bound(x, 1)

julia> x_val = value(x)
0.0

julia> set_lower_bound(x, x_val)

julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1
```

!!! warning
    Be careful doing this. If your particular solver does not support
    querying solution information after modification, it may silently return
    incorrect solutions or throw an error.

```@meta
# TODO: How to accurately measure the solve time.
```

## Accessing attributes

[MathOptInterface](@ref moi_documentation) defines many model attributes that
can be queried. Some attributes can be directly accessed by getter functions.
These include:
- [`solve_time`](@ref)
- [`relative_gap`](@ref)
- [`simplex_iterations`](@ref)
- [`barrier_iterations`](@ref)
- [`node_count`](@ref)

## Sensitivity analysis for LP

Given an LP problem and an optimal solution corresponding to a basis, we can
question how much an objective coefficient or standard form right-hand side
coefficient (c.f., [`normalized_rhs`](@ref)) can change without violating primal
or dual feasibility of the basic solution.

Note that not all solvers compute the basis, and for sensitivity analysis, the
solver interface must implement [`MOI.ConstraintBasisStatus`](@ref).

!!! tip
    Read the [Sensitivity analysis of a linear program](@ref) for more
    information on sensitivity analysis.

To give a simple example, we could analyze the sensitivity of the optimal
solution to the following (non-degenerate) LP problem:

```jldoctest solutions_sensitivity
julia> begin
           model = Model(HiGHS.Optimizer)
           set_silent(model)
           @variable(model, x[1:2])
           set_lower_bound(x[2], -0.5)
           set_upper_bound(x[2], 0.5)
           @constraint(model, c1, x[1] + x[2] <= 1)
           @constraint(model, c2, x[1] - x[2] <= 1)
           @objective(model, Max, x[1])
           print(model)
       end
Max x[1]
Subject to
 c1 : x[1] + x[2] ≤ 1
 c2 : x[1] - x[2] ≤ 1
 x[2] ≥ -0.5
 x[2] ≤ 0.5
```

To analyze the sensitivity of the problem we could check the allowed
perturbation ranges of, for example, the cost coefficients and the right-hand side
coefficient of the constraint `c1` as follows:

```jldoctest solutions_sensitivity
julia> optimize!(model)

julia> value.(x)
2-element Vector{Float64}:
  1.0
 -0.0

julia> report = lp_sensitivity_report(model);

julia> x1_lo, x1_hi = report[x[1]]
(-1.0, Inf)

julia> println("The objective coefficient of x[1] could decrease by $(x1_lo) or increase by $(x1_hi).")
The objective coefficient of x[1] could decrease by -1.0 or increase by Inf.

julia> x2_lo, x2_hi = report[x[2]]
(-1.0, 1.0)

julia> println("The objective coefficient of x[2] could decrease by $(x2_lo) or increase by $(x2_hi).")
The objective coefficient of x[2] could decrease by -1.0 or increase by 1.0.

julia> c_lo, c_hi = report[c1]
(-1.0, 1.0)

julia> println("The RHS of c1 could decrease by $(c_lo) or increase by $(c_hi).")
The RHS of c1 could decrease by -1.0 or increase by 1.0.
```

The range associated with a variable is the range of the allowed perturbation of
the corresponding objective coefficient. Note that the current primal solution
remains optimal within this range; however the corresponding dual solution might
change since a cost coefficient is perturbed. Similarly, the range associated
with a constraint is the range of the allowed perturbation of the corresponding
right-hand side coefficient. In this range the current dual solution remains
optimal, but the optimal primal solution might change.

If the problem is degenerate, there are multiple optimal bases and hence these
ranges might not be as intuitive and seem too narrow, for example, a larger cost
coefficient perturbation might not invalidate the optimality of the current
primal solution. Moreover, if a problem is degenerate, due to finite precision,
it can happen that, for example, a perturbation seems to invalidate a basis even though
it doesn't (again providing too narrow ranges). To prevent this, increase the
`atol` keyword argument to [`lp_sensitivity_report`](@ref). Note that this might
make the ranges too wide for numerically challenging instances. Thus, do not
blindly trust these ranges, especially not for highly degenerate or numerically
unstable instances.

## Conflicts

When the model you input is infeasible, some solvers can help you find the
cause of this infeasibility by offering a conflict, that is, a subset of the
constraints that create this infeasibility. Depending on the solver,
this can also be called an IIS (irreducible inconsistent subsystem).

If supported by the solver, use [`compute_conflict!`](@ref) to trigger the
computation of a conflict. Once this process is finished, query the
[`MOI.ConflictStatus`](@ref) attribute to check if a conflict was found.

If found, copy the IIS to a new model using [`copy_conflict`](@ref), which you
can then print or write to a file for easier debugging:
```julia
julia> using JuMP

julia> import Gurobi

julia> model = Model(Gurobi.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0)
x

julia> @constraint(model, c1, x >= 2)
c1 : x ≥ 2.0

julia> @constraint(model, c2, x <= 1)
c2 : x ≤ 1.0

julia> optimize!(model)

julia> compute_conflict!(model)

julia> if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
           iis_model, _ = copy_conflict(model)
           print(iis_model)
       end
Feasibility
Subject to
 c1 : x ≥ 2.0
 c2 : x ≤ 1.0
```

If you need more control over the list of constraints that appear in the
conflict, iterate over the list of constraints and query the
[`MOI.ConstraintConflictStatus`](@ref) attribute:
```julia
julia> list_of_conflicting_constraints = ConstraintRef[]
ConstraintRef[]

julia> for (F, S) in list_of_constraint_types(model)
           for con in all_constraints(model, F, S)
               if get_attribute(con, MOI.ConstraintConflictStatus()) == MOI.IN_CONFLICT
                   push!(list_of_conflicting_constraints, con)
               end
           end
       end

julia> list_of_conflicting_constraints
2-element Vector{ConstraintRef}:
 c1 : x ≥ 2.0
 c2 : x ≤ 1.0
```

## Multiple solutions

Some solvers support returning multiple solutions. You can check how many
solutions are available to query using [`result_count`](@ref).

Functions for querying the solutions, for example, [`primal_status`](@ref),
[`dual_status`](@ref), [`value`](@ref), [`dual`](@ref), and [`solution_summary`](@ref)
all take an additional keyword argument `result` which can be used to specify
which result to return.

!!! warning
    Even if [`termination_status`](@ref) is [`OPTIMAL`](@ref), some of the
    returned solutions may be suboptimal. However, if the solver found at least
    one optimal solution, then `result = 1` will always return an optimal
    solution. Use [`objective_value`](@ref) to assess the quality of the
    remaining solutions.

```jldoctest; filter=[r"Solve time.+", r"Dual objective value.+"]
julia> using JuMP

julia> import MultiObjectiveAlgorithms as MOA

julia> import HiGHS

julia> model = Model(() -> MOA.Optimizer(HiGHS.Optimizer));

julia> set_attribute(model, MOA.Algorithm(), MOA.Dichotomy())

julia> set_silent(model)

julia> @variable(model, x1 >= 0)
x1

julia> @variable(model, 0 <= x2 <= 3)
x2

julia> @objective(model, Min, [3x1 + x2, -x1 - 2x2])
2-element Vector{AffExpr}:
 3 x1 + x2
 -x1 - 2 x2

julia> @constraint(model, 3x1 - x2 <= 6)
3 x1 - x2 ≤ 6

julia> optimize!(model)

julia> solution_summary(model; result = 1)
* Solver : MOA[algorithm=MultiObjectiveAlgorithms.Dichotomy, optimizer=HiGHS]

* Status
  Result count       : 3
  Termination status : OPTIMAL
  Message from the solver:
  "Solve complete. Found 3 solution(s)"

* Candidate solution (result #1)
  Primal status      : FEASIBLE_POINT
  Dual status        : NO_SOLUTION
  Objective value    : [0.00000e+00,0.00000e+00]
  Objective bound    : [0.00000e+00,-9.00000e+00]
  Relative gap       : Inf
  Dual objective value : -9.00000e+00

* Work counters
  Solve time (sec)   : 1.82720e-03
  Simplex iterations : 0
  Barrier iterations : 0
  Node count         : -1

julia> for i in 1:result_count(model)
           println("Solution $i")
           println("   x = ", value.([x1, x2]; result = i))
           println(" obj = ", objective_value(model; result = i))
       end
Solution 1
   x = [0.0, 0.0]
 obj = [0.0, 0.0]
Solution 2
   x = [0.0, 3.0]
 obj = [3.0, -6.0]
Solution 3
   x = [3.0, 3.0]
 obj = [12.0, -9.0]
```

!!! tip
    The [Multi-objective knapsack](@ref) tutorial provides more examples of
    querying multiple solutions.

## Checking feasibility of solutions

To check the feasibility of a primal solution, use
[`primal_feasibility_report`](@ref), which takes a `model`, a dictionary mapping
each variable to a primal solution value (defaults to the last solved solution),
and a tolerance `atol` (defaults to `0.0`).

The function returns a dictionary which maps the infeasible constraint
references to the distance between the primal value of the constraint and the
nearest point in the corresponding set. A point is classed as infeasible if the
distance is greater than the supplied tolerance `atol`.

```@meta
# Add a filter here because the output of the dictionary is not ordered, and
# changes in printing order will cause the doctest to fail.
```
```jldoctest feasibility; filter=[r"x.+?\=\> 0.1", r"c1.+? \=\> 0.01"]
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 1, Int);

julia> @variable(model, y);

julia> @constraint(model, c1, x + y <= 1.95);

julia> point = Dict(x => 1.9, y => 0.06);

julia> primal_feasibility_report(model, point)
Dict{Any, Float64} with 2 entries:
  x integer         => 0.1
  c1 : x + y ≤ 1.95 => 0.01

julia> primal_feasibility_report(model, point; atol = 0.02)
Dict{Any, Float64} with 1 entry:
  x integer => 0.1
```

If the point is feasible, an empty dictionary is returned:
```jldoctest feasibility
julia> primal_feasibility_report(model, Dict(x => 1.0, y => 0.0))
Dict{Any, Float64}()
```

To use the primal solution from a solve, omit the `point` argument:
```jldoctest feasibility
julia> optimize!(model)

julia> primal_feasibility_report(model; atol = 0.0)
Dict{Any, Float64}()
```
Calling [`primal_feasibility_report`](@ref) without the `point` argument is
useful when [`primal_status`](@ref) is [`FEASIBLE_POINT`](@ref) or
[`NEARLY_FEASIBLE_POINT`](@ref), and you want to assess the solution quality.

!!! warning
    To apply [`primal_feasibility_report`](@ref) to infeasible models, you must
    also provide a candidate point (solvers generally do not provide one). To
    diagnose the source of infeasibility, see [Conflicts](@ref).

Pass `skip_mising = true` to skip constraints which contain variables that are
not in `point`:
```jldoctest feasibility
julia> primal_feasibility_report(model, Dict(x => 2.1); skip_missing = true)
Dict{Any, Float64} with 1 entry:
  x integer => 0.1
```

You can also use the functional form, where the first argument is a function
that maps variables to their primal values:
```jldoctest feasibility
julia> optimize!(model)

julia> primal_feasibility_report(v -> value(v), model)
Dict{Any, Float64}()
```
