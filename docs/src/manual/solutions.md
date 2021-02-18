```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# Solutions

This section of the manual describes how to access a solved solution to a
problem. It uses the following model as an example:
```jldoctest solutions
model = Model(GLPK.Optimizer)
@variable(model, x >= 0)
@variable(model, y[[:a, :b]] <= 1)
@objective(model, Max, -12x - 20y[:a])
@expression(model, my_expr, 6x + 8y[:a])
@constraint(model, my_expr >= 100)
@constraint(model, c1, 7x + 12y[:a] >= 120)
optimize!(model)
print(model)

# output

Max -12 x - 20 y[a]
Subject to
 6 x + 8 y[a] ≥ 100.0
 c1 : 7 x + 12 y[a] ≥ 120.0
 x ≥ 0.0
 y[a] ≤ 1.0
 y[b] ≤ 1.0
```

## Why did the solver stop?

Use[`termination_status`](@ref) to understand why the solver stopped.

```jldoctest solutions
julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1
```

The [`MOI.TerminationStatusCode`](@ref) enum describes the full list of statuses
that could be returned.

Common return values include `MOI.OPTIMAL`, `MOI.LOCALLY_SOLVED`,
`MOI.INFEASIBLE`, `MOI.DUAL_INFEASIBLE`, and `MOI.TIME_LIMIT`.

!!! info
    A return status of `MOI.OPTIMAL` means the solver found (and proved) a
    globally optimal solution. A return status of `MOI.LOCALLY_SOLVED` means the
    solver found a locally optimal solution (which may also be globally
    optimal, but it could not prove so).

!!! warning
    A return status of `MOI.DUAL_INFEASIBLE` does not guarantee that the primal
    is unbounded. When the dual is infeasible, the primal is unbounded if
    there exists a feasible primal solution.

Use [`raw_status`](@ref) to get a solver-specific string explaining why the
optimization stopped:
```jldoctest solutions
julia> raw_status(model)
"Solution is optimal"
```

## Primal solutions

### Primal solution status

Use [`primal_status`](@ref) to return an [`MOI.ResultStatusCode`](@ref) enum
describing the status of the primal solution.
```jldoctest solutions
julia> primal_status(model)
FEASIBLE_POINT::ResultStatusCode = 1
```
Other common returns are `MOI.NO_SOLUTION`, and `MOI.INFEASIBILITY_CERTIFICATE`.
The first means that the solver doesn't have a solution to return, and the
second means that the solution has a certificate of dual infeasbility (a primal
unbounded ray).

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

julia> objective_bound(model)  # GLPK only implements objective bound for MIPs
Inf

julia> dual_objective_value(model)
-205.1428571428571
```

### Primal solution values

If the solver has a primal solution to return, use [`value`](@ref) to access it:
```julia solutions
julia> value(x)
15.428571428571429
```

Broadcast [`value]`(@ref) over containers:
```julia solutions
julia> value.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
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
Other common returns are `MOI.NO_SOLUTION`, and `MOI.INFEASIBILITY_CERTIFICATE`.
The first means that the solver doesn't have a solution to return, and the
second means that the solution has a certificate of primal infeasbility (a dual
unbounded ray).

You can also use [`has_duals`](@ref), which returns `true` if there is a
solution that can be queried, and `false` otherwise.
```jldoctest solutions
julia> has_duals(model)
true
```

### Dual solution values

If the solver has a dual solution to return, use [`dual`](@ref) to access it:
```julia solutions
julia> dual(c1)
1.7142857142857142
```

Query the duals of variable bounds using [`LowerBoundRef`](@ref),
[`UpperBoundRef`](@ref), and [`FixRef`](@ref):
```julia solutions
julia> dual(LowerBoundRef(x))
0.0

julia> dual.(UpperBoundRef.(y))
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
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

```julia solutions
julia> shadow_price(c1)
1.7142857142857142

julia> reduced_cost(x)
0.0

julia> reduced_cost.(y)
1-dimensional DenseAxisArray{Float64,1,...} with index sets:
    Dimension 1, Symbol[:a, :b]
And data, a 2-element Array{Float64,1}:
  0.5714285714285694
 -0.0
```

## Recommended workflow

The following is a recommended workflow for solving a model and querying the
solution is something like the following:
```jldoctest solutions
if termination_status(model) == MOI.OPTIMAL
    println("Solution is optimal")
elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
    println("Solution is suboptimal")
else
    error("The model was not solved correctly.")
end
println("  objective value = ", objective_value(model))
if primal_status(model) == MOI.FEASIBLE_POINT
    println("  primal solution: x = ", value(x))
end
if dual_status(model) == MOI.FEASIBLE_POINT
    println("  dual solution: c1 = ", dual(c1))
end

# output

Solution is optimal
  objective value = -205.14285714285714
  primal solution: x = 15.428571428571429
  dual solution: c1 = 1.7142857142857142
```

!!! warning
    Querying solution information after modifying a solved model is undefined
    behavior, and solvers may throw an error or return incorrect results.
    Modifications include adding, deleting, or modifying any variable,
    objective, or constraint. Instead of modify-then-query, query the results
    first, then modify the problem. For example:
    ```julia
    model = Model(GLPK.Optimizer)
    @variable(model, x >= 0)
    optimize!(model)
    # Bad:
    set_lower_bound(x, 1)
    @show value(x)
    # Good:
    x_val = value(x)
    set_lower_bound(x, 1)
    @show x_val
    ```

```@meta
# TODO: How to accurately measure the solve time.
```

## Accessing attributes

[MathOptInterface](https://jump.dev/MathOptInterface.jl/v0.9.10/) defines a large
number of model attributes that can be queried. Some attributes can be directly
accessed by getter functions. These include:
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
solver interface must implement `MOI.ConstraintBasisStatus`.

To give a simple example, we could analyze the sensitivity of the optimal
solution to the following (non-degenerate) LP problem:

```jldoctest solutions_sensitivity
model = Model(GLPK.Optimizer)
@variable(model, x[1:2])
set_lower_bound(x[2], -0.5)
set_upper_bound(x[2], 0.5)
@constraint(model, c1, x[1] + x[2] <= 1)
@constraint(model, c2, x[1] - x[2] <= 1)
@objective(model, Max, x[1])
print(model)

# output

Max x[1]
Subject to
 c1 : x[1] + x[2] ≤ 1.0
 c2 : x[1] - x[2] ≤ 1.0
 x[2] ≥ -0.5
 x[2] ≤ 0.5
```

To analyze the sensitivity of the problem we could check the allowed
perturbation ranges of, e.g., the cost coefficients and the right-hand side
coefficient of the constraint `c1` as follows:

```jldoctest solutions_sensitivity
julia> optimize!(model)

julia> value.(x)
2-element Array{Float64,1}:
 1.0
 0.0

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
ranges might not be as intuitive and seem too narrow. E.g., a larger cost
coefficient perturbation might not invalidate the optimality of the current
primal solution. Moreover, if a problem is degenerate, due to finite precision,
it can happen that, e.g., a perturbation seems to invalidate a basis even though
it doesn't (again providing too narrow ranges). To prevent this, increase the
`atol` keyword argument to [`lp_sensitivity_report`](@ref). Note that this might
make the ranges too wide for numerically challenging instances. Thus, do not
blindly trust these ranges, especially not for highly degenerate or numerically
unstable instances.

## Conflicts

When the model you input is infeasible, some solvers can help you find the
cause of this infeasibility by offering a conflict, i.e., a subset of the
constraints that create this infeasibility. Depending on the solver,
this can also be called an IIS (irreducible inconsistent subsystem).

The function [`compute_conflict!`](@ref) is used to trigger the computation of
a conflict. Once this process is finished, the attribute
[`MOI.ConflictStatus`](@ref) returns a [`MOI.ConflictStatusCode`](@ref).

If there is a conflict, you can query from each constraint whether it
participates in the conflict or not using the attribute
[`MOI.ConstraintConflictStatus`](@ref), which returns a
[`MOI.ConflictParticipationStatusCode`](@ref).

For instance, this is how you can use this functionality:

```julia
using JuMP
model = Model() # You must use a solver that supports conflict refining/IIS computation, like CPLEX or Gurobi
@variable(model, x >= 0)
@constraint(model, c1, x >= 2)
@constraint(model, c2, x <= 1)
optimize!(model)

# termination_status(model) will likely be MOI.INFEASIBLE,
# depending on the solver

compute_conflict!(model)
if MOI.get(model, MOI.ConflictStatus()) != MOI.CONFLICT_FOUND
    error("No conflict could be found for an infeasible model.")
end

# Both constraints should participate in the conflict.
MOI.get(model, MOI.ConstraintConflictStatus(), c1)
MOI.get(model, MOI.ConstraintConflictStatus(), c2)
```

## Multiple solutions

Some solvers support returning multiple solutions. You can check how many
solutions are available to query using [`result_count`](@ref).

Functions for querying the solutions, e.g., [`primal_status`](@ref) and
[`value`](@ref), all take an additional keyword argument `result` which can be
used to specify which result to return.

!!! warning
    Even if [`termination_status`](@ref) is `MOI.OPTIMAL`, some of the returned
    solutions may be suboptimal! However, if the solver found at least one
    optimal solution, then `result = 1` will always return an optimal solution.
    Use [`objective_value`](@ref) to assess the quality of the remaining
    solutions.

```julia
using JuMP
model = Model()
@variable(model, x[1:10] >= 0)
# ... other constraints ...
optimize!(model)

if termination_status(model) != MOI.OPTIMAL
    error("The model was not solved correctly.")
end

an_optimal_solution = value.(x; result = 1)
optimal_objective = objective_value(model; result = 1)
for i in 2:result_count(model)
    @assert has_values(model; result = i)
    println("Solution $(i) = ", value.(x; result = i))
    obj = objective_value(model; result = i)
    println("Objective $(i) = ", obj)
    if isapprox(obj, optimal_objective; atol = 1e-8)
        print("Solution $(i) is also optimal!")
    end
end
```
