## Working with solutions

After you've modeled your problem, and passed it to your solver, the next step is to understand what the solver came up with. Alternatively, we may want to provide a starting solution to guide the solver. In this section we will discuss the various scenarios, the information that is available, and how to get it.

### Solver statuses

The reported status of the solver is returned from `solve`. A solver can return any value, but the following statuses are the most common and are encouraged:

 * `:Optimal`
 * `:Infeasible`
 * `:Unbounded`
 * `:UserLimit` (e.g., iteration limit, time limit)
 * `:Error`

The precise meaning of each status will vary between solvers. For example, for integer optimization problems it is common to terminate when the bound gap falls below a certain value. While this is technically "termination due to a user limit", it will normally be described as `:Optimal`. See your solver documentation for further clarification.

### Primal and dual solutions

The primal and dual solutions can be accessed with the `getValue` and `getDual` functions. However, the exact behavior depends on the solver status and the solver itself. Here we will describe three cases: optimal, infeasible, and unbounded.

##### Optimal

The simplest case is when the solver reports an optimal solution. To obtain the primal solution, we can use `getValue` on a variable. To obtain the duals for constraints (also known as "shadow prices"), we need to retain a reference to the constraints when we construct them. Finally, to get reduced costs (or equivalently, the dual for variables bounds), we can call `getDual` on a variable.
```julia
m = Model()
@defVar(m, x)
@defVar(m, y <= 3)
@defVar(m, -1 <= z <= 1)
@setObjective(m, Max, 3x + 7y - 5z)
π = @addConstraint(m, x <= 2)
solve(m)
@show getValue(x) # =  2 (primal solution)
@show getDual(π)  # =  3 (dual on constraint)
@show getDual(y)  # =  7 (dual on bound)
@show getDual(z)  # = -5 (dual on lower bound)
```
To understand all four numbers:
* Variable `x` is `2` (`getValue(x)`) because we are maximizing `3x`, and the constraint `π` limits the value to be no more than `2`. If we increase the right-hand-side of the constraint by 1, the objective function value would increase by `(3-2)*3=3` (`getDual(π)`).
* Variable `y` is at its upper bound, `3`, so `getDual(y)` is `7` because if we increased the bound by `1` (to `4`) then the objective function value would change by `(4 - 3)*7 = 7`.
* Variable `z` is at its lower bound, `-1`, so `getDual(z)` is `-5` because if we increased the bound by `1` (to `0`), then the objective function would change by `(0 - -1)*-5 = -5`.

If we call `getValue` or `getDual` on an indexed variable or constraint, we will obtain an object that is indexed the same way. If the variables or constraints are indexed "simply" (from `1` to some other integer, in steps of 1), then a normal Julia array will be returned that can be used easily in further computations.
```julia
m = Model()
@defVar(m, x[1:10] <= 1)    # "Simply" indexed
@defVar(m, y[2:2:8] <= 1)   # Not indexed from 1
@setObjective(m, Max, sum(x) + sum(y))
@addConstraint(m, xcon[i=1:10], x[i] <= i/10)
@addConstraint(m, ycon[j=2:2:8], y[j] <= j/10)
solve(m)
@show isa(getValue(x),   Vector{Float64})   # true
@show isa(getDual(xcon), Vector{Float64})   # true
@show isa(getValue(y),   Vector{Float64})   # false
@show isa(getDual(ycon), Vector{Float64})   # false
```

##### Infeasible

If the problem is `:Infeasible`, this normally means that there is no primal feasible solution. However, `getDual` can still return an *infeasibility ray*, or *Farkas proof*, assuming the solver provides one. If we add this ray to a dual-feasible solution, then the dual objective function value will improve.
```julia
primal = Model()
@defVar(primal, x >= 0)
@setObjective(primal, Min, x)
con = @addConstraint(primal, x == -1)
solve(primal)
@show getDual(con)  # -1
# Dual problem:
#   Max -π subject to π <= 1
# Feasible solutions:
#   1 + (-1)*k = 1 - k, for all k >= 0
```

##### Unbounded

If the problem is `:Unbounded`, we can obtain an *unbounded ray* if the solver provides one. Adding this ray to a primal-feasible solution will improve the primal objective function value.
```julia
m = Model()
@defVar(m, x >= 0)
@defVar(m, y >= 0)
@setObjective(m, Max, x + y)
@addConstraint(m, x == 2y)
solve(m)
@show getValue(x)   # 2.0
@show getValue(y)   # 1.0
# Feasible solutions:
#  (x + 2k, y + k) for all k >= 0
# Objective function values:
#  x + y + 3k
```
