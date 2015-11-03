## Advanced model building concepts and techniques

In this chapter we discuss some more advanced model building concepts and
techniques that build on the material in the Getting Started chapter.

### Variable categories

JuMP supports N different variable *categories*:

* Continuous (`Cont`, the default)
* Integer (`Int`)
* Binary (`Bin`, `x ∈ {0,1}`)
* Semi-continuous (`SemiCont`, `x ∈ {0} ∪ [L,U]`)
* Semi-integer (`SemiInt`, `x ∈ {0} ∪ {L,…,U}`)
* Fixed value (`Fixed`, `x = L`)
* Positive semidefinite (`SDP`)

Solver support for these categories varies, and JuMP does not attempt to
automatically transform or otherwise handle categories that are not supported
(for example, JuMP will not do branch-and-bound for binary variables,
or reformulate semi-continuous variables).

##### Semi-continuous, semi-integer
Many of the open-source mixed-integer solvers do not support the semi-integer
and semi-continuous categories. However, it is possible to simulate these
categories with the use of an auxiliary binary variable:
```julia
# If the solver supports semi-continuous:
@defVar(m, 10 <= real_semicont <= 20, SemiCont)
# If the solver doesn't support semi-continuous
@defVar(m,  0 <= fake_semicont <= 20, Cont)
@defVar(m, auxiliary_var, Bin)
@addConstraint(m, fake_semicont >= 10*auxiliary_var)
@addConstraint(m, fake_semicont <= 20*auxiliary_var)
```

##### Positive Semidefinite
The `SDP` category differs from the other categories in that it applies only to
two-dimensional square matrices of variables, that are indexed from `1`:
```julia
# Not OK, not a 2D matrix
@defVar(m, three_d[1:5,1:5,1:5], SDP)
# Not OK, square 2D but not 1-based
@defVar(m, not_one[2:4,2:4], SDP)
# OK, square 2D, 1-based
@defVar(m, correct[1:5,1:5], SDP)
```

##### Fixed
The `Fixed` category is unlike the other categories in that it is not explicitly
passed to `@defVar`. Instead, we use `==` to indicate a variable is fixed to a
given value:
```julia
@defVar(m, one_fix == 5)
@defVar(m, many_fix[i=1:5] == i)
```
We can change the value the variable is fixed to using `setValue`.
The primary use of `Fixed` variables is to facilitate sensitivity analysis.
For example, we could solve a model for a variety of possible parameters:
```julia
knap = Model()
@defVar(knap, take[1:N], Bin)
@defVar(knap, capacity == 0)
@setObjective(knap, Max, dot(rewards,take))
@addConstraint(knap, dot(weights,take) <= capacity)
for C in 1:10
    setValue(capacity, C)
    solve(knap)
    @show C, getObjectiveValue(knap)
end
```

### Special ordered sets (SOS)

*Special ordered sets* are a different way to express constraints to integer solvers. They can be expressed in other ways, but if provided explicitly as an SOS then the solver *may* be able to solve the problem more efficiently. For more information about SOS, see the [`lpsolve` documentation](http://lpsolve.sourceforge.net/5.5/SOS.htm). There are two types of SOS constraint that are supported by JuMP. If a solver doesn't support SOS constraints then the user should manually reformulate using constraints and possibly auxiliary variables. Each variable in the set has a *weight*, and solvers expect the weights to be unique. The weights on the variables determine the ordering. If there is no natural weighting/ordering on the variables, an SOS constraint is probably unnecessary.

##### SOS Type 1
An SOS Type 1 constraint says that, of the variables contained in the set, at most one may be non-zero. Weights are provided as the coefficients of variables.
```julia
@defVar(m, x[1:5], Bin)
addSOS1(m, [i*x[i] for i in 1:5])
# Equivalent for binaries if SOS1 not supported
@addConstraint(m, sum(x) <= 1)
```

##### SOS Type 2
An SOS Type 2 constraint says that, of the variables contained in the set, no more than two adjacent (in the order) variables may be non-zero.
```julia
N = 10
@defVar(m, x[1:N], Bin)
addSOS2(m, [i*x[i] for i in 1:N])
# Equivalent for binaries if SOS2 not supported
for i in 1:N, j in i+2:N
    @addConstraint(m, x[i] + x[j] <= 1)
end
```

### Building up expressions manually

Some complex models can be difficult to express in a single expression, or have repeated components that appear in multiple places. For example, consider a stylized inventory control problem, with `T` time stages, where we have demand `d[t]` (data) and stock orders `x[t]` (decisons):
```julia
T = 10
demand = rand(50:250, T)
m = Model()
@defVar(m, x[1:T] >= 0)
```
In inventory control problems we often want to use the total inventory at each time period in multiple places. For example, we might put it in the objective function to penalize holding too much stock, or in a constraint to ensure we never have a negative amount of inventory. The total inventory at the end of time period `t` is the sum of stock ordered up to that time, less the sum of demand up to that time. To manually express the objective and constraint we just described, we could write our model as
```julia
holdingcost = rand(1:5, T)
@setObjective(m, Min, sum{holdingcost[t]*sum{x[s]-demand[s], s=1:t}, t=1:T})
for t in 1:T
    @addConstraint(m, sum{x[s]-demand[s], s=1:t} >= 0)
end
```
Alternatively, we could avoid repeating ourselves by creating an expression. There are three general ways to do this.
1. Define auxiliary variables and fix their values using equality constraints. How this is handled at the solver level varies - the solver may eliminate the auxiliaries in presolve, or leave it as is.
2. Use the `@defExpr` macro. The `@defExpr` macro creates expressions much as `@defVar` creates variables. Expressions can be indexed, and the expression itself can depend on the index. Unlike `@defExpr`, we do not need to pass the model as the first argument.
```julia
@defExpr(totalinv[t=1:T], sum{x[s]-demand[s], s=1:t})
@setObjective(m, Min, sum{holdingcost[t]*totalinv[t], t=1:T})
for t in 1:T
    @addConstraint(m, totalinv[t] >= 0)
end
```
3. Manually constructing the *affine expression*. Internally, JuMP stores a combination of coefficients and variables as an `AffExpr`. Users can manually create, construct, and combine `AffExpr`s. However, in general this should be avoided, as it can be far less efficient than use `@defExpr` for larger expressions. Operators like `+` and `*` can be used on expressions, but a more efficient method is to use `push!` and `append!`, e.g.,
```julia
@defVar(m, x[1:5])
foo = AffExpr()     # 0 (empty expression)
foo += x[1]         # 1 x[1]
bar = 5 * x[5]      # 5 x[5]
# append!(foo::AffExpr, bar::AffExpr)
append!(foo, bar)   # foo = 1 x[1] + 5 x[5]
# push!(foo::AffExpr, c::Number, x::Variable)
push!(foo, 3, x[3]) # foo = 1 x[1] + 3 x[3] + 5 x[5]
totalinv = zeros(AffExpr, T)
totalinv = x[1] - demand[1]
for t in 2:T
    totalinv[t] = total[t-1] + x[t] - demand[t]
end
```

### TODO/IDEAS
Performance
Array of SDPs
Passing variables into functions
