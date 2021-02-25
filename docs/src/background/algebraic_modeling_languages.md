```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Algebraic modeling languages](@id algebraic-modeling-language)

## What is an algebraic modeling language?

If you have taken a class in mixed-integer linear programming, you will have
seen a formulation like:
```math
\begin{aligned}
\min \;     & c^\top x \\
\text{s.t.} & A x = b  \\
            & x \ge 0  \\
            & x_i \in \mathbb{Z}, \quad \forall i \in \mathcal{I}
\end{aligned}
```
where `c`, `A`, and `b` are appropriately sized vectors and matrices of data,
and $\mathcal{I}$ denotes the set of variables that are integer.

Solvers expect problems in a _standard form_ like this because it limits the
types of constraints that they need to consider. This makes writing a solver
much easier.

!!! info "What is a solver?"
    A solver is a software package that computes solutions to one or more
    classes of problems.

    For example, GLPK is a solver for linear programming (LP) and mixed integer
    programming (MIP) problems. It incorporates algorithms such as the simplex
    method and the interior-point method.

    JuMP currently supports a number of open-source and commercial solvers,
    which can be viewed in the [Supported-solvers](@ref) table.


However, you probably formulated problems algebraically like so:
```math
\begin{aligned}
\min \;     & \sum\limits_{i = 1}^n c_i x_i                   \\
\text{s.t.} & \sum\limits_{i = 1}^n w_i x_i \le b             \\
            & x_i \ge 0 \quad \forall i = 1,\ldots,n          \\
            & x_i \in \mathbb{Z} \quad \forall i = 1,\ldots,n.
\end{aligned}
```
!!! info
    Do you recognize this formulation? It's the knapsack problem.

Users prefer to write problems in _algebraic form_ because it is more
convenient. For example, we just used $\le b$, even though the standard form
only supported constraints of the form $Ax = b$.

We could convert our knapsack problem into the standard form by adding a new
slack variable $x_0$ like so:
```math
\begin{aligned}
\min \;     & \sum\limits_{i = 1}^n c_i x_i            \\
\text{s.t.} & x_0 + \sum\limits_{i = 1}^n w_i x_i = b  \\
            & x_i \ge 0 \quad \forall i = 0,\ldots,n   \\
            & x_i \in \mathbb{Z} \quad \forall i = 1,\ldots,n.
\end{aligned}
```
However, as models get more complicated, this manual conversion becomes more and
more error-prone.

An algebraic modeling language is a tool that simplifies the translation between
the algebraic form of the modeler, and the standard form of the solver.

Each algebraic modeling language has two main parts:

 1. A domain specific language for the user to write down problems in algebraic
    form.
 2. A converter from the algebraic form into a standard form supported by the
    solver (and back again).

## Part I: writing in algebraic form

JuMP provides the first part of an algebraic modeling language using the
[`@variable`](@ref), [`@objective`](@ref), and [`@constraint`](@ref) macros.

For example, here's how we would write the knapsack problem in JuMP:
```jldoctest
julia> function algebraic_knapsack(c, w, b)
           n = length(c)
           model = Model()
           @variable(model, x[1:n] >= 0, Int)
           @objective(model, Min, sum(c[i] * x[i] for i = 1:n))
           @constraint(model, sum(w[i] * x[i] for i = 1:n) <= b)
           return print(model)
       end
algebraic_knapsack (generic function with 1 method)

julia> algebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
Min x[1] + 2 x[2]
Subject to
 0.5 x[1] + 0.5 x[2] ≤ 1.25
 x[1] ≥ 0.0
 x[2] ≥ 0.0
 x[1] integer
 x[2] integer
```
This formulation is compact, and it closely matches the algebraic formulation of
the model we wrote out above.

Here's what the JuMP code would look like if we didn't use macros:
```jldoctest
julia> function nonalgebraic_knapsack(c, w, b)
           n = length(c)
           model = Model()
           x = [VariableRef(model) for i = 1:n]
           for i = 1:n
               set_lower_bound(x[i], 0)
               set_integer(x[i])
               set_name(x[i], "x[$i]")
           end
           obj = AffExpr(0.0)
           for i = 1:n
               add_to_expression!(obj, c[i], x[i])
           end
           set_objective(model, MOI.MIN_SENSE, obj)
           lhs = AffExpr(0.0)
           for i = 1:n
               add_to_expression!(lhs, w[i], x[i])
           end
           con = build_constraint(error, lhs, MOI.LessThan(b))
           add_constraint(model, con)
           return print(model)
       end
nonalgebraic_knapsack (generic function with 1 method)

julia> nonalgebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
Min x[1] + 2 x[2]
Subject to
 0.5 x[1] + 0.5 x[2] ≤ 1.25
 x[1] ≥ 0.0
 x[2] ≥ 0.0
 x[1] integer
 x[2] integer
```

Hopefully you agree that the macro version is much easier to read!

## Part II: talking to solvers

Now that we have the algebraic problem from the user, we need a way of
communicating the problem to the solver, and a way of returning the solution
from the solver back to the user.

This is less trivial than it might seem, because each solver has a unique
application programming interface (API) and data structures for representing
optimization models and obtaining results.

JuMP uses the
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)
package to abstract these differences between solvers.

### What is MathOptInterface?

MathOptInterface (MOI) is an abstraction layer designed to provide an
interface to mathematical optimization solvers so that users do not need to
understand multiple solver-specific APIs. MOI can be used directly, or through
a higher-level modeling interface like JuMP.

There are three main parts to MathOptInterface:

 1. A solver-independent API that abstracts concepts such as adding and deleting
    variables and constraints, settting and getting parameters, and querying
    results. For more information on the MathOptInterface API, read the
    [documentation]([https://jump.dev/MathOptInterface.jl/stable/])

 2. An automatic rewriting system based on equivalent formulations of a
    constraint. For more information on this rewriting system, read the
    [LazyBridgeOptimizer](@ref) section of the manual, and our
    [paper on arXiv](https://arxiv.org/abs/2002.03447).

 3. Utilities for managing how and when models are copied to solvers. For more
    information on this, read the [CachingOptimizer](@ref) section of the
    manual.
