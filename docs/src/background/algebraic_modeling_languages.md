```@meta
CurrentModule = JuMP
DocTestSetup = quote
    using JuMP, GLPK
end
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", r"MathOptInterface|MOI"]
```

# [Algebraic modeling languages](@id algebraic-modeling-language)

JuMP is an algebraic modeling language for mathematical optimization written in
the [Julia language](https://julialang.org). In this page, we explain what an
algebraic modeling language actually is.

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

    For example, [GLPK](https://www.gnu.org/software/glpk/) is a solver for
    linear programming (LP) and mixed integer programming (MIP) problems. It
    incorporates algorithms such as the simplex method and the interior-point
    method.

    JuMP currently supports a number of open-source and commercial solvers,
    which can be viewed in the [Supported-solvers](@ref) table.


Despite the textbook view of a linear program, you probably formulated problems
algebraically like so:
```math
\begin{aligned}
\max \;     & \sum\limits_{i = 1}^n c_i x_i                   \\
\text{s.t.} & \sum\limits_{i = 1}^n w_i x_i \le b             \\
            & x_i \ge 0 \quad \forall i = 1,\ldots,n          \\
            & x_i \in \mathbb{Z} \quad \forall i = 1,\ldots,n.
\end{aligned}
```
!!! info
    Do you recognize this formulation? It's the knapsack problem.

Users prefer to write problems in _algebraic form_ because it is more
convenient. For example, we used $\le b$, even though the standard form
only supported constraints of the form $Ax = b$.

We could convert our knapsack problem into the standard form by adding a new
slack variable $x_0$:
```math
\begin{aligned}
\max \;     & \sum\limits_{i = 1}^n c_i x_i            \\
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

Part 2 is less trivial than it might seem, because each solver has a unique
application programming interface (API) and data structure for representing
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
    variables and constraints, setting and getting parameters, and querying
    results. For more information on the MathOptInterface API, read the
    [documentation](@ref moi_documentation).

 2. An automatic rewriting system based on equivalent formulations of a
    constraint. For more information on this rewriting system, read the
    [LazyBridgeOptimizer](@ref) section of the manual, and our
    [paper on arXiv](https://arxiv.org/abs/2002.03447).

 3. Utilities for managing how and when models are copied to solvers. For more
    information on this, read the [CachingOptimizer](@ref) section of the
    manual.

## From user to solver

This section provides a brief summary of the steps that happen in order to
translate the model that the user writes into a model that the solver
understands.

### Step I: writing in algebraic form

JuMP provides the first part of an algebraic modeling language using the
[`@variable`](@ref), [`@objective`](@ref), and [`@constraint`](@ref) macros.

For example, here's how we write the knapsack problem in JuMP:
```jldoctest
julia> using JuMP, GLPK

julia> function algebraic_knapsack(c, w, b)
           n = length(c)
           model = Model(GLPK.Optimizer)
           @variable(model, x[1:n] >= 0, Int)
           @objective(model, Max, sum(c[i] * x[i] for i = 1:n))
           @constraint(model, sum(w[i] * x[i] for i = 1:n) <= b)
           optimize!(model)
           return value.(x)
       end
algebraic_knapsack (generic function with 1 method)

julia> algebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
This formulation is compact, and it closely matches the algebraic formulation of
the model we wrote out above.

### Step II: algebraic to functional

For the next step, JuMP's macros re-write the variables and constraints into a
functional form. Here's what the JuMP code looks like after this step:
```jldoctest
julia> using JuMP, GLPK

julia> function nonalgebraic_knapsack(c, w, b)
           n = length(c)
           model = Model(GLPK.Optimizer)
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
           set_objective(model, MAX_SENSE, obj)
           lhs = AffExpr(0.0)
           for i = 1:n
               add_to_expression!(lhs, w[i], x[i])
           end
           con = build_constraint(error, lhs, MOI.LessThan(b))
           add_constraint(model, con)
           optimize!(model)
           return value.(x)
       end
nonalgebraic_knapsack (generic function with 1 method)

julia> nonalgebraic_knapsack([1, 2], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```

Hopefully you agree that the macro version is much easier to read!

### Part III: JuMP to MathOptInterface

In the third step, JuMP converts the functional form of the problem, i.e.,
`nonalgebraic_knapsack`, into the MathOptInterface API:
```jldoctest
julia> using MathOptInterface, GLPK

julia> const MOI = MathOptInterface;

julia> function mathoptinterface_knapsack(optimizer, c, w, b)
           n = length(c)
           model = MOI.instantiate(optimizer)
           x = MOI.add_variables(model, n)
           for i in 1:n
               MOI.add_constraint(model, x[i], MOI.GreaterThan(0.0))
               MOI.add_constraint(model, x[i], MOI.Integer())
               MOI.set(model, MOI.VariableName(), x[i], "x[$i]")
           end
           MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
           obj = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0)
           MOI.set(model, MOI.ObjectiveFunction{typeof(obj)}(), obj)
           MOI.add_constraint(
               model,
               MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(w, x), 0.0),
               MOI.LessThan(b),
           )
           MOI.optimize!(model)
           return MOI.get.(model, MOI.VariablePrimal(), x)
       end
mathoptinterface_knapsack (generic function with 1 method)

julia> mathoptinterface_knapsack(GLPK.Optimizer, [1.0, 2.0], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
The code is becoming more verbose and looking less like the mathematical
formulation that we started with.

### Step IV: MathOptInterface to GLPK

As a final step, the [GLPK.jl](https://github.com/jump-dev/GLPK.jl) package
converts the MathOptInterface form, i.e., `mathoptinterface_knapsack`, into a
GLPK-specific API:
```jldoctest
julia> using GLPK

julia> function glpk_knapsack(c, w, b)
           n = length(c)
           model = glp_create_prob()
           glp_add_cols(model, n)
           for i in 1:n
               glp_set_col_bnds(model, i, GLP_LO, 0.0, GLP_DBL_MAX)
               glp_set_col_kind(model, i, GLP_IV)
               glp_set_col_name(model, i, "x[$i]")
           end
           glp_set_obj_dir(model, GLP_MAX)
           for i in 1:n
               glp_set_obj_coef(model, i, c[i])
           end
           glp_set_obj_coef(model, 0, 0.0)
           glp_add_rows(model, 1)
           glp_set_mat_row(
               model,
               1,
               length(w),
               Ref(Cint.(1:n), 0),
               Ref(w, 0),
           )
           glp_set_row_bnds(model, 1, GLP_UP, -GLP_DBL_MAX, b)
           simplex_options = glp_smcp()
           glp_init_smcp(simplex_options)
           simplex_options.msg_lev = GLP_MSG_ERR
           glp_simplex(model, simplex_options)
           options = glp_iocp()
           glp_init_iocp(options)
           options.msg_lev = GLP_MSG_ERR
           glp_intopt(model, options)
           x = glp_mip_col_val.(model, 1:n)
           glp_delete_prob(model)
           return x
       end
glpk_knapsack (generic function with 1 method)

julia> glpk_knapsack([1.0, 2.0], [0.5, 0.5], 1.25)
2-element Vector{Float64}:
 0.0
 2.0
```
We've now gone from a algebraic model that looked identical to the mathematical
model we started with, to a verbose function that uses GLPK-specific
functionality.

The difference between `algebraic_knapsack` and `glpk_knapsack` highlights the
benefit that algebraic modeling languages provide to users. Moreover, if we used
a different solver, the solver-specific function would be entirely different. A
key benefit of an algebraic modeling language is that you can change the solver
without needing to rewrite the model.
