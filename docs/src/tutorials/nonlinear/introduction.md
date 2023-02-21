# Introduction

[Nonlinear programs (NLPs)](https://en.wikipedia.org/wiki/Nonlinear_programming)
are a class of optimization problems in which some of the constraints or the
objective function are nonlinear:
```math
\begin{align}
    \min_{x \in \mathbb{R}^n} & f_0(x) \\
    \;\;\text{s.t.} & l_j \le f_j(x) \le u_j & j = 1 \ldots m \\
    & l_i \le x_i \le u_i & i = 1 \ldots n.
\end{align}
```

Mixed-integer nonlinear linear programs (MINLPs) are extensions of nonlinear
programs in which some (or all) of the decision variables take discrete values.

## How to choose a solver

JuMP supports a range of nonlinear solvers; look for "NLP" in the list
of [Supported solvers](@ref). However, very few solvers support mixed-integer
nonlinear linear programs. Solvers supporting discrete variables start with
"(MI)" in the list of [Supported solvers](@ref).

If the only nonlinearities in your model are quadratic terms (that is,
multiplication between two decision variables), you can also use second-order
cone solvers, which are indicated by "SOCP." In most cases, these solvers are
restricted to convex quadratic problems and will error if you pass a nonconvex
quadratic function; however, Gurobi has the ability to solve nonconvex quadratic
terms.

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The following tutorials are worked examples that present a problem in words,
   then formulate it in mathematics, and then solve it in JuMP. This usually
   involves some sort of visualization of the solution. Start here if you are
   new to JuMP.
   * [Rocket Control](@ref)
   * [Optimal control for a Space Shuttle reentry trajectory](@ref)
   * [Portfolio optimization](@ref)
 * The [Tips and tricks](@ref nonlinear_tips_and_tricks) tutorial contains a
   number of helpful reformulations and tricks you can use when modeling
   nonlinear programs. Look here if you are stuck trying to formulate a problem
   as a nonlinear program.
 * The [Computing Hessians](@ref) is an advanced tutorial which explains how to
   compute the Hessian of the Lagrangian of a nonlinear program. This is useful
   only in particular cases.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
