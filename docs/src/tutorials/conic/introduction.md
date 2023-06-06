# Introduction

[Conic programs](https://en.wikipedia.org/wiki/Conic_optimization) are a class
of convex nonlinear optimization problems which use cones to represent the
nonlinearities. They have the form:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & f_0(x) \\
    & \;\;\text{s.t.} & f_j(x) \in \mathcal{S}_j & \;\; j = 1 \ldots m
\end{align}
```

Mixed-integer conic programs (MICPs) are extensions of conic programs in which
some (or all) of the decision variables take discrete values.

## How to choose a solver

JuMP supports a range of conic solvers, although support differs on what types
of cones each solver supports. In the list of [Supported solvers](@ref), "SOCP"
denotes solvers supporting second-order cones and "SDP" denotes solvers
supporting semidefinite cones. In addition, solvers such as SCS and Mosek have
support for the exponential cone. Moreover, due to the bridging system in
MathOptInterface, many of these solvers support a much wider range of exotic
cones than they natively support. Solvers supporting discrete variables start
with "(MI)" in the list of [Supported solvers](@ref).

!!! tip
    Duality plays a large role in solving conic optimization models. Depending
    on the solver, it can be more efficient to solve the dual instead of the
    primal. If performance is an issue, see the [Dualization](@ref) tutorial for
    more details.

## How these tutorials are structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * The following tutorials are worked examples that present a problem in words,
   then formulate it in mathematics, and then solve it in JuMP. This usually
   involves some sort of visualization of the solution. Start here if you are
   new to JuMP.
   * [Experiment design](@ref)
   * [Logistic regression](@ref)
 * The [Tips and tricks](@ref conic_tips_and_tricks) tutorial contains a
   number of helpful reformulations and tricks you can use when modeling
   conic programs. Look here if you are stuck trying to formulate a problem
   as a conic program.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
