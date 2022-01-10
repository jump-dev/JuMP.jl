# Introduction

[Linear programs (LPs)](https://en.wikipedia.org/wiki/Linear_programming) are a
fundamental class of optimization problems:
```math
\begin{align}
    & \min_{x \in \mathbb{R}^n} & \sum\limits_{i=1}^n c_i x_i
    \\
    & \;\;\text{s.t.} & l_j \le \sum\limits_{i=1}^n a_{ij} x_i \le u_j & j = 1 \ldots m
    & & l_i \le x_i \le u_i.
\end{align}
```
The most important thing to note is that all terms are of the form
`coefficient * variable`, and that there are no nonlinear terms or
multiplications between variables.

Mixed-integer linear programs (MILPs) are extensions of linear programs in which
some (or all) of the decision variables take discrete values.

## How to choose a solver

Almost all solvers support linear programs; look for "LP" in the list of
[Supported solvers](@ref). However, fewer solvers support mixed-integer linear
programs. Solvers supporting discrete variables start with "(MI)" in the list of
[Supported solvers](@ref).

## How this section is structured

Having a high-level overview of how this part of the documentation is structured
will help you know where to look for certain things.

 * Tutorial-style worked examples present a problem in words, then formulate it
   in mathematics, and then solve it in JuMP. This usually involves some sort of
   visualization of the solution. Start here if you are new to JuMP.
   * [The diet problem](@ref)
   * [The cannery problem](@ref)
   * [The facility location problem](@ref)
   * [Financial modeling problems](@ref)
   * [Network flow problems](@ref)
   * [N-Queens](@ref)
   * [Sudoku](@ref)
 * The [Tips and tricks](@ref linear_tips_and_tricks) tutorial contains a number
   of helpful reformulations and tricks you can use when modeling linear
   programs. Look here if you are stuck trying to formulate a problem as a
   linear program.
 * The [Callbacks](@ref callbacks_tutorial) tutorial explains how to write a
   variety of solver-independent callbacks. Look here if you want to write a
   callback.
 * The remaining tutorials are less verbose and styled in the form of short code
   examples. These tutorials have less explanation, but may contain useful
   code snippets, particularly if they are similar to a problem you are trying
   to solve.
