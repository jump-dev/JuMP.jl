# Copyright (c) 2021 Oscar Dowson and contributors                               #src
#                                                                                #src
# Permission is hereby granted, free of charge, to any person obtaining a copy   #src
# of this software and associated documentation files (the "Software"), to deal  #src
# in the Software without restriction, including without limitation the rights   #src
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      #src
# copies of the Software, and to permit persons to whom the Software is          #src
# furnished to do so, subject to the following conditions:                       #src
#                                                                                #src
# The above copyright notice and this permission notice shall be included in all #src
# copies or substantial portions of the Software.                                #src
#                                                                                #src
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #src
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       #src
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    #src
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         #src
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  #src
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE  #src
# SOFTWARE.                                                                      #src

# # Debugging

# Dealing with bugs is an unavoidable part of coding optimization models in
# any framework, including JuMP. Sources of bugs include not only generic coding
# errors (method errors, typos, off-by-one issues), but also semantic mistakes
# in the formulation of an optimization problem and the incorrect use of a
# solver.

# This tutorial explains some common sources of bugs and modeling issues that
# you might encounter when writing models in JuMP, and it suggests a variety of
# strategies to deal with them.

# !!! tip
#     This tutorial is more advanced than the other "Getting started" tutorials.
#     It's in the "Getting started" section to give you an early preview of how
#     to debug JuMP models. However, if you are new to JuMP, you may want to
#     briefly skim the tutorial, and come back to it once you have written a few
#     JuMP models.

using JuMP
import HiGHS

# ## Getting help

# Debugging can be a frustrating part of modeling, particularly if you're new to
# optimization and programming. If you're stuck, join the
# [community forum](https://discourse.julialang.org/c/domain/opt/13) to search
# for answers to commonly asked questions.

# Before asking a new question, make sure to read the post
# [Make it easier to help you](https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757),
# which contains a number of tips on how to ask a good question.

# Above all else, take time to simplify your code as much as possible. The
# fewer lines of code you can post that reproduces the same issue, the faster
# someone can answer your question.

# ## Debugging Julia code

# Read the [Debugging chapter](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html#chap21)
# in the book [ThinkJulia.jl](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html).
# It has a number of great tips and tricks for debugging Julia code.

# ## Solve failures

# When a solver experiences an issue that prevents it from finding an optimal
# solution (or proving that one does not exist), JuMP may return one of a number
# of [`termination_status`](@ref)es.

# For example, if the solver found a solution, but experienced numerical
# imprecision, it may return a status such as `ALMOST_OPTIMAL` or
# `ALMOST_LOCALLY_SOLVED` indicating that the problem was solved to a relaxed
# set of tolerances. Alternatively, the solver may return a problematic status
# such as `NUMERICAL_ERROR`, `SLOW_PROGRESS`, or `OTHER_ERROR`, indicating that
# it could not find a solution to the problem.

# Most solvers can experience numerical imprecision because they use
# [floating-point arithmetic](https://en.wikipedia.org/wiki/Floating-point_arithmetic)
# to perform operations such as addition, subtraction, and multiplication. These
# operations aren't exact, and small errors can accrue between the theoretical
# value and the value that the computer computes. For example:

0.1 * 3 == 0.3

# !!! tip
#     Read the [Guidelines for numerical issues](https://www.gurobi.com/documentation/9.5/refman/guidelines_for_numerical_i.html)
#     section of the Gurobi documentation, along with the
#     [Debugging numerical problems](https://yalmip.github.io/inside/debuggingnumerics/)
#     section of the YALMIP documentation.

# ### Common sources

# Common sources of solve failures are:
#
#  * Very large numbers and very small numbers as problem coefficients. Exactly
#    what "large" is depends on the solver and the problem, but in general,
#    values above 1e6 or smaller than 1e-6 cause problems.
#  * Nonlinear problems with functions that are not defined in parts of their
#    domain. For example, minimizing `log(x)` where `x >= 0` is undefined when
#    `x = 0` (a common starting value).

# ### Strategies

# Strategies to debug sources of solve failures include:
#
#  * Rescale variables in the problem and their associated coefficients to
#    make the magnitudes of all coefficients in the 1e-4 to 1e4 range. For
#    example, that might mean rescaling a variable from measuring distance in
#    centimeters to kilometers.
#  * Try a different solver. Some solvers might be more robust than others for a
#    particular problem.
#  * Read the documentation of your solver, and try settings that encourage
#    numerical robustness.
#  * Set bounds or add constraints so that all nonlinear functions are defined
#    across all of the feasible region. This particularly applies for functions
#    like `1 / x` and `log(x)` which are not defined for `x = 0`.

# ## Incorrect results

# Sometimes, you might find that the solver returns an "optimal" solution that
# is incorrect according to the model you are trying to solve (perhaps the
# solution is suboptimal, or it doesn't satisfy some of the constraints).

# Incorrect results can be hard to detect and debug, because the solver gives no
# hints that there is a problem. Indeed, the [`termination_status`](@ref) will
# likely be `OPTIMAL` and a solution will be available.

# ### Common sources

# Common sources of incorrect results are:
#
#  * A modeling error, so that your JuMP model does not match the formulation
#    you have on paper
#  * Not accounting for the tolerances that solvers use (for example, if `x` is
#    binary, a value like `x = 1.0000001` may still be considered feasible)
#  * A bug in JuMP or the solver.

# The probability of the issue being a bug in JuMP or the solver is much smaller
# than a modeling error. When in doubt, first assume there is a bug in your code
# before assuming that there is a bug in JuMP.

# ### Strategies

# Strategies to debug sources of incorrect results include:

#  * Print your JuMP model to see if it matches the formulation you have on
#    paper. Look out for incorrect signs `+` instead of `-`, and off-by-one
#    errors such as `x[t]` instead of `x[t-1]`.
#  * Check that you are not using exact comparisons like `value(x) == 1.0`;
#    always use `isapprox(value(x), 1.0; atol = 1e-6)` where you manually
#    specify the comparison tolerance.
#  * Try a different solver. If one solver succeeds where another doesn't this
#    is a sign that the problem is a numerical issue or a bug in the solver.

# ## Debugging an infeasible model

# A model is infeasible if there is no primal solution that satisfies all of
# the constraints. In general, an infeasible model means one of two things:
#
#  * Your problem really has no feasible solution
#  * There is a mistake in your model.

# ### Example

# A simple example of an infeasible model is:

model = Model(HiGHS.Optimizer);
set_silent(model)
@variable(model, x >= 0)
@objective(model, Max, 2x + 1)
@constraint(model, con, 2x - 1 <= -2)

# because the bound says that `x >= 0`, but we can rewrite the constraint to be
# `x <= -1/2`. When the problem is infeasible, JuMP may return one of a number
# of statuses. The most common is `INFEASIBLE`:

optimize!(model)
termination_status(model)

# Depending on the solver, you may also receive `INFEASIBLE_OR_UNBOUNDED` or
# `LOCALLY_INFEASIBLE`.

# A termination status of `INFEASIBLE_OR_UNBOUNDED` means that the solver could
# not prove if the solver was infeasible or unbounded, only that the model does
# not have a finite feasible optimal solution.

# Nonlinear optimizers such as Ipopt may return the status `LOCALLY_INFEASIBLE`.
# This does not mean that the solver _proved_ no feasible solution exists, only
# that it could not find one. If you know a primal feasible point, try providing
# it as a starting point using [`set_start_value`](@ref) and re-optimize.

# ### Common sources

# Common sources of infeasibility are:
#
#  * Incorrect units, for example, using a lower bound of megawatts and an upper
#    bound of kilowatts
#  * Using `+` instead of `-` in a constraint
#  * Off-by-one and related errors, for example, using `x[t]` instead of
#    `x[t-1]` in part of a constraint
#  * Otherwise invalid mathematical formulations

# ### Strategies

# Strategies to debug sources of infeasibility include:

#  * Iteratively comment out a constraint (or block of constraints) and re-solve
#    the problem. When you find a constraint that makes the problem infeasible
#    when added, check the constraint carefully for errors.
#  * If the problem is still infeasible with all constraints commented out,
#    check all variable bounds. Do they use the right data?
#  * If you have a known feasible solution, use [`primal_feasibility_report`](@ref)
#    to evaluate the constraints and check for violations. You'll probably find
#    that you have a typo in one of the constraints.
#  * Try a different solver. Sometimes, solvers have bugs, and they can
#    incorrectly report a problem as infeasible when it isn't. If you find such
#    a case where one solver reports the problem is infeasible and another can
#    find an optimal solution, please report it by opening an issue on the
#    GitHub repository of the solver that reports infeasibility.

# !!! tip
#     Some solvers also have specialized support for debugging sources of
#     infeasibility via an [irreducible infeasible subsystem](@ref Conflicts).
#     To see if your solver has support, try calling [`compute_conflict!`](@ref):
#     ```julia
#     julia> compute_conflict!(model)
#     ERROR: ArgumentError: The optimizer HiGHS.Optimizer does not support `compute_conflict!`
#     ```
#     In this case, HiGHS does not support computing conflicts, but other
#     solvers such as Gurobi and CPLEX do. If the solver does support computing
#     conflicts, read [Conflicts](@ref) for more details.

# ### Penalty relaxation

# Another strategy to debug sources of infeasibility is the
# [`relax_with_penalty!`](@ref) function.
#
# The penalty relaxation modifies constraints of the form ``f(x) \in S`` into
# ``f(x) + y - z \in S``, where ``y, z \ge 0``, and then it introduces a
# penalty term into the objective of ``a \times (y + z)`` (if minimizing, else
# ``-a``), where ``a`` is a penalty.

map = relax_with_penalty!(model)

# Here `map` is a dictionary which maps constraint indices to an affine
# expression representing ``(y + z)``.

# If we optimize the relaxed model, this time we get a feasible solution:

optimize!(model)
termination_status(model)

# Iterate over the contents of `map` to see which constraints are violated:

for (con, penalty) in map
    violation = value(penalty)
    if violation > 0
        println("Constraint `$(name(con))` is violated by $violation")
    end
end

# Once you find a violated constraint in the relaxed problem, take a look to see
# if there is a typo or other common mistake in that particular constraint.

# Consult the docstring [`relax_with_penalty!`](@ref) for information on how to
# modify the penalty cost term `a`, either for every constraint in the model or
# a particular subset of the constraints.

# When using [`relax_with_penalty!`](@ref), you should be aware that:
#
#  * Variable bounds and integrality restrictions are not relaxed. If the
#    problem is still infeasible after calling [`relax_with_penalty!`](@ref),
#    check the variable bounds.
#  * You cannot undo the penalty relaxation. If you need an unmodified model,
#    rebuild the problem, or call [`copy_model`](@ref) before calling
#    [`relax_with_penalty!`](@ref).

# ## Debugging an unbounded model

# A model is unbounded if there is no limit on how good the objective value can
# get. Most often, an unbounded model means that you have an error in your
# modeling, because all physical systems have limits. (You cannot make an
# infinite amount of profit.)

# ### Example

# A simple example of an unbounded model is:

model = Model(HiGHS.Optimizer);
set_silent(model)
@variable(model, x >= 0)
@objective(model, Max, 2x + 1)

# because we can increase `x` without limit, and the objective value `2x + 1`
# gets better as `x` increases.

# When the problem is unbounded, JuMP may return one of a number of statuses.
# The most common is `DUAL_INFEASIBLE`:

optimize!(model)
termination_status(model)

# Depending on the solver, you may also receive `INFEASIBLE_OR_UNBOUNDED` or
# an error code like `NORM_LIMIT`.

# ### Common sources

# Common sources of unboundedness are:
#
#  * Using `Max` instead of `Min`
#  * Omitting variable bounds, such as `0 <= x <= 1`
#  * Using `+` instead of `-` in a term of the objective function.

# ### Strategies

# Strategies to debug sources of unboundedness include:

#  * Double check whether you intended `Min` or `Max` in the [`@objective`](@ref)
#    line.
#  * Print the objective function with `print(objective_function(model))` and
#    verify that the value and sign of each coefficient is as you expect.
#  * Add large bounds to all variables that are free or have one-sided bounds,
#    then re-solve the problem. Because all variables are now bounded, the
#    problem will have a finite optimal solution. Look at the value of each
#    variable in the optimal solution to see if it is at one of the new bounds.
#    If it is, you either need to specify a better bound for that variable, or
#    there might be a mistake in the objective function associated with that
#    variable (for example, a `+` instead of a `-`).

# If there are too many variables to add bounds to, or there are too many terms
# to examine by hand, another strategy is to create a new variable with a large
# upper bound (if maximizing, lower bound if minimizing) and a constraint that
# the variable must be less-than or equal to the expression of the objective
# function. For example:

model = Model(HiGHS.Optimizer);
set_silent(model)
@variable(model, x >= 0)
## @objective(model, Max, 2x + 1)
@variable(model, objective <= 10_000)
@constraint(model, objective <= 2x + 1)
@objective(model, Max, objective)

# This new model has a finite optimal solution, so we can solve it and then look
# for variables with large positive or negative values in the optimal solution.

optimize!(model)
for var in all_variables(model)
    if var == objective
        continue
    end
    if abs(value(var)) > 1e3
        println("Variable `$(name(var))` may be unbounded")
    end
end
