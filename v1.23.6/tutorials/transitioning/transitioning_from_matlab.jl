# Copyright (c) 2024 Mateus Araújo and contributors                              #src
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

# # Transitioning from MATLAB

# [YALMIP](https://yalmip.github.io/) and [CVX](https://cvxr.com/cvx/) are two
# packages for mathematical optimization in [MATLAB®](https://mathworks.com/products/matlab.html).
# They are independently developed and are in no way affiliated with JuMP.

# The purpose of this tutorial is to help new users to JuMP who have previously
# used YALMIP or CVX by comparing and contrasting their different features.

# !!! tip
#     If you have not used Julia before, read the [Getting started with Julia](@ref)
#     tutorial.

# ## Namespaces

# Julia has namespaces, which MATLAB lacks. Therefore one needs to either use
# the command:

using JuMP

# in order bring all names exported by JuMP into scope, or:

import JuMP

# in order to merely make the JuMP package available. `import` requires
# prefixing everything you use from JuMP with `JuMP.`. In this tutorial we use
# the former.

# ## Models

# YALMIP and CVX have a single, implicit optimization model that you build by
# defining variables and constraints.

# In JuMP, we create an explicit model first, and then, when you declare
# variables, constraints, or the objective function, you specify to which model
# they are being added.

# Create a new JuMP model with the command:

model = Model()

# ## Variables

# In most cases there is a direct translation between variable declarations.
# The following table shows some common examples:

# | JuMP                                              | YALMIP                                  | CVX                         |
# | :------------------------------------------------ | :-------------------------------------- | :-------------------------- |
# | `@variable(model, x)`                             | `x = sdpvar`                            | `variable x`                |
# | `@variable(model, x, Int)`                        | `x = intvar`                            | `variable x integer`        |
# | `@variable(model, x, Bin)`                        | `x = binvar`                            | `variable x binary`         |
# | `@variable(model, v[1:d])`                        | `v = sdpvar(d, 1)`                      | `variable v(d)`             |
# | `@variable(model, m[1:d, 1:d])`                   | `m = sdpvar(d,d,'full')`                | `variable m(d, d)`          |
# | `@variable(model, m[1:d, 1:d] in ComplexPlane())` | `m = sdpvar(d,d,'full','complex')`      | `variable m(d,d) complex`   |
# | `@variable(model, m[1:d, 1:d], Symmetric)`        | `m = sdpvar(d)`                         | `variable m(d,d) symmetric` |
# | `@variable(model, m[1:d, 1:d], Hermitian)`        | `m = sdpvar(d,d,'hermitian','complex')` | `variable m(d,d) hermitian` |

# Like CVX, but unlike YALMIP, JuMP can also constrain variables upon creation:

# | JuMP                                                  | CVX                                    |
# | :---------------------------------------------------- | :------------------------------------- |
# | `@variable(model, v[1:d] >= 0)`                       | `variable v(d) nonnegative`            |
# | `@variable(model, m[1:d, 1:d], PSD)`                  | `variable m(d,d) semidefinite`         |
# | `@variable(model, m[1:d, 1:d] in PSDCone())`          | `variable m(d,d) semidefinite`         |
# | `@variable(model, m[1:d, 1:d] in HermitianPSDCone())` | `variable m(d,d) complex semidefinite` |

# JuMP can additionally set variable bounds, which may be handled more
# efficiently by a solver than an equivalent linear constraint. For example:

@variable(model, -1 <= x[i in 1:3] <= i)
upper_bound.(x)

# A more interesting case is when you want to declare, for example, `n` real
# symmetric matrices. Both YALMIP and CVX allow you to put the matrices as the
# slices of a 3-dimensional array, via the commands `m = sdpvar(d, d, n)` and
# `variable m(d, d, n) symmetric`, respectively. With JuMP this is not possible.
# Instead, to achieve the same result one needs to declare a vector of `n`
# matrices:

d, n = 3, 2
m = [@variable(model, [1:d, 1:d], Symmetric) for _ in 1:n]

#-

m[1]

#-

m[2]

# The analogous construct in MATLAB would be a cell array containing the
# optimization variables, which every discerning programmer avoids as cell
# arrays are rather slow. This is not a problem in Julia: a vector of matrices
# is almost as fast as a 3-dimensional array.

# ## [Constraints](@id matlab_constraints)

# As in the case of variables, in most cases there is a direct translation
# between the packages:

# | JuMP                                                     | YALMIP               | CVX                                      |
# | :------------------------------------------------------- | :------------------- | :--------------------------------------- |
# | `@constraint(model, v == c)`                             | `v == c`             | `v == c`                                 |
# | `@constraint(model, v >= 0)`                             | `v >= 0`             | `v >= 0`                                 |
# | `@constraint(model, m >= 0, PSDCone())`                  | `m >= 0`             | `m == semidefinite(length(m))`           |
# | `@constraint(model, m >= 0, HermitianPSDCone())`         | `m >= 0`             | `m == hermitian_semidefinite(length(m))` |
# | `@constraint(model, [t; v] in SecondOrderCone())`        | `cone(v, t)`         | `{v, t} == lorentz(length(v))`           |
# | `@constraint(model, [x, y, z] in MOI.ExponentialCone())` | `expcone([x, y, z])` | `{x, y, z} == exponential(1)`            |

# Like YALMIP and CVX, JuMP is smart enough to not generate redundant
# constraints when declaring equality constraints between `Symmetric` or
# `Hermitian` matrices. In these cases `@constraint(model, m == c)` will not
# generate constraints for the lower diagonal and the imaginary part of the
# diagonal (in the complex case).

# Experienced MATLAB users will probably be relieved to see that you must pass
# `PSDCone()` or `HermitianPSDCone()` to make a matrix positive semidefinite, because
# the `>=` ambiguity in YALMIP and CVX is common source of bugs.

# ## Setting the objective

# Like CVX, but unlike YALMIP, JuMP has a specific command for setting an
# objective function:

@objective(model, Min, sum(i * x[i] for i in 1:3))

# Here the third argument is any expression you want to optimize, and `Min` is
# an objective sense (the other possibility is `Max`).

# ## Setting solver and options

#  In order to set an optimizer with JuMP, do:

import Clarabel
set_optimizer(model, Clarabel.Optimizer)

# where "Clarabel" is an example solver. See the list of [Supported solvers](@ref)
# for other choices.

# To configure the solver options you use the command:

set_attribute(model, "verbose", true)

# where `verbose` is an option specific to Clarabel.

# A crucial difference is that with JuMP you must explicitly choose a solver
# before optimizing. Both YALMIP and CVX allow you to leave it empty and will
# try to guess an appropriate solver for the problem.

# ## Optimizing

# Like YALMIP, but unlike CVX, with JuMP you need to explicitly start the
# optimization, with the command:

optimize!(model)

# The exclamation mark here is a Julia-ism that means the function is modifying
# its argument, `model`.

# ## Querying solution status

# After the optimization is done, you should check for the solution status to
# see what solution (if any) the solver found.

# Like YALMIP and CVX, JuMP provides a solver-independent way to check it, via
# the command:

is_solved_and_feasible(model)

# If the return value is `false`, you should investigate with [`termination_status`](@ref),
# [`primal_status`](@ref), and [`raw_status`](@ref), See [Solutions](@ref jump_solutions)
# for more details on how to query and interpret solution statuses.

# ## Extracting variables

# Like YALMIP, but unlike CVX, with JuMP you need to explicitly ask for the value
# of your variables after optimization is done, with the function call `value(x)`
# to obtain the value of variable `x`.

value.(m[1][1, 1])

# A subtlety is that, unlike YALMIP, the function `value` is only defined for
# scalars. For vectors and matrices you need to use Julia broadcasting:
# `value.(v)`.

value.(m[1])

# There is also a specialized function for extracting the value of the objective,
# `objective_value(model)`, which is useful if your objective doesn't have a
# convenient expression.

objective_value(model)

# ## Dual variables

# Like YALMIP and CVX, JuMP allows you to recover the dual variables. In order
# to do that, the simplest method is to name the constraint you're interested in,
# for example, `@constraint(model, bob, sum(v) == 1)` and then, after the
# optimzation is done, call `dual(bob)`. See [Duality](@ref) for more details.

# ## Reformulating problems

# Perhaps the biggest difference between JuMP and YALMIP and CVX is how far the
# package is willing to go in reformulating the problems you give to it.

# CVX is happy to reformulate anything it can, even using approximations if your
# solver cannot handle the problem.

# YALMIP will only do exact reformulations, but is still fairly adventurous,
# for example, being willing to reformulate a nonlinear objective in terms of
# conic constraints.

# JuMP does no such thing: it only reformulates objectives into objectives, and
# constraints into constraints, and is fairly conservative at that. As a result,
# you might need to do some reformulations manually, for which a good guide is
# the [Modeling with cones](@ref) tutorial.

# ## Vectorization

# In MATLAB, it is absolutely essential to "vectorize" your code to obtain
# acceptable performance. This is because MATLAB is a slow interpreted
# language, which sends your commands to fast libraries. When you "vectorize"
# your code you are minimizing the MATLAB part of the work and sending it to the
# fast libraries instead.

# There's no such duality with Julia.

# Everything you write and most libraries you use will compile down to LLVM, so
# "vectorization" has no effect.

# For example, if you are writing a linear program in MATLAB and instead of the
# usual `constraints = [v >= 0]` you write:
# ```matlab
# for i = 1:n
#    constraints = [constraints, v(i) >= 0];
# end
# ```
# performance will be poor.

# With Julia, on the other hand, there is hardly any difference between
# ```julia
# @constraint(model, v >= 0)
# ```
#  and
# ```julia
# for i in 1:n
#     @constraint(model, v[i] >= 0)
# end
# ```

# ## Symmetric and Hermitian matrices

# Julia has specialized support for symmetric and Hermitian matrices in the
# `LinearAlgebra` package:

import LinearAlgebra

# If you have a matrix that is numerically symmetric:

x = [1 2; 2 3]

#-

LinearAlgebra.issymmetric(x)

# then you can wrap it in a `LinearAlgebra.Symmetric` matrix to tell Julia's
# type system that the matrix is symmetric.

LinearAlgebra.Symmetric(x)

# Using a `Symmetric` matrix lets Julia and JuMP use more efficient algorithms
# when they are working with symmetric matrices.

# If you have a matrix that is nearly but not exactly symmetric:

x = [1.0 2.0; 2.001 3.0]
LinearAlgebra.issymmetric(x)

# then you could, as you might do in MATLAB, make it numerically symmetric as
# follows:

x_sym = 0.5 * (x + x')

# In Julia, you can explicitly choose whether to use the lower or upper triangle
# of the matrix:

x_sym = LinearAlgebra.Symmetric(x, :L)

#-

x_sym = LinearAlgebra.Symmetric(x, :U)

# The same applies for Hermitian matrices, using `LinearAlgebra.Hermitian` and
# `LinearAlgebra.ishermitian`.

# ## Primal versus dual form

# When you translate some optimization problems from YALMIP or CVX to JuMP, you
# might be surprised to see it get much faster or much slower, even if you're
# using exactly the same solver. The most likely reason is that YALMIP will
# always interpret the problem as the dual form, whereas CVX and JuMP will try to
# interpret the problem in the form most appropriate to the solver. If the
# problem is more naturally formulated in the primal form it is likely that
# YALMIP's performance will suffer, or if JuMP gets it wrong, its performance will
# suffer. It might be worth trying both primal and dual forms if you're having
# trouble, which can be done automatically with the package [Dualization.jl](@ref).
#
# For an in-depth explanation of this issue, see the [Dualization](@ref) tutorial.
# ## Rosetta stone

# In this section, we show a complete example of the same optimization problem
# being solved with JuMP, YALMIP, and CVX. It is a semidefinite program that
# computes a lower bound on the random robustness of entanglement using the
# partial transposition criterion.

# The code is complete, apart from the function that does partial transposition.
# With both YALMIP and CVX we use the function `PartialTranspose` from
# [QETLAB](https://github.com/nathanieljohnston/QETLAB). With JuMP, we could use
# the function `Convex.partialtranspose` from [Convex.jl](https://jump.dev/Convex.jl/stable/),
# but we reproduce it here for simplicity:

function partial_transpose(x::AbstractMatrix, sys::Int, dims::Vector)
    @assert size(x, 1) == size(x, 2) == prod(dims)
    @assert 1 <= sys <= length(dims)
    n = length(dims)
    s = n - sys + 1
    p = collect(1:2n)
    p[s], p[n+s] = n + s, s
    r = reshape(x, (reverse(dims)..., reverse(dims)...))
    return reshape(permutedims(r, p), size(x))
end

# ### JuMP

# The JuMP code to solve this problem is:

using JuMP
import Clarabel
import LinearAlgebra

function random_state_pure(d)
    x = randn(Complex{Float64}, d)
    y = x * x'
    return LinearAlgebra.Hermitian(y / LinearAlgebra.tr(y))
end

function robustness_jump(d)
    rho = random_state_pure(d^2)
    id = LinearAlgebra.Hermitian(LinearAlgebra.I(d^2))
    rhoT = LinearAlgebra.Hermitian(partial_transpose(rho, 1, [d, d]))
    model = Model()
    @variable(model, λ)
    @constraint(model, PPT, rhoT + λ * id in HermitianPSDCone())
    @objective(model, Min, λ)
    set_optimizer(model, Clarabel.Optimizer)
    set_attribute(model, "verbose", true)
    optimize!(model)
    if is_solved_and_feasible(model)
        WT = dual(PPT)
        return value(λ), real(LinearAlgebra.dot(WT, rhoT))
    else
        return "Something went wrong: $(raw_status(model))"
    end
end

robustness_jump(3)

# ### YALMIP

# The corresponding YALMIP code is:

# ```matlab
# function robustness_yalmip(d)
#     rho = random_state_pure(d^2);
#     % PartialTranspose from https://github.com/nathanieljohnston/QETLAB
#     rhoT = PartialTranspose(rho, 1, [d d]);
#     lambda = sdpvar;
#     constraints = [(rhoT + lambda*eye(d^2) >= 0):'PPT'];
#     ops = sdpsettings(sdpsettings, 'verbose', 1, 'solver', 'sedumi');
#     sol = optimize(constraints, lambda, ops);
#     if sol.problem == 0
#         WT = dual(constraints('PPT'));
#         value(lambda)
#         real(WT(:)' * rhoT(:))
#     else
#         display(['Something went wrong: ', sol.info])
#     end
# end
#
# function rho = random_state_pure(d)
#     x = randn(d, 1) + 1i * randn(d, 1);
#     y = x * x';
#     rho = y / trace(y);
# end
# ```

# ### CVX

# The corresponding CVX code is:

# ```matlab
# function robustness_cvx(d)
#     rho = random_state_pure(d^2);
#     % PartialTranspose from https://github.com/nathanieljohnston/QETLAB
#     rhoT = PartialTranspose(rho, 1, [d d]);
#     cvx_begin
#         variable lambda
#         dual variable WT
#         WT : rhoT + lambda * eye(d^2) == hermitian_semidefinite(d^2)
#         minimise lambda
#     cvx_end
#     if strcmp(cvx_status, 'Solved')
#         lambda
#         real(WT(:)' * rhoT(:))
#     else
#         display('Something went wrong.')
#     end
# end
#
# function rho = random_state_pure(d)
#     x = randn(d, 1) + 1i * randn(d, 1);
#     y = x * x';
#     rho = y / trace(y);
# end
# ```
