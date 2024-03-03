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

# The purpose of this tutorial is to explain the similarities and differences of
# JuMP and [YALMIP](https://yalmip.github.io/) or [CVX](https://cvxr.com/cvx/).

# ## Namespaces

# Julia has namespaces, which MATLAB lacks. Therefore one needs to either use
# the command:

using JuMP

# in order bring all names exported by JuMP into scope, or:

import JuMP

# in order to merely make the JuMP package available. `import` requires
# prefixing everything you use from JuMP with `JuMP.`. In this tutorial we use
# the former.

# ## Starting a problem

# With YALMIP and CVX you have a single, implicit optimization problem you're
# working on. With JuMP you have to create one explicitly, and when you
# declare variables, constraints, or the objective function you have to
# specify to which problem they belong. You create a problem with the command:

model = Model()

# You can optionally also set the optimizer and options as an argument to the
# function `Model()`, but we'll do this later.

# ## Declaring variables

# In most cases there is a direct translation between variable declarations.
# The following table shows some common examples:

# |                          | JuMP                                              | YALMIP                                  | CVX                         |
# | ------------------------ | ------------------------------------------------- | --------------------------------------- | --------------------------- |
# | real variable            | `@variable(model, x)`                             | `x = sdpvar`                            | `variable x`                |
# | real vector              | `@variable(model, v[1:d])`                        | `v = sdpvar(d, 1)`                      | `variable v(d)`             |
# | real matrix              | `@variable(model, m[1:d, 1:d])`                   | `m = sdpvar(d,d,'full')`                | `variable m(d, d)`          |
# | complex matrix           | `@variable(model, m[1:d, 1:d] in ComplexPlane())` | `m = sdpvar(d,d,'full','complex')`      | `variable m(d,d) complex`   |
# | real symmetric matrix    | `@variable(model, m[1:d, 1:d], Symmetric)`        | `m = sdpvar(d)`                         | `variable m(d,d) symmetric` |
# | complex Hermitian matrix | `@variable(model, m[1:d, 1:d], Hermitian)`        | `m = sdpvar(d,d,'hermitian','complex')` | `variable m(d,d) hermitian` |

# A more interesting case is when you want to declare for example `n` real
# symmetric matrices. Both YALMIP and CVX allow you to put the matrices as the
# slices of a 3-dimensional array, via the commands `m = sdpvar(d,d,n)` and
# `variable m(d,d,n) symmetric`, respectively. With JuMP this is not possible.
# Instead, to achieve the same result one needs to declare a vector of `n`
# matrices: `m = [@variable(model, [1:d,1:d], Symmetric) for i=1:n]`.
# The analogous construct in MATLAB would be a cell array containing the
# optimization variables, which every discerning programmer avoids as cell
# arrays are rather slow. This is not a problem in Julia: a vector of matrices
# is almost as fast as a 3-dimensional array.

# ## Declaring constraints

# As in the case of variables, in most cases there is a direct translation
# between the syntaxes:

# |                     | JuMP                                                     | YALMIP               | CVX                                      |
# | ------------------- | -------------------------------------------------------- | -------------------- | ---------------------------------------- |
# | equality constraint | `@constraint(model, v == c)`                             | `v == c`             | `v == c`                                 |
# | nonnegative cone    | `@constraint(model, v >= 0)`                             | `v >= 0`             | `v >= 0`                                 |
# | real PSD cone       | `@constraint(model, m in PSDCone())`                     | `m >= 0`             | `m == semidefinite(length(m))`           |
# | complex PSD cone    | `@constraint(model, m in HermitianPSDCone())`            | `m >= 0`             | `m == hermitian_semidefinite(length(m))` |
# | second order cone   | `@constraint(model, [t; v] in SecondOrderCone())`        | `cone(v, t)`         | `{v, t} == lorentz(length(v))`           |
# | exponential cone    | `@constraint(model, [x, y, z] in MOI.ExponentialCone())` | `expcone([x, y, z])` | `{x, y, z} == exponential(1)`            |

# A subtlety appears when declaring equality constraints for matrices. In
# general, JuMP uses `@constraint(model, m .== c)`, with the dot meaning
# broadcasting in Julia, except when `m` is `Symmetric` or `Hermitian`: in this
# case `@constraint(model, m == c)` is allowed, and is much better, as JuMP is
# smart enough to not generate redundant constraints for the lower diagonal and
# the imaginary part of the diagonal (in the complex case). Both YALMIP and CVX
# are also smart enough to do this and the syntax is always just `m == c`.

# Experienced YALMIP users will probably be relieved to see that `>=` is only
# ever used to make a vector nonnegative, never to make a matrix positive
# semidefinite, as this ambiguity is reliable source of bugs.

# Like CVX, but unlike YALMIP, JuMP can also constrain variables upon creation:

# |                    | JuMP                                                 | CVX                                    |
# | ------------------ | ---------------------------------------------------- | -------------------------------------- |
# | nonnegative vector | `@variable(model, v[1:d] >= 0)`                      | `variable v(d) nonnegative`            |
# | real PSD matrix    | `@variable(model, m[1:d,1:d] in PSDCone())`          | `variable m(d,d) semidefinite`         |
# | complex PSD matrix | `@variable(model, m[1:d,1:d] in HermitianPSDCone())` | `variable m(d,d) complex semidefinite` |

# ## Setting the objective

# Like CVX, but unlike YALMIP, JuMP has a specific command for setting an
# objective function: `@objective(model, Min, obj)`. Here `obj` is any
# expression you want to optimize, and `Min` is an objective sense (the other
# possibility is `Max`).

# ## Setting solver and options

#  In order to set an optimizer with JuMP you can use at any point the command
# `set_optimizer(model, Hypatia.Optimizer)`, where "Hypatia" is an example
# solver. See the list of [Supported solvers](@ref) for more.

# To configure the solver options you use the command `set_attribute(model,
# "verbose", true)` after you set the optimizer, where as an example we were
# setting the option verbose to true.
#
# A crucial difference is that with JuMP you always have to explicitly choose a
# solver before optimizing. Both YALMIP and CVX allow you to leave it empty and
# will try to guess an appropriate solver for the problem.

# ## Optimizing

# Like YALMIP, but unlike CVX, with JuMP you need to explicitly start the
# optimization, with the command `optimize!(model)`. The exclamation mark here
# is a Julia-ism that means the function is modifying its argument, `model`.

# ## Extracting variables

# Like YALMIP, but unlike CVX, with JuMP you need to explicitly ask for the value
# of your variables after optimization is done, with the function call `value(x)`
# to obtain the value of variable `x`. A subtlety is that unlike YALMIP the
# function `value` is only defined for scalars, so for vectors and matrices you
# need to use Julia broadcasting: `value.(v)`. There's also a specialized function
# for extracting the value of the objective, `objective_value(model)`, which is
# useful if your objective doesn't have a convenient expression.

# ## Dual variables

# Like YALMIP and CVX, JuMP allows you to recover the dual variables. In order
# to do that, the simplest method is to name the constraint you're interested in,
# for example, `bob = @constraint(model, sum(v) == 1)` and then, after the
# optimzation is done, call `dual(bob)`.

# See [Duals of variable bounds](@ref) for more.

# ## Reformulating problems

# Perhaps the biggest difference between JuMP and YALMIP and CVX is how far the
# modeller is willing to go in reformulating the problems you give to it. CVX is
# happy to reformulate anything it can, even using approximations if your solver
# cannot handle the problem. YALMIP will only do exact reformulations, but is
# still fairly adventurous, being willing to reformulate a nonlinear objective
# in terms of conic constraints, for example. JuMP does no such thing: it only
# reformulates objectives into objectives, and constraints into constraints, and
# is fairly conservative at that. As a result, you might need to do some
# reformulations manually, for which a good guide is available [here](@ref conic_tips_and_tricks).

# An interesting possibility that JuMP offers is turning off reformulations
# altogether, which will reduce latency if they're in fact not needed, or
# produce an error if they are. To do that, set your optimizer with the option
# `add_bridges = false`, as for example, in `set_optimizer(model, Hypatia.Optimizer, add_bridges=false)`

# ## Vectorization

# In MATLAB it is absolutely essential to "vectorize" your code to obtain
# acceptable performance. This is because MATLAB is a very slow interpreted
# language, which sends your commands to extremely fast libraries. When you
# "vectorize" your code you're minimizing the MATLAB part of the work and sending
# it to the libraries instead. There's no such duality with Julia. Everything you
# write and most libraries you use will compile down to LLVM, so "vectorization"
# has no effect.

# For example, if you are writing a linear program in MATLAB and instead of the
# usual `constraints = [v >= 0]` you write:
# ```matlab
# for i = 1:n
#    constraints = [constraints, v(i) >= 0];
# end
# ```
# performance will be horrible. With Julia, on the other hand, there's hardly
# any difference between `@constraint(model, v >= 0)` and

for i in 1:n
    @constraint(model, v[i] >= 0)
end

# ## Rosetta stone

# To finish this tutorial, we show a complete example of the same optimization
# problem being solved with JuMP, YALMIP, and CVX. It is an SDP computing a
# simple lower bound on the random robustness of entanglement using the partial
# transposition criterion. The code is complete apart from the function that does
# partial transposition. With JuMP we use the function `partialtranspose` from
# [Convex.jl](https://jump.dev/Convex.jl/stable/). With both YALMIP and CVX we
# use the function `PartialTranspose` from [QETLAB](https://github.com/nathanieljohnston/QETLAB).

# ### JuMP

using JuMP
import Convex
import Hypatia
import LinearAlgebra

function random_state_pure(d)
    x = randn(ComplexF64, d)
    y = x * x'
    return LinearAlgebra.Hermitian(y / LinearAlgebra.tr(y))
end

function robustness_jump(d)
    model = Model()
    @variable(model, λ)
    ρ = random_state_pure(d^2)
	id = LinearAlgebra.Hermitian(LinearAlgebra.I(d^2))
	ρT = LinearAlgebra.Hermitian(Convex.partialtranspose(ρ, 1,[d ,d]))
    PPT = @constraint(model, ρT + λ * id in HermitianPSDCone())
    @objective(model, Min, λ)
    set_optimizer(model, Hypatia.Optimizer; add_bridges = false)
    set_attribute(model, "verbose", true)
    optimize!(model)
	WT = dual(PPT)
	return value(λ), LinearAlgebra.dot(WT, ρT)
end

# ### YALMIP

# ```matlab
# function rho = random_state_pure(d)
#     x = randn(d, 1) + 1i * randn(d, 1);
#     y = x * x';
#     rho = y / trace(y);
# end

# function robustness_yalmip(d)
#     rho = random_state_pure(d^2);
#     rhoT = PartialTranspose(rho, 1, [d d]);
#     lambda = sdpvar;
#     constraints = [(rhoT + lambda*eye(d^2) >= 0):'PPT'];
#     ops = sdpsettings(sdpsettings, 'verbose', 1, 'solver', 'sedumi');
#     optimize(constraints, lambda, ops);
#     WT = dual(constraints('PPT'));
#     value(lambda)
#     real(WT(:).'*rhoT(:))
# end
# ```

# ### CVX

# ```matlab
# function rho = random_state_pure(d)
#     x = randn(d, 1) + 1i * randn(d, 1);
#     y = x * x';
#     rho = y / trace(y);
# end

# function robustness_cvx(d)
#     rho = random_state_pure(d^2);
#     rhoT = PartialTranspose(rho, 1, [d d]);
#     cvx_begin
#         variable lambda
#         dual variable WT
#         WT : rhoT + lambda*eye(d^2) == hermitian_semidefinite(d^2)
#         minimise lambda
#     cvx_end
#     lambda
#     real(WT(:)'*rhoT(:))
# end
# ```
