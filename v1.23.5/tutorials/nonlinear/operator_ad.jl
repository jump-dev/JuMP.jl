# Copyright (c) 2024 Oscar Dowson and contributors                               #src
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

# # Automatic differentiation of user-defined operators

# The purpose of this tutorial is to demonstrate how to apply automatic
# differentiation to [User-defined operators](@ref jump_user_defined_operators).

# !!! tip
#     This tutorial is for advanced users. As an alternative, consider using
#     [Function tracing](@ref) instead of creating an operator, and if an
#     operator is necessary, consider using the default of `@operator(model, op_f, N, f)`
#     instead of passing explicit [Gradients and Hessians](@ref).

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import DifferentiationInterface
import Enzyme
import ForwardDiff
import Ipopt
import Test

# ## Primal function

# As a simple example, we consider the Rosenbrock function:

f(x::T...) where {T} = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2

# Here's the value at a random point:

x = rand(2)

#-

f(x...)

# ## Analytic derivative

# If expressions are simple enough, you can provide analytic functions for the
# gradient and Hessian.

# ### Gradient

# The Rosenbrock function has the gradient vector:

function analytic_∇f(g::AbstractVector, x...)
    g[1] = 400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2
    g[2] = 200 * (x[2] - x[1]^2)
    return
end

# Let's evaluate it at the same vector `x`:

analytic_g = zeros(2)
analytic_∇f(analytic_g, x...)
analytic_g

# ### Hessian

# The Hessian matrix is:

function analytic_∇²f(H::AbstractMatrix, x...)
    H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
    ## H[1, 2] = -400 * x[1] <-- not needed because Hessian is symmetric
    H[2, 1] = -400 * x[1]
    H[2, 2] = 200.0
    return
end

# Note that because the Hessian is symmetric, JuMP requires that we fill in only
# the lower triangle.

analytic_H = zeros(2, 2)
analytic_∇²f(analytic_H, x...)
analytic_H

# ### JuMP example

# Putting our analytic functions together, we get:

function analytic_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @operator(model, op_rosenbrock, 2, f, analytic_∇f, analytic_∇²f)
    @objective(model, Min, op_rosenbrock(x[1], x[2]))
    optimize!(model)
    Test.@test is_solved_and_feasible(model)
    return value.(x)
end

analytic_rosenbrock()

# ## ForwardDiff

# Instead of analytic functions, you can use [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
# to compute derivatives.

# !!! info
#     If you do not specify a gradient or Hessian, JuMP will use ForwardDiff.jl
#     to compute derivatives by default. We provide this section as a worked
#     example of what is going on under the hood.

# ### Pros and cons

# The main benefit of ForwardDiff is that it is simple, robust, and works with a
# broad range of Julia syntax.

# The main downside is that `f` must be a function that accepts arguments of
# `x::Real...`. See [Common mistakes when writing a user-defined operator](@ref)
# for more details.

# ### Gradient

# The gradient can be computed using `ForwardDiff.gradient!`. Note that
# ForwardDiff expects a single `Vector{T}` argument, but we receive `x` as a
# tuple, so we need `y -> f(y...)` and `collect(x)` to get things in the right
# format.

function fdiff_∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
    ForwardDiff.gradient!(g, y -> f(y...), collect(x))
    return
end

# Let's check that we find the analytic solution:

fdiff_g = zeros(2)
fdiff_∇f(fdiff_g, x...)
Test.@test ≈(analytic_g, fdiff_g)

# ### Hessian

# The Hessian is a bit more complicated, but code to implement it is:

function fdiff_∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
    h = ForwardDiff.hessian(y -> f(y...), collect(x))
    for i in 1:N, j in 1:i
        H[i, j] = h[i, j]
    end
    return
end

# Let's check that we find the analytic solution:

fdiff_H = zeros(2, 2)
fdiff_∇²f(fdiff_H, x...)
Test.@test ≈(analytic_H, fdiff_H)

# ### JuMP example

# The code for computing the gradient and Hessian using ForwardDiff can be
# re-used for many operators. Thus, it is helpful to encapsulate it into the
# function:

"""
    fdiff_derivatives(f::Function) -> Tuple{Function,Function}

Return a tuple of functions that evaluate the gradient and Hessian of `f` using
ForwardDiff.jl.
"""
function fdiff_derivatives(f::Function)
    function ∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
        ForwardDiff.gradient!(g, y -> f(y...), collect(x))
        return
    end
    function ∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
        h = ForwardDiff.hessian(y -> f(y...), collect(x))
        for i in 1:N, j in 1:i
            H[i, j] = h[i, j]
        end
        return
    end
    return ∇f, ∇²f
end

# Here's an example using `fdiff_derivatives`:

function fdiff_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @operator(model, op_rosenbrock, 2, f, fdiff_derivatives(f)...)
    @objective(model, Min, op_rosenbrock(x[1], x[2]))
    optimize!(model)
    Test.@test is_solved_and_feasible(model)
    return value.(x)
end

fdiff_rosenbrock()

# ## Enzyme

# Another library for automatic differentiation in Julia is [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl).

# ### Pros and cons

# The main benefit of Enzyme is that it can produce fast derivatives for
# functions with many input arguments.

# The main downsides are that it may throw unusual errors if your code uses an
# unsupported feature of Julia and that it may have large compile times.

# !!! warning
#     The JuMP developers cannot help you debug error messages related to
#     Enzyme. If the operator works, it works. If not, we suggest you try a
#     different automatic differentiation library. See [juliadiff.org](https://juliadiff.org/)
#     for details.

# ### Gradient

# The gradient can be computed using `Enzyme.autodiff` with the
# `Enzyme.Reverse` mode. We need to wrap `x` in `Enzyme.Active` to indicate that
# we want to compute the derivatives with respect to these arguments.

function enzyme_∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
    g .= Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active.(x)...)[1]
    return
end

# Let's check that we find the analytic solution:

enzyme_g = zeros(2)
enzyme_∇f(enzyme_g, x...)
Test.@test ≈(analytic_g, enzyme_g)

# ### Hessian

# We can compute the Hessian in Enzyme using forward-over-reverse automatic
# differentiation.

# The code to implement the Hessian in Enzyme is complicated, so we will not
# explain it in detail; see the [Enzyme documentation](https://enzymead.github.io/Enzyme.jl/stable/generated/autodiff/#Vector-forward-over-reverse).

function enzyme_∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
    ## direction(i) returns a tuple with a `1` in the `i`'th entry and `0`
    ## otherwise
    direction(i) = ntuple(j -> Enzyme.Active(T(i == j)), N)
    ## As the inner function, compute the gradient using Reverse mode
    ∇f(x...) = Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, x...)[1]
    ## For the outer autodiff, use Forward mode.
    hess = Enzyme.autodiff(
        Enzyme.Forward,
        ∇f,
        ## Compute multiple evaluations of Forward mode, each time using `x` but
        ## initializing with a different direction.
        Enzyme.BatchDuplicated.(Enzyme.Active.(x), ntuple(direction, N))...,
    )[1]
    ## Unpack Enzyme's `hess` data structure into the matrix `H` expected by
    ## JuMP.
    for j in 1:N, i in 1:j
        H[j, i] = hess[j][i]
    end
    return
end

# Let's check that we find the analytic solution:

enzyme_H = zeros(2, 2)
enzyme_∇²f(enzyme_H, x...)
Test.@test ≈(analytic_H, enzyme_H)

# ### JuMP example

# The code for computing the gradient and Hessian using Enzyme can be re-used
# for many operators. Thus, it is helpful to encapsulate it into the function:

"""
    enzyme_derivatives(f::Function) -> Tuple{Function,Function}

Return a tuple of functions that evaluate the gradient and Hessian of `f` using
Enzyme.jl.
"""
function enzyme_derivatives(f::Function)
    function ∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
        g .= Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active.(x)...)[1]
        return
    end
    function ∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
        direction(i) = ntuple(j -> Enzyme.Active(T(i == j)), N)
        ∇f(x...) = Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, x...)[1]
        hess = Enzyme.autodiff(
            Enzyme.Forward,
            ∇f,
            Enzyme.BatchDuplicated.(Enzyme.Active.(x), ntuple(direction, N))...,
        )[1]
        for j in 1:N, i in 1:j
            H[j, i] = hess[j][i]
        end
        return
    end
    return ∇f, ∇²f
end

# Here's an example using `enzyme_derivatives`:

function enzyme_rosenbrock()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @operator(model, op_rosenbrock, 2, f, enzyme_derivatives(f)...)
    @objective(model, Min, op_rosenbrock(x[1], x[2]))
    optimize!(model)
    Test.@test is_solved_and_feasible(model)
    return value.(x)
end

enzyme_rosenbrock()

# ## DifferentiationInterface

# Julia offers [many different autodiff packages](https://juliadiff.org/).
# [DifferentiationInterface.jl](https://github.com/gdalle/DifferentiationInterface.jl)
# is a package that provides an abstraction layer across a few underlying
# autodiff libraries.

# !!! warning
#     The JuMP developers cannot help you debug error messages related to
#     DifferentiationInterface. If the operator works, it works. If not, we
#     suggest you directly try using a different automatic differentiation
#     library rather than the DI wrapper. See [juliadiff.org](https://juliadiff.org/)
#     for details.

# All the necessary information about your choice of underlying autodiff
# package is encoded in a "backend object" like this one:

DifferentiationInterface.AutoForwardDiff()

# This type comes from another package called [ADTypes.jl](https://github.com/SciML/ADTypes.jl),
# but DifferentiationInterface re-exports it. Other options include
# `AutoZygote()` and `AutoFiniteDiff()`.

# ### Gradient

# Apart from providing the backend object, the syntax below remains very
# similar:

function di_∇f(
    g::AbstractVector{T},
    x::Vararg{T,N};
    backend = DifferentiationInterface.AutoForwardDiff(),
) where {T,N}
    DifferentiationInterface.gradient!(splat(f), g, backend, collect(x))
    return
end

# Let's check that we find the analytic solution:

di_g = zeros(2)
di_∇f(di_g, x...)
Test.@test ≈(analytic_g, di_g)

# ### Hessian

# The Hessian follows exactly the same logic, except we need only the lower
# triangle.

function di_∇²f(
    H::AbstractMatrix{T},
    x::Vararg{T,N};
    backend = DifferentiationInterface.AutoForwardDiff(),
) where {T,N}
    H_dense = DifferentiationInterface.hessian(splat(f), backend, collect(x))
    for i in 1:N, j in 1:i
        H[i, j] = H_dense[i, j]
    end
    return
end

# Let's check that we find the analytic solution:

di_H = zeros(2, 2)
di_∇²f(di_H, x...)
Test.@test ≈(analytic_H, di_H)

# ### JuMP example

# The code for computing the gradient and Hessian using DifferentiationInterface
# can be re-used for many operators. Thus, it is helpful to encapsulate it into
# the function:

"""
    di_derivatives(f::Function; backend) -> Tuple{Function,Function}

Return a tuple of functions that evaluate the gradient and Hessian of `f` using
DifferentiationInterface.jl with any given `backend`.
"""
function di_derivatives(f::Function; backend)
    function ∇f(g::AbstractVector{T}, x::Vararg{T,N}) where {T,N}
        DifferentiationInterface.gradient!(splat(f), g, backend, collect(x))
        return
    end
    function ∇²f(H::AbstractMatrix{T}, x::Vararg{T,N}) where {T,N}
        H_dense =
            DifferentiationInterface.hessian(splat(f), backend, collect(x))
        for i in 1:N, j in 1:i
            H[i, j] = H_dense[i, j]
        end
        return
    end
    return ∇f, ∇²f
end

# Here's an example using `di_derivatives`:

function di_rosenbrock(; backend)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @operator(model, op_rosenbrock, 2, f, di_derivatives(f; backend)...)
    @objective(model, Min, op_rosenbrock(x[1], x[2]))
    optimize!(model)
    Test.@test is_solved_and_feasible(model)
    return value.(x)
end

di_rosenbrock(; backend = DifferentiationInterface.AutoForwardDiff())
