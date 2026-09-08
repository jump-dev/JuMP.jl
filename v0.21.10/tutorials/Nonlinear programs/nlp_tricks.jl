# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Nonlinear tips and tricks

# This example collates some tips and tricks you can use when formulating
# nonlinear programs. It uses the following packages:

using JuMP
import Ipopt
import Test

# ## User-defined functions with vector outputs

# A common situation is to have a user-defined function like the following that
# returns multiple outputs (we define `function_calls` to keep track of how
# many times we call this method):

function_calls = 0
function foo(x, y)
    global function_calls += 1
    common_term = x^2 + y^2
    term_1 = sqrt(1 + common_term)
    term_2 = common_term
    return term_1, term_2
end

# For example, the first term might be used in the objective, and the second
# term might be used in a constraint, and often they share share work that is
# expensive to evaluate.

# This is a problem for JuMP, because it requires user-defined functions to
# return a single number. One option is to define two separate functions, the
# first returning the first argument, and the second returning the second
# argument.

foo_1(x, y) = foo(x, y)[1]
foo_2(x, y) = foo(x, y)[2]

# However, if the common term is expensive to compute, this approach is wasteful
# because it will evaluate the expensive term twice. Let's have a look at how
# many times we evaluate `x^2 + y^2` during a solve:

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[1:2] >= 0, start = 0.1)
register(model, :foo_1, 2, foo_1; autodiff = true)
register(model, :foo_2, 2, foo_2; autodiff = true)
@NLobjective(model, Max, foo_1(x[1], x[2]))
@NLconstraint(model, foo_2(x[1], x[2]) <= 2)
function_calls = 0
optimize!(model)
Test.@test objective_value(model) ≈ √3 atol = 1e-4
Test.@test value.(x) ≈ [1.0, 1.0] atol = 1e-4
println("Naive approach: function calls = $(function_calls)")
naive_approach = function_calls  #src

# An alternative approach is to use _memoization_, which uses a cache to store
# the result of function evaluations. We can write a memoization function as
# follows:

"""
    memoize(foo::Function, n_outputs::Int)

Take a function `foo` and return a vector of length `n_outputs`, where each
element is a function that returns the `i`'th output of `foo`.

To avoid duplication of work, cache the most-recent evaluations of `foo`.
Because `foo_i` is auto-differentiated with ForwardDiff, our cache needs to
work when `x` is a `Float64` and a `ForwardDiff.Dual`.
"""
function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx = nothing, nothing
    function foo_i(i, x::T...) where {T<:Real}
        if T == Float64
            if x != last_x
                last_x, last_f = x, foo(x...)
            end
            return last_f[i]::T
        else
            if x != last_dx
                last_dx, last_dfdx = x, foo(x...)
            end
            return last_dfdx[i]::T
        end
    end
    return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end

# Let's see how it works. First, construct the memoized versions of `foo`:

memoized_foo = memoize(foo, 2)

# Now try evaluating the first element of `memoized_foo`.

function_calls = 0
memoized_foo[1](1.0, 1.0)
Test.@test function_calls == 1  #src
println("function_calls = ", function_calls)

# As expected, this evaluated the function once. However, if we call the
# function again, we hit the cache instead of needing to re-compute `foo` and
# `function_calls` is still `1`!

memoized_foo[1](1.0, 1.0)
Test.@test function_calls == 1  #src
println("function_calls = ", function_calls)

# Now let's see how this works during a real solve:

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[1:2] >= 0, start = 0.1)
register(model, :foo_1, 2, memoized_foo[1]; autodiff = true)
register(model, :foo_2, 2, memoized_foo[2]; autodiff = true)
@NLobjective(model, Max, foo_1(x[1], x[2]))
@NLconstraint(model, foo_2(x[1], x[2]) <= 2)
function_calls = 0
optimize!(model)
Test.@test objective_value(model) ≈ √3 atol = 1e-4
Test.@test value.(x) ≈ [1.0, 1.0] atol = 1e-4
println("Memoized approach: function_calls = $(function_calls)")
Test.@test function_calls <= naive_approach / 2 + 1  #src

# Compared to the naive approach, the memoized approach requires half as many
# function evaluations!
