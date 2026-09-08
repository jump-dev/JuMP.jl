# Copyright (c) 2019 Arpit Bhatia and contributors                               #src
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

# # Tips and tricks

# **Originally Contributed by**: Arpit Bhatia

# !!! tip
#     A good source of tips is the [Mosek Modeling Cookbook](https://docs.mosek.com/modeling-cookbook/mio.html).

# This tutorial collates some tips and tricks you can use when formulating
# mixed-integer programs. It uses the following packages:

using JuMP

# ## Boolean Operators on Binary Variables

# Binary variables can be used to construct logical operators. Here are some
# example.

# ### OR

# $$x_3 = x_1 \lor x_2$$

model = Model()
@variable(model, x[1:3], Bin)
@constraints(model, begin
    x[1] <= x[3]
    x[2] <= x[3]
end)

# ### And

# $$x_3 = x_1 \land x_2$$

model = Model()
@variable(model, x[1:3], Bin)
@constraints(model, begin
    x[3] <= x[1]
    x[3] <= x[2]
    x[3] >= x[1] + x[2] - 1
end)

# ### Not

# $$x_1 \neg x_2$$

model = Model()
@variable(model, x[1:2], Bin)
@constraint(model, x[1] == 1 - x[2])

# ### Implies

# $$x_1 \implies x_2$$

model = Model()
@variable(model, x[1:2], Bin)
@constraint(model, x[1] <= x[2])

# ## Big-M Disjunctive Constraints (OR)

# **Problem** Suppose that we have two constraints $a^\top x \leq b$ and
# $c^\top x \leq d$, and we want at least one to hold.

# **Trick** Introduce a "big-M" multiplied by a binary variable to relax one of
# the constraints.

# **Example** Either $x_1 \leq 1$ and/or $x_2 \leq 2$.

model = Model()
@variable(model, x[1:2])
@variable(model, y, Bin)
M = 100
@constraint(model, x[1] <= 1 + M * y)
@constraint(model, x[2] <= 2 + M * (1 - y))

# !!! warning
#     If `M` is too small, the solution may be suboptimal. If `M` is too big,
#     the solver may encounter numerical issues. Try to use domain knowledge to
#     choose an `M` that is just right. Gurobi has a [good documentation section](https://www.gurobi.com/documentation/9.1/refman/dealing_with_big_m_constra.html)
#     on this topic.

# ## Indicator Constraints ($\implies$)

# ### Problem

# Suppose we want to model that a certain linear inequality must be satisfied
# when some other event occurs, i.e., for a binary variable $z$, we want to
# model the implication:

# $$z = 1 \implies a^Tx \leq b$$

# ### Trick 1

# Some solvers have native support for indicator constraints.

# **Example** $x_1 + x_2 \leq 1$ if $z = 1$.

model = Model()
@variable(model, x[1:2])
@variable(model, z, Bin)
@constraint(model, z => {sum(x) <= 1})

# **Example** $x_1 + x_2 \leq 1$ if $z = 0$.

model = Model()
@variable(model, x[1:2])
@variable(model, z, Bin)
@constraint(model, !z => {sum(x) <= 1})

# ### Trick 2

# If the solver doesn't support indicator constraints, you an use the big-M
# trick.

# **Example** $x_1 + x_2 \leq 1$ if $z = 1$.

model = Model()
@variable(model, x[1:2])
@variable(model, z, Bin)
M = 100
@constraint(model, sum(x) <= 1 + M * (1 - z))

# **Example** $x_1 + x_2 \leq 1$ if $z = 0$.

model = Model()
@variable(model, x[1:2])
@variable(model, z, Bin)
M = 100
@constraint(model, sum(x) <= 1 + M * z)

# ## Semi-Continuous Variables

# A semi-continuous variable is a continuous variable between bounds $[l,u]$
# that also can assume the value zero. ie.
# $$x \in \{0\} \cup [l,u].$$

# **Example** $$x \in \{0\}\cup [1, 2]$$

model = Model()
@variable(model, x in MOI.Semicontinuous(1.0, 2.0))

# ## Semi-Integer Variables

# A semi-integer variable is a variable which assumes integer values between
# bounds $[l,u]$ and can also assume the value zero:
# $$x \in \{0\} \cup [l, u] \cap \mathbb{Z}.$$

model = Model()
@variable(model, x in MOI.Semiinteger(5.0, 10.0))

# ## Special Ordered Sets of Type I (SOS1)

# A Special Ordered Set of Type I is a set of variables, at most one of which
# can take a non-zero value, all others being at 0.
#
# They most frequently apply where a set of variables are actually binary
# variables. In other words, we have to choose at most one from a set of
# possibilities.

model = Model()
@variable(model, x[1:3], Bin)
@constraint(model, x in SOS1())

# You can optionally pass `SOS1` a weight vector like

@constraint(model, x in SOS1([0.2, 0.5, 0.3]))

# If the decision variables are related and have a physical ordering, then the
# weight vector, although not used directly in the constraint, can help the
# solver make a better decision in the solution process.

# ## [Special Ordered Sets of Type II (SOS2)](@id tip_sos2)

# A Special Ordered Set of type 2 is a set of non-negative variables, of which
# at most two can be non-zero, and if two are non-zero these must be consecutive
# in their ordering.

model = Model()
@variable(model, x[1:3])
@constraint(model, x in MOI.SOS2([3.0, 1.0, 2.0]))

# The ordering provided by the weight vector is more important in this case as
# the variables need to be consecutive according to the ordering.
# For example, in the above constraint, the possible pairs are:
# * Consecutive
#   * (`x[1]` and `x[3]`) as they correspond to 3 and 2 resp. and thus can be non-zero
#   * (`x[2]` and `x[3]`) as they correspond to 1 and 2 resp. and thus can be non-zero
# * Non-consecutive
#   * (`x[1]` and `x[2]`) as they correspond to 3 and 1 resp. and thus cannot be non-zero

# ## Piecewise linear approximations

# [SOSII constraints](@ref tip_sos2) are most often used to form piecewise
# linear approximations of a function.

# Given a set of points for `x`:
x̂ = -1:0.5:2
# and a set of corresponding points for `y`:
ŷ = x̂ .^ 2
# the piecewise linear approximation is constructed by representing `x` and `y`
# as convex combinations of `x̂` and `ŷ`.

N = length(x̂)
model = Model()
@variable(model, -1 <= x <= 2)
@variable(model, y)
@variable(model, 0 <= λ[1:N] <= 1)
@objective(model, Max, y)
@constraints(model, begin
    x == sum(x̂[i] * λ[i] for i in 1:N)
    y == sum(ŷ[i] * λ[i] for i in 1:N)
    sum(λ) == 1
    λ in SOS2()
end)
