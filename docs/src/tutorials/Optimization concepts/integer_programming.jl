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

# # Integer Programming

# **Originally Contributed by**: Arpit Bhatia

# While we already know how to set a variable as integer or binary in the
# [`@variable`](@ref) macro, this tutorial covers other JuMP features for
# integer programming along with some modelling techniques.

using JuMP
import Random

Random.seed!(1234)

# ## Modelling Logical Conditions

# Generally, in a mathematical programming problem, all constraints must hold.
# However, we might want to have conditions where we have some logical
# conditions between constraints. In such cases, we can use binary variables for
# modelling logical conditions between constraints.

# ### Disjunctive Constraints (OR)

# Suppose that we are given two constraints $a'x \geq b$ and $c' x \geq d$,
# in which all components of $a$ and $c$ are non-negative.

# We would like to model a requirement that at least one of the two constraints
# is satisfied. For this, we defined a binary variable $y$ and impose the
# constraints:

# ```math
# \begin{aligned}
# a' x \geq y b \\
# c' x \geq (1 - y) d \\
# y \in \{0,1\}
# \end{aligned}
# ```

a = rand(1:100, 5, 5)
c = rand(1:100, 5, 5)
b = rand(1:100, 5)
d = rand(1:100, 5)

model = Model()
@variable(model, x[1:5])
@variable(model, y, Bin)
@constraint(model, a * x .>= y .* b)
@constraint(model, c * x .>= (1 - y) .* d);

# ### Conditional Constraints ($\implies$)

# Suppose we want to model that a certain linear inequality must be satisfied
# when some other event occurs, i.e., for a binary variable $z$, we want to
# model the implication:

# ```math
# \begin{aligned}
# z = 1 \implies a^Tx\leq b
# \end{aligned}
# ```

# If we know in advance an upper bound $a^Tx\leq b$. Then we can write the above
# as a linear inequality

# ```math
# \begin{aligned}
# a^Tx\leq b + M(1-z)
# \end{aligned}
# ```

a = rand(1:100, 5, 5)
b = rand(1:100, 5)
m = rand(10000:11000, 5)

model = Model()
@variable(model, x[1:5])
@variable(model, z, Bin)
@constraint(model, a * x .<=  b .+ (m .* (1 - z)));

# If z was a regular Julia variable, we would not have had to use the vectorized
# dot operator

# ### Boolean Operators on Binary Variables

# The following table is useful when we want to model boolean operators in the
# form of linear inequalities that can be given to a solver.

# | Boolean Expression | Constraint                           |
# |:----------         |                           ----------:|
# | $z=x \lor y$       | $x \leq z,  y \leq z,  z \leq x+y$   |
# | $z=x \land y$      | $x \geq z,  y \geq z,  z+1 \geq x+y$ |
# | $z= \neg x$        | $z = 1 − x$                          |
# | $x \implies y$     | $x \leq y$                           |
# | $x \iff y$         | $x = y$                              |

# ## Modelling Integer Variables

# ### Integer Variables using Constraints
# We can add binary and integer restrictions to the domain of each variable
# using the [`@constraint`](@ref) macro as well.

model = Model()
@variable(model, x)
@variable(model, y)
@constraint(model, x in MOI.ZeroOne())
@constraint(model, y in MOI.Integer())

# ### Semi-Continuous Variables

# A semi-continuous variable is a continuous variable between bounds $[l,u]$
# that also can assume the value zero. ie.
# $$x \in \{0\} \cup \{l,u\}.$$

l = 7.45
u = 22.22
@variable(model, a)
@constraint(model, a in MOI.Semicontinuous(l, u))

# ### Semi-Integer Variables

# A semi-integer variable is a variable which asummes integer values between
# bounds $[l,u]$ and can also assume the value zero:
# $$x \in \{0\} \cup (\{l,u\} \cap \mathbb{Z}).$$

l = 5
u = 34
@variable(model, b)
@constraint(model, b in MOI.Semiinteger(l, u))

# Note that the bounds specified in `MOI.Semiinteger` must be integral otherwise
# it would throw an error.

# ## Special Ordered Sets

# ### Special Ordered Sets of Type I (SOS1)

# A Special Ordered Set of type 1 is a set of variables, at most one of which
# can take a non-zero value, all others being at 0.
#
# They most frequently apply where a set of variables are actually binary
# variables. In other words, we have to choose at most one from a set of
# possibilities.

@variable(model, u[1:3])
@constraint(model, u in MOI.SOS1([1.0, 2.0, 3.0]))

# Note that we have to pass `MOI.SOS1` a weight vector which is essentially an
# ordering on the variables.

# If the decision variables are related and have a physical ordering, then the
# weight vector, although not used directly in the constraint, can help the
# solver make a better decision in the solution process.

# ### Special Ordered Sets of Type II (SOS2)

# A Special Ordered Set of type 2 is a set of non-negative variables, of which
# at most two can be non-zero, and if two are non-zero these must be consecutive
# in their ordering.

@variable(model, v[1:3])
@constraint(model, v in MOI.SOS2([3.0, 1.0, 2.0]))

# The ordering provided by the weight vector is more important in this case as
# the variables need to be consecutive according to the ordering.
# For example, in the above constraint, the possible pairs are:
# * Consecutive
#   * (`x[1]` and `x[3]`) as they correspond to 3 and 2 resp. and thus can be non-zero
#   * (`x[2]` and `x[3]`) as they correspond to 1 and 2 resp. and thus can be non-zero
# * Non-consecutive
#   * (`x[1]` and `x[2]`) as they correspond to 3 and 1 resp. and thus cannot be non-zero
