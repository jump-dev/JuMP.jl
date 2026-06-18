# Copyright (c) 2019 Matthew Help and contributors                               #src
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

# # N-Queens

# **This tutorial was originally contributed by Matthew Helm and Mathieu Tanneau.**

# This tutorial solves the N-Queens problem—placing N non-attacking queens on
# an N×N chessboard—as a binary integer program in JuMP. It is a compact
# example of how combinatorial feasibility problems can be modelled with row,
# column, and diagonal constraints.
#
# **Learning intentions:**
# * Model an N×N binary decision grid and enforce row and column exclusivity
#   as equality constraints
# * Enforce at-most-one-queen constraints on every diagonal by iterating over
#   the diagonals of both the original and row-reversed matrix
# * Recognise this as a pure feasibility problem with no objective, and extract
#   the integer solution by rounding the continuous values returned by the solver

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import HiGHS
import LinearAlgebra

# ## Formulation

# Here is an example of an N-Queens problem with four queens.

# ![Four Queens](../../assets/n_queens4.png)

# *Note that none of the queens above are able to attack any other as a result
# of their careful placement.*

# N-Queens
N = 8

model = Model(HiGHS.Optimizer)
set_silent(model)

# Next, let's create an N x N chessboard of binary values. 0 will represent an
# empty space on the board and 1 will represent a space occupied by one of our
# queens:

@variable(model, x[1:N, 1:N], Bin);

# Now we can add our constraints:

# There must be exactly one queen in a given row/column
for i in 1:N
    @constraint(model, sum(x[i, :]) == 1)
    @constraint(model, sum(x[:, i]) == 1)
end

# There can only be one queen on any given diagonal
for i in (-(N-1)):(N-1)
    @constraint(model, sum(LinearAlgebra.diag(x, i)) <= 1)
    @constraint(model, sum(LinearAlgebra.diag(reverse(x; dims = 1), i)) <= 1)
end

# We are ready to put our model to work and see if it is able to find
# a feasible solution:

optimize!(model)
assert_is_solved_and_feasible(model)

# We can now review the solution that our model found:

solution = round.(Int, value.(x))
