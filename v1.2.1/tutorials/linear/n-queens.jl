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

# **Originally Contributed by**: Matthew Helm ([with help from @mtanneau on Julia Discourse](https://discourse.julialang.org/t/which-jump-jl-solver-for-this-problem/43350/17?u=mthelm85))

# The N-Queens problem involves placing N queens on an N x N chessboard such
# that none of the queens attacks another. In chess, a queen can move
# vertically, horizontally, and diagonally so there cannot be more than one
# queen on any given row, column, or diagonal.

# ![Four Queens](../../assets/n_queens4.png)

# *Note that none of the queens above are able to attack any other as a result
# of their careful placement.*

using JuMP
import HiGHS
import LinearAlgebra

# N-Queens
N = 8

model = Model(HiGHS.Optimizer)

# Next, let's create an N x N chessboard of binary values. 0 will represent an
# empty space on the board and 1 will represent a space occupied by one of our
# queens:

@variable(model, x[1:N, 1:N], Bin)

# Now we can add our constraints:

# There must be exactly one queen in a given row/column
for i in 1:N
    @constraint(model, sum(x[i, :]) == 1)
    @constraint(model, sum(x[:, i]) == 1)
end

# There can only be one queen on any given diagonal
for i in -(N - 1):(N-1)
    @constraint(model, sum(LinearAlgebra.diag(x, i)) <= 1)
    @constraint(model, sum(LinearAlgebra.diag(reverse(x; dims = 1), i)) <= 1)
end

# That's it! We are ready to put our model to work and see if it is able to find
# a feasible solution:

optimize!(model)

# We can now review the solution that our model found:

solution = round.(Int, value.(x))
