# Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors    #src
# This Source Code Form is subject to the terms of the Mozilla Public License   #src
# v.2.0. If a copy of the MPL was not distributed with this file, You can       #src
# obtain one at https://mozilla.org/MPL/2.0/.                                   #src

# # Chordal decomposition

# The purpose of this tutorial is to show how to use [MathOptChordalDecomposition.jl](@ref)
# to improve the performance of models with PSD constraints.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import Downloads
import LinearAlgebra
import MathOptChordalDecomposition
import SCS
import SparseArrays

# ## Background

# Chordal decomposition is a technique for decomposing a large PSD constraint
# into a set of smaller PSD constraints and some linear equality constraints.

# If the original PSD constraint is sparse, the decomposed problem can be faster
# to solve than the original.

# For more information on chordal decomposition, watch Michael Garstka's talk at
# [JuMP-dev 2019](https://www.youtube.com/watch?v=H4Q0ZXDqB70).

# Some solvers, such as [Clarabel.jl](@ref) and [COSMO.jl](@ref) implement
# chordal decomposition internally. Others, such as [SCS.jl](@ref) do not
# implement chordal decomposition.

# The Julia package [MathOptChordalDecomposition.jl](@ref) is a MathOptInterface
# layer that implements chordal decomposition of sparse semidefinite constraints.
# It can be used to wrap any solver which supports PSD constraints.

# ## JuMP Model

# To demonstrate the benefits of chordal decomposition, we use the `mcp124-1`
# model from [SDPLIB](http://euler.nmt.edu/~brian/sdplib/sdplib.html).

file = "mcp124-1.dat-s"
dir = mktempdir()
filename = joinpath(dir, file)
Downloads.download(
    "https://raw.githubusercontent.com/vsdp/SDPLIB/refs/heads/master/data/$file",
    filename,
)
## This line is needed to work-around a bug in MathOptInterface v1.40 and
## earlier
write(filename, replace(read(filename, String), r"[,{}]" => " "))
model = read_from_file(filename)

# This model has 124 decision variables and one PSD constraint. This PSD
# constraint is sparse, which means that many elements of the matrix are zero.

# To view the matrix, use [`all_constraints`](@ref) to get a list of the
# constraints, then use [`constraint_object`](@ref) to get the function and set
# form of the constraint:

S = MOI.PositiveSemidefiniteConeTriangle
constraints = all_constraints(model, Vector{AffExpr}, S)
con = constraint_object(constraints[1]);
con.set
con.func

# The constraint function is given in vectorized form. Use [`reshape_vector`](@ref)
# to convert it into a matrix:

F = reshape_vector(con.func, SymmetricMatrixShape(con.set.side_dimension))

# The `F` matrix is dense, but many elements are zero. Use `SparseArrays.sparse`
# to turn it into a sparse matrix:

A = SparseArrays.sparse(F)

# The sparse matrix has 422 nonzeros, which is a density of 2.7%:

SparseArrays.nnz(A) / size(A, 1)^2

# ## Solution speed

# [SCS.jl](@ref) is a first-order solver that does not exploit the sparsity of
# PSD constraints. Let's solve it and see how long it took:

set_optimizer(model, SCS.Optimizer)
@time optimize!(model)

# In comparison, if we wrap `SCS.Optimizer` in `MathOptChordalDecomposition.Optimizer`,
# then the problem takes less than 1 second to solve:

set_optimizer(model, () -> MathOptChordalDecomposition.Optimizer(SCS.Optimizer))
@time optimize!(model)

# The difference in performance is because of the chordal decomposition. The
# decomposed problem introduced new variables (there are now 1,155 variables
# instead of 124) and constraints (there are now 115 PSD constraints instead of
# one), but each PSD constraint is smaller than the original.

decom = unsafe_backend(model)

# With a bit of effort, we can compute the number of PSD constraints of each
# size:

count_by_size = Dict{Int,Int}()
for ci in MOI.get(decom, MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S}())
    set = MOI.get(decom, MOI.ConstraintSet(), ci)
    n = set.side_dimension
    count_by_size[n] = get(count_by_size, n, 0) + 1
end
count_by_size
