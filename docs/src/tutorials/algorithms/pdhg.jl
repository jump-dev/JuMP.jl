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

# # Writing a solver interface

# The purpose of this tutorial is to demonstrate how to implement a basic solver
# interface to MathOptInterface. As a motivating example, we implement the
# Primal Dual Hybrid Gradient (PDHG) method. PDHG is a first-order method that
# can solve convex optimization problems.
#
# Google has a [good introduction to the math behind PDLP](https://developers.google.com/optimization/lp/pdlp_math),
# which is a variant of PDHG specialized for linear programs.

# ## Required packages

# This tutorial requires the following packages:

using JuMP
import LinearAlgebra
import MathOptInterface as MOI
import Printf
import SparseArrays

# ## Primal Dual Hybrid Gradient

# The following function is a pedagogical implementation of PDHG that solves the
# linear program:
# ```math
# \begin{aligned}
# \text{min}        \ & c^\top x   \\
# \text{subject to} \ & Ax = b     \\
#                     & x \ge 0.
# \end{aligned}
# ```
#
# Note that this implementation is intentionally kept simple. It is not robust
# nor efficient, and it does not incorporate the theoretical improvements in the
# PDLP paper.

function solve_pdhg(
    A::SparseArrays.SparseMatrixCSC{Float64,Int},
    b::Vector{Float64},
    c::Vector{Float64};
    maximum_iterations::Int = 100_000,
    tol::Float64 = 1e-4,
    verbose::Bool = true,
    log_frequency::Int = 1_000,
)
    printf(x::Float64) = Printf.@sprintf("% 1.6e", x)
    printf(x::Int) = Printf.@sprintf("%6d", x)
    m, n = size(A)
    η = τ = 1 / LinearAlgebra.norm(A) - 1e-6
    x, y, k, status = zeros(n), zeros(m), 0, MOI.OTHER_ERROR
    if verbose
        println(
            "  iter      pobj          dobj         pfeas         dfeas       objfeas",
        )
    end
    while status == MOI.OTHER_ERROR
        x_next = max.(0.0, x - η * (A' * y + c))
        y += τ * (A * (2 * x_next - x) - b)
        x = x_next
        k += 1
        pfeas = LinearAlgebra.norm(A * x - b)
        dfeas = LinearAlgebra.norm(max.(0.0, -A' * y - c))
        objfeas = abs(c' * x - -b' * y)
        if pfeas <= tol && dfeas <= tol && objfeas <= tol
            status = MOI.OPTIMAL
        elseif k == maximum_iterations
            status = MOI.ITERATION_LIMIT
        end
        if verbose && (mod(k, log_frequency) == 0 || status != MOI.OTHER_ERROR)
            logs = printf.((k, c' * x, -b' * y, pfeas, dfeas, objfeas))
            println(join(logs, " "))
        end
    end
    return status, k, x, y
end

# Here's an example:

A = [0.0 -1.0 -1.0 0.0 0.0; 6.0 8.0 0.0 -1.0 0.0; 7.0 12.0 0.0 0.0 -1.0]
b = [-3.0, 100.0, 120.0]
c = [12.0, 20.0, 0.0, 0.0, 0.0]
status, k, x, y = solve_pdhg(SparseArrays.sparse(A), b, c);

# The termination status is:

status

# The solve took the following number of iterations:

k

# The primal solution is:

x

# The dual multipliers are:

y

# ## The MOI interface

# Converting a linear program from the modeler's form into the `A`, `b`, and `c`
# matrices of the standard form required by our implementation of PDHG is
# tedious and error-prone. This section walks through how to implement a basic
# interface to MathOptInterface, so that we can use our algorithm from JuMP.

# For a more comprehensive guide, see [Implementing a solver interface](@ref).

# ### The Optimizer type

# Create an optimizer by subtyping [`MOI.AbstractOptimizer`](@ref). By
# convention, the name of this type is `Optimizer`, and most optimizers are
# available as `PackageName.Optimizer`.

# The fields inside the optimizer are arbitrary. Store whatever is useful.

"""
    Optimizer()

Create a new optimizer for PDHG.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    ## A mapping from variable to column
    x_to_col::Dict{MOI.VariableIndex,Int}
    ## A mapping from constraint to rows
    ci_to_rows::Dict{
        MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},MOI.Zeros},
        Vector{Int},
    }
    ## Information from solve_pdhg
    status::MOI.TerminationStatusCode
    iterations::Int
    x::Vector{Float64}
    y::Vector{Float64}
    ## Other useful quantities
    solve_time::Float64
    obj_value::Float64

    function Optimizer()
        F = MOI.VectorAffineFunction{Float64}
        return new(
            Dict{MOI.VariableIndex,Int}(),
            Dict{MOI.ConstraintIndex{F,MOI.Zeros},Vector{Int}}(),
            MOI.OPTIMIZE_NOT_CALLED,
            0,
            Float64[],
            Float64[],
            0.0,
            0.0,
        )
    end
end

# Now that we have an `Optimizer`, we need to implement two methods:
# [`MOI.is_empty`](@ref) and [`MOI.empty!`](@ref). These are called whenever MOI
# needs to ensure that the optimizer is in a clean state.

function MOI.is_empty(model::Optimizer)
    ## You might want to check every field, not just a few
    return isempty(model.x_to_col) && model.status == MOI.OPTIMIZE_NOT_CALLED
end

function MOI.empty!(model::Optimizer)
    empty!(model.x_to_col)
    empty!(model.ci_to_rows)
    model.status = MOI.OPTIMIZE_NOT_CALLED
    model.iterations = 0
    model.solve_time = 0.0
    model.obj_value = 0.0
    empty!(model.x)
    empty!(model.y)
    return
end

# Next, we need to define what constraints the optimizer supports. Since our
# standard form was $Ax = b$, we support only $Ax + b \in \{0\}$, which is a
# [`MOI.VectorAffineFunction`](@ref) in [`MOI.Zeros`](@ref) constraint. Note
# that you might have expected $Ax - b \in \{0\}$. We'll address the difference
# in the sign of `b` in a few places later on.

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{Float64}},
    ::Type{MOI.Zeros},
)
    return true
end

# By default, MOI assumes that it can add free variables. This isn't true for
# our standard form, because we support only $x \ge 0$. Let's tell MOI that:

MOI.supports_add_constrained_variables(::Optimizer, ::Type{MOI.Reals}) = false

function MOI.supports_add_constrained_variables(
    ::Optimizer,
    ::Type{MOI.Nonnegatives},
)
    return true
end

# The objective function that we support is [`MOI.ScalarAffineFunction`](@ref):

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

# Finally, we'll implement [`MOI.SolverName`](@ref) so that MOI knows how to
# print the name of our optimizer:

MOI.get(::Optimizer, ::MOI.SolverName) = "PDHG"

# ### GenericModel

# The simplest way to solve a problem with your optimizer is to implement the
# method `MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)`, where `src` is an
# input model and `dest` is your empty optimizer.

# To implement this method you would need to query the variables and constraints
# in `src` and the convert these into the matrix data expected by `solve_pdhg`.
# Since matrix input is a common requirement of solvers, MOI includes utilities
# to simplify the process.

# The downside of the utilities is that they involve a highly parameterized type
# with a large number of possible configurations.The upside of the utilities is
# that, once setup, they requires few lines of code to extract the problem
# matrices.

# First, we need to define the set of sets that our standard form supports. For
# PDHG, we support only `Ax + b in {0}`:

MOI.Utilities.@product_of_sets(SetOfZeros, MOI.Zeros)

# Then, we define a [`MOI.Utilities.GenericModel`](@ref). This is the highly
# parameterized type that can be customized.

const CacheModel = MOI.Utilities.GenericModel{
    ## The coefficient type is Float64
    Float64,
    ## We use the default objective container
    MOI.Utilities.ObjectiveContainer{Float64},
    ## We use the default variable container
    MOI.Utilities.VariablesContainer{Float64},
    ## We use a Matrix of Constraints to represent `A * x + b in K`
    MOI.Utilities.MatrixOfConstraints{
        ## The number type is Float64
        Float64,
        ## The matrix type `A` is a sparse matrix
        MOI.Utilities.MutableSparseMatrixCSC{
            ## ... with Float64 coefficients
            Float64,
            ## ... Int64 row and column indices
            Int,
            ## ... and it uses one-based indexing
            MOI.Utilities.OneBasedIndexing,
        },
        ## The vector type `b` is a Julia `Vector`
        Vector{Float64},
        ## The set type `K` is the SetOfZeros type we defined above
        SetOfZeros{Float64},
    },
}

# As one example of possible alternate configuration, if you were interfacing
# with a solver written in C that expected zero-based indices, you might use
# instead:

MOI.Utilities.MutableSparseMatrixCSC{
    Cdouble,
    Cint,
    MOI.Utilities.ZeroBasedIndexing,
}

# !!! tip
#     The best place to look at how to configure `GenericModel` is to find an
#     existing solver with the same input standard form that you require.

# We need to make one modification to `CacheModel` to tell MOI that
# $x \in \mathbb{R}_+$ is equivalent to adding variables in
# [`MOI.GreaterThan`](@ref):

function MOI.add_constrained_variables(model::CacheModel, set::MOI.Nonnegatives)
    x = MOI.add_variables(model, MOI.dimension(set))
    MOI.add_constraint.(model, x, MOI.GreaterThan(0.0))
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Nonnegatives}(x[1].value)
    return x, ci
end

# ### The optimize method

# Now we define the most important method for our optimizer.

function MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)
    ## In addition to the values returned by `solve_pdhg`, it may be useful to
    ## record other attributes, such as the solve time.
    start_time = time()
    ## Construct a cache to store our problem data:
    cache = CacheModel()
    ## MOI includes a utility to copy an arbitrary `src` model into `cache`. The
    ## return, `index_map`, is a mapping from indices in `src` to indices in
    ## `dest`.
    index_map = MOI.copy_to(cache, src)
    ## Now we can access the `A` matrix:
    A = convert(
        SparseArrays.SparseMatrixCSC{Float64,Int},
        cache.constraints.coefficients,
    )
    ## and the b vector (note that MOI models Ax = b as Ax + b in {0}, so b
    ## differs by -):
    b = -cache.constraints.constants
    ## The `c` vector is more involved, because we need to account for the
    ## objective sense:
    sense = ifelse(cache.objective.sense == MOI.MAX_SENSE, -1, 1)
    F = MOI.ScalarAffineFunction{Float64}
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())
    c = zeros(size(A, 2))
    for term in obj.terms
        c[term.variable.value] += sense * term.coefficient
    end
    ## Now we can solve the problem with PDHG and record the solution:
    dest.status, dest.iterations, dest.x, dest.y = solve_pdhg(A, b, c)
    ## To help assign the values of the x and y vectors to the appropriate
    ## variables and constrats, we need a map of the constraint indices to their
    ## row in the `dest` matrix and a map of the variable indices to their
    ## column in the `dest` matrix:
    F, S = MOI.VectorAffineFunction{Float64}, MOI.Zeros
    for src_ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        dest.ci_to_rows[index_map[src_ci]] =
            MOI.Utilities.rows(cache.constraints.sets, index_map[src_ci])
    end
    for (i, src_x) in enumerate(MOI.get(src, MOI.ListOfVariableIndices()))
        dest.x_to_col[index_map[src_x]] = i
    end
    ## We can now record two derived quantities: the primal objective value and
    ## the solve time.
    dest.obj_value = obj.constant + sense * c' * dest.x
    dest.solve_time = time() - start_time
    ## We need to return the index map, and `false`, to indicate to MOI that we
    ## do not support incremental modification of the model.
    return index_map, false
end

# ## Solutions

# Now that we know how to solve a model, let's implement the required solution
# attributes.

# First, we need to tell MOI how many solutions we found via
# [`MOI.ResultCount`](@ref):

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return model.status == MOI.OPTIMAL ? 1 : 0
end

# and implement [`MOI.RawStatusString`](@ref) to provide a user-readable string
# that describes what happened:

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.status == MOI.OPTIMAL
        return "found a primal-dual optimal solution (subject to tolerances)"
    end
    return "failed to solve"
end

# Then, we need to implement the three types of problem status:
# [`MOI.TerminationStatus`](@ref), [`MOI.PrimalStatus`](@ref) and
# [`MOI.DualStatus`](@ref):

MOI.get(model::Optimizer, ::MOI.TerminationStatus) = model.status

function MOI.get(model::Optimizer, attr::Union{MOI.PrimalStatus,MOI.DualStatus})
    if attr.result_index == 1 && model.status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

# Now we can implement [`MOI.ObjectiveValue`](@ref), [`MOI.VariablePrimal`](@ref),
# and [`MOI.ConstraintDual`](@ref):

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return model.obj_value
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    return model.x[model.x_to_col[x]]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},MOI.Zeros},
)
    MOI.check_result_index_bounds(model, attr)
    ## MOI models Ax = b as Ax + b in {0}, so the dual differs by -
    return -model.y[model.ci_to_rows[ci]]
end

# Some other useful result quantities are [`MOI.SolveTimeSec`](@ref) and
# [`MOI.BarrierIterations`](@ref):

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time
MOI.get(model::Optimizer, ::MOI.BarrierIterations) = model.iterations

# ## A JuMP example

# Now we can solve an arbitrary linear program with JuMP. Here's the same
# standard form as before:

model = Model(Optimizer)
@variable(model, x[1:5] >= 0)
@objective(model, Min, c' * x)
@constraint(model, c3, A * x == b)
optimize!(model)

#-

solution_summary(model; verbose = true)

# But we could also have written:

model = Model(Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
optimize!(model)

#-

solution_summary(model; verbose = true)

# Other variations are also possible:

model = Model(Optimizer)
@variable(model, x[1:5] >= 0)
@objective(model, Max, -c' * x)
@constraint(model, c4, A * x .== b)
optimize!(model)

#-

solution_summary(model; verbose = true)

# Behind the scenes, JuMP and MathOptInterface reformulate the problem from the
# modeller's form into the standard form defined by our `Optimizer`.
