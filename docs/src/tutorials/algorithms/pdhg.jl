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

# The purpose of this tutorial is to demonstrate

using JuMP
import LinearAlgebra
import MathOptInterface as MOI
import Printf
import SparseArrays

# ## PDHG

"""
    solve_pdhg(
        A::SparseArrays.SparseMatrixCSC{Float64,Int},
        b::Vector{Float64},
        c::Vector{Float64};
        maximum_iterations::Int = 10_000,
        tol::Float64 = 1e-4,
        verbose::Bool = true,
        log_frequency::Int = 1_000,
    ) -> status, iterations, x, y

A pedagogical implementation of PDHG that solves the linear program:

```math
\\begin{aligned}
\\text{min}        \\ & c^\\top x \\
\\text{subject to} \\ & Ax = b \\
                      & x \\ge 0.
\\end{aligned}
```

Note that this implementation is intentionally kept simple. It is not robust nor
efficient, and it does not incorporate the theoretical improvements in the PDLP
paper.

## Keyword arguments

 * `maximum_iterations::Int = 10_000`: the maximum number of iterations before
   termination.
 * `tol::Float64 = 1e-4`: the combined primal, dual, and strong duality
   tolerance.
 * `verbose::Bool = true`: print iteration logs
 * `log_frequency::Int = 1_000`: print iteration logs every N iterations
"""
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
        objfeas = abs(c' * x + b' * y)
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
status, k, x, y = solve_pdhg(SparseArrays.sparse(A), b, c)

# ## The MOI interface

# Converting example linear program from the modeler's form into the standard
# form required by PDHG is tedious and error-prone. This section walks through
# how to implement a basic interface to MathOptInterface, so that we can use our
# algorithm from JuMP.

"""
    Optimizer()

Create a new optimizer for PDHG.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    ## Information from solve_pdhg
    status::MOI.TerminationStatusCode
    iterations::Int
    x::Vector{Float64}
    y::Vector{Float64}
    ## Other useful quantities
    solve_time::Float64
    obj_value::Float64

    function Optimizer()
        return new(MOI.OPTIMIZE_NOT_CALLED, 0, Float64[], Float64[], 0.0, 0.0)
    end
end

# First, we need to implement two methods: [`MOI.is_empty`](@ref) and
# [`MOI.empty!`](@ref). These are called whenever MOI needs to ensure that the
# optimizer is in a clean state.

function MOI.is_empty(model::Optimizer)
    ## You might want to check every field, not just the status
    return model.status == MOI.OPTIMIZE_NOT_CALLED
end

function MOI.empty!(model::Optimizer)
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

# Finally, the objective function that we support is
# [`MOI.ScalarAffineFunction`](@ref):

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

# Finally, we'll implement [`MOI.SolverName`](@ref) so that MOI knows how to
# print the name of our optimizer:

MOI.get(::Optimizer, ::MOI.SolverName) = "PDHG"

MOI.Utilities.@product_of_sets(LinearZero, MOI.Zeros)

const Cache = MOI.Utilities.GenericModel{
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
        ## The set type `K` is the LinearZero set we defined above
        LinearZero{Float64},
    },
}

# If you were interfacing with a solver written in C that expected zero-based
# indices, you might use instead:
MOI.Utilities.MutableSparseMatrixCSC{
    Cdouble,
    Cint,
    MOI.Utilities.ZeroBasedIndexing,
}

function MOI.add_constrained_variables(model::Cache, set::MOI.Nonnegatives)
    x = MOI.add_variables(model, MOI.dimension(set))
    MOI.add_constraint.(model, x, MOI.GreaterThan(0.0))
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.Nonnegatives}(x[1].value)
    return x, ci
end

function MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)
    start_time = time()
    cache = Cache()
    index_map = MOI.copy_to(cache, src)
    A = convert(
        SparseArrays.SparseMatrixCSC{Float64,Int},
        cache.constraints.coefficients,
    )
    ## MOI models Ax = b as Ax + b in {0}, so b differs by -
    b = -cache.constraints.constants
    c = zeros(size(A, 2))
    sense = ifelse(cache.objective.sense == MOI.MAX_SENSE, -1, 1)
    F = MOI.ScalarAffineFunction{Float64}
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())
    for term in obj.terms
        c[term.variable.value] += sense * term.coefficient
    end
    dest.status, dest.iterations, dest.x, dest.y = solve_pdhg(A, b, c)
    dest.obj_value = obj.constant + c' * dest.x
    dest.solve_time = time() - start_time
    return index_map, false
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    return model.status == MOI.OPTIMAL ? 1 : 0
end

MOI.get(model::Optimizer, ::MOI.RawStatusString) = string(model.status)

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
    return model.x[x.value]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex,
)
    MOI.check_result_index_bounds(model, attr)
    ## MOI models Ax = b as Ax + b in {0}, so the dual differs by -
    return -model.y[ci.value]
end

# Some other useful result quantities are [`MOI.SolveTimeSec`](@ref) and
# [`MOI.BarrierIterations`](@ref):

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time
MOI.get(model::Optimizer, ::MOI.BarrierIterations) = model.iterations

# ## A JuMP example

# Now we can solve an arbitrary linear program with JuMP:

model = Model(Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
optimize!(model)
solution_summary(model; verbose = true)
