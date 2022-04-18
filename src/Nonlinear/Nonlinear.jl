#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module Nonlinear

import Base.Meta: isexpr
import ForwardDiff
import MathOptInterface
import OrderedCollections: OrderedDict
import SparseArrays

using SpecialFunctions

const MOI = MathOptInterface

# Override basic math functions to return NaN instead of throwing errors.
# This is what NLP solvers expect, and sometimes the results aren't needed
# anyway, because the code may compute derivatives wrt constants.
import NaNMath:
    sin,
    cos,
    tan,
    asin,
    acos,
    acosh,
    atanh,
    log,
    log2,
    log10,
    lgamma,
    log1p,
    pow,
    sqrt

include("univariate_expressions.jl")
include("operators.jl")
include("types.jl")
include("parse.jl")

"""
    set_objective(data::NonlinearData, obj)::Nothing

Parse `obj` into a `NonlinearExpression` and set as the objective function of
`data`.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
set_objective(data, :(\$x^2 + 1))
set_objective(data, x)
```
"""
function set_objective(data::NonlinearData, obj)
    data.objective = parse_expression(data, obj)
    return
end

"""
    add_expression(data::NonlinearData, expr)::ExpressionIndex

Parse `expr` into a `NonlinearExpression` and add to `data`. Returns an
`ExpressionIndex` that can be interpolated into other input expressions.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
ex = add_expression(data, :(\$x^2 + 1))
set_objective(data, :(sqrt(\$ex)))
```
"""
function add_expression(data::NonlinearData, expr)
    push!(data.expressions, parse_expression(data, expr))
    return ExpressionIndex(length(data.expressions))
end

function Base.getindex(data::NonlinearData, index::ExpressionIndex)
    return data.expressions[index.value]
end

"""
    add_constraint(data::NonlinearData, input::Expr)::ConstraintIndex

Parse `input` into a `NonlinearConstraint` and add to `data`. Returns an
`ConstraintIndex` that can be used to delete the constraint or query solution
information.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
c = add_constraint(data, :(\$x^2 <= 1))
```
"""
function add_constraint(data::NonlinearData, input::Expr)
    expr, set = _expr_to_constraint(input)
    f = parse_expression(data, expr)
    data.last_constraint_index += 1
    index = ConstraintIndex(data.last_constraint_index)
    data.constraints[index] = NonlinearConstraint(f, set)
    return index
end

"""
    delete(data::NonlinearData, c::ConstraintIndex)::Nothing

Delete the constraint index `c` from `data`.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
c = add_constraint(data, :(\$x^2 <= 1))
delete(data, c)
```
"""
function delete(data::NonlinearData, c::ConstraintIndex)
    delete!(data.constraints, c)
    return
end

"""
    row(data::NonlinearData, c::ConstraintIndex)::Int

Return the row of the constraint index `c` in `data`.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
c = add_constraint(data, :(\$x^2 <= 1))
row(data, c)  # Returns 1
```
"""
function row(data::NonlinearData, c::ConstraintIndex)
    # TODO(odow): replace with a cache that maps indices to their 1-indexed
    # row in the constraint matrix. But failing that, since we know that
    # constraints are added in increasing order and that they can be deleted, we
    # know that index `i` must appear as constraint `1` to `i`. So we start at
    # `i` and backtrack (to account for deleted constraints) until we find it.
    # In the typical case with no deletion, there should be no overhead.
    for i in min(c.value, length(data.ordered_constraints)):-1:1
        if data.ordered_constraints[i] == c
            return i
        end
    end
    return error("Invalid constraint index $(c)")
end

function Base.getindex(data::NonlinearData, index::ConstraintIndex)
    return data.constraints[index]
end

function MOI.is_valid(data::NonlinearData, index::ConstraintIndex)
    return haskey(data.constraints, index)
end

"""
    add_parameter(data::NonlinearData, value::Float64)::ParameterIndex

Add a new parameter to `data` with the default value `value`. Returns a
`ParameterIndex` that can be interpolated into other input expressions and used
to modify the value of the parameter.

## Examples

```julia
data = NonlinearData()
x = MOI.VariableIndex(1)
p = add_parameter(data, 1.2)
c = add_constraint(data, :(\$x^2 <= \$p))
```
"""
function add_parameter(data::NonlinearData, value::Float64)
    push!(data.parameters, value)
    return ParameterIndex(length(data.parameters))
end

function Base.getindex(data::NonlinearData, p::ParameterIndex)
    return data.parameters[p.value]
end

function Base.setindex!(data::NonlinearData, value::Real, p::ParameterIndex)
    return data.parameters[p.value] = convert(Float64, value)::Float64
end

"""
    register_operator(
        data::NonlinearData,
        op::Symbol,
        nargs::Int,
        f::Function,
        [∇f::Function],
        [∇²f::Function],
    )

Register the user-defined operator `op` with `nargs` input arguments in `data`.

## Univariate functions

 * `f(x::T)::T` must be a function that takes a single input argument `x` and
   returns the function evaluated at `x`. If `∇f` and `∇²f` are not provided,
   `f` must support any `Real` input type `T`.
 * `∇f(x::T)::T` is a function that takes a single input argument `x` and
   returns the first derivative of `f` with respect to `x`. If `∇²f` is not
   provided, `∇f` must support any `Real` input type `T`.
 * `∇²f(x::T)::T` is a function that takes a single input argument `x` and
   returns the second derivative of `f` with respect to `x`.

## Multivariate functions

* `f(x::T...)::T` must be a function that takes a `nargs` input arguments `x`
  and returns the function evaluated at `x`. If `∇f` and `∇²f` are not provided,
  `f` must support any `Real` input type `T`.
* `∇f(g::AbstractVector{T}, x::T...)::T` is a function that takes a cache vector
  `g` of length `length(x)`, and fills each element `g[i]` with the partial
  derivative of `f` with respect to `x[i]`.

Hessian are not supported for multivariate functions.
"""
function register_operator(
    data::NonlinearData,
    op::Symbol,
    nargs::Int,
    f::Function...,
)
    return register_operator(data.operators, op, nargs, f...)
end

function Base.copy(::NonlinearData)
    return error("Copying nonlinear problems not yet implemented")
end

function MOI.features_available(data::NonlinearData)
    features = Symbol[]
    if data.inner !== nothing
        append!(features, MOI.features_available(data.inner))
    end
    if !(:ExprGraph in features)
        push!(features, :ExprGraph)
    end
    return features
end

function Base.show(io::IO, data::NonlinearData)
    Base.print(io, "NonlinearData with available features:")
    for feature in MOI.features_available(data)
        print(io, "\n  * :", feature)
    end
    return
end

function MOI.initialize(data::NonlinearData, features::Vector{Symbol})
    empty!(data.ordered_constraints)
    empty!(data.julia_expressions)
    data.eval_objective_timer = 0.0
    data.eval_objective_gradient_timer = 0.0
    data.eval_constraint_timer = 0.0
    data.eval_constraint_jacobian_timer = 0.0
    data.eval_hessian_lagrangian_timer = 0.0
    append!(data.ordered_constraints, keys(data.constraints))
    if :ExprGraph in features
        for i in 1:length(data.expressions)
            push!(
                data.julia_expressions,
                _to_expr(data, data.expressions[i]; use_x_ref = true),
            )
        end
        filter!(f -> f != :ExprGraph, features)
    end
    if data.inner !== nothing
        MOI.initialize(data.inner, features)
    end
    return
end

function MOI.objective_expr(data::NonlinearData)
    @assert data.objective !== nothing
    return _to_expr(data, data.objective; use_x_ref = true)
end

function MOI.constraint_expr(data::NonlinearData, i::Int)
    constraint = data[data.ordered_constraints[i]]
    f = _to_expr(data, constraint.expression; use_x_ref = true)
    if constraint.set isa MOI.LessThan
        return :($f <= $(constraint.set.upper))
    elseif constraint.set isa MOI.GreaterThan
        return :($f >= $(constraint.set.lower))
    elseif constraint.set isa MOI.EqualTo
        return :($f == $(constraint.set.value))
    else
        @assert constraint.set isa MOI.Interval
        return :($(constraint.set.lower) <= $f <= $(constraint.set.upper))
    end
end

function MOI.eval_objective(data::NonlinearData, x)
    start = time()
    obj = MOI.eval_objective(data.inner, x)
    data.eval_objective_timer += time() - start
    return obj
end

function MOI.eval_objective_gradient(data::NonlinearData, g, x)
    start = time()
    MOI.eval_objective_gradient(data.inner, g, x)
    data.eval_objective_gradient_timer += time() - start
    return
end

function MOI.eval_constraint(data::NonlinearData, g, x)
    start = time()
    MOI.eval_constraint(data.inner, g, x)
    data.eval_constraint_timer += time() - start
    return
end

function MOI.jacobian_structure(data::NonlinearData)
    return MOI.jacobian_structure(data.inner)
end

function MOI.eval_constraint_jacobian(data::NonlinearData, J, x)
    start = time()
    MOI.eval_constraint_jacobian(data.inner, J, x)
    data.eval_constraint_jacobian_timer += time() - start
    return
end

function MOI.hessian_lagrangian_structure(data::NonlinearData)
    return MOI.hessian_lagrangian_structure(data.inner)
end

function MOI.eval_hessian_lagrangian(data::NonlinearData, H, x, σ, μ)
    start = time()
    MOI.eval_hessian_lagrangian(data.inner, H, x, σ, μ)
    data.eval_hessian_lagrangian_timer += time() - start
    return
end

function MOI.eval_constraint_jacobian_product(data::NonlinearData, y, x, w)
    start = time()
    MOI.eval_constraint_jacobian_product(data.inner, y, x, w)
    data.eval_constraint_jacobian_timer += time() - start
    return
end

function MOI.eval_constraint_jacobian_transpose_product(
    data::NonlinearData,
    y,
    x,
    w,
)
    start = time()
    MOI.eval_constraint_jacobian_transpose_product(data.inner, y, x, w)
    data.eval_constraint_jacobian_timer += time() - start
    return
end

function MOI.eval_hessian_lagrangian_product(data::NonlinearData, H, x, v, σ, μ)
    start = time()
    MOI.eval_hessian_lagrangian_product(data.inner, H, x, v, σ, μ)
    data.eval_hessian_lagrangian_timer += time() - start
    return
end

"""
    adjacency_matrix(nodes::Vector{Node})

Compute the sparse adjacency matrix describing the parent-child relationships in
`nodes`.

The element `(i, j)` is `true` if there is an edge *from* `node[j]` to
`node[i]`. Since we get a column-oriented matrix, this gives us a fast way to
look up the edges leaving any node (i.e., the children).
"""
function adjacency_matrix(nodes::Vector{Node})
    N = length(nodes)
    I, J = Vector{Int}(undef, N), Vector{Int}(undef, N)
    numnz = 0
    for (i, node) in enumerate(nodes)
        if node.parent < 0
            continue
        end
        numnz += 1
        I[numnz] = i
        J[numnz] = node.parent
    end
    resize!(I, numnz)
    resize!(J, numnz)
    return SparseArrays.sparse(I, J, ones(Bool, numnz), N, N)
end

"""
    evaluate(
        f::AbstractDict,
        data::NonlinearData,
        index::ExpressionIndex,
    )

Evaluate the nonlinear expression `index`, where `f[x]` returns the primal value
of decision variable `x::MOI.VariableIndex`.
"""
function evaluate(
    f::AbstractDict,
    data::NonlinearData,
    index::ExpressionIndex;
    kwargs...,
)
    return evaluate(f, data, data[index]; kwargs...)
end

"""
    evaluate(
        f::AbstractDict,
        data::NonlinearData,
        expr::NonlinearExpression,
    )

Evaluate the nonlinear expression `expr`, where `f[x]` returns the primal value
of decision variable `x::MOI.VariableIndex`.
"""
function evaluate(
    f::AbstractDict,
    data::NonlinearData,
    expr::NonlinearExpression;
    evaluated_expressions = Dict{Int,Float64}(),
)
    storage = zeros(length(expr.nodes))
    adj = adjacency_matrix(expr.nodes)
    children_arr = SparseArrays.rowvals(adj)
    # An arbitrary limit on the potential input size of a multivariate
    # operation. This will get resized if need-be.
    input_cache = zeros(10)
    for k in length(expr.nodes):-1:1
        node = expr.nodes[k]
        if node.type == NODE_MOI_VARIABLE
            storage[k] = f[MOI.VariableIndex(node.index)]
        elseif node.type == NODE_VALUE
            storage[k] = expr.values[node.index]
        elseif node.type == NODE_SUBEXPRESSION
            if !haskey(evaluated_expressions, node.index)
                evaluated_expressions[node.index] = evaluate(
                    f,
                    data,
                    ExpressionIndex(node.index);
                    evaluated_expressions = evaluated_expressions,
                )
            end
            storage[k] = evaluated_expressions[node.index]
        elseif node.type == NODE_PARAMETER
            storage[k] = data.parameters[node.index]
        elseif node.type == NODE_CALL_MULTIVARIATE
            children_indices = SparseArrays.nzrange(adj, k)
            N = length(children_indices)
            if length(input_cache) < N
                resize!(input_cache, N)
            end
            f_input = view(input_cache, 1:N)
            for (r, i) in enumerate(children_indices)
                f_input[r] = storage[children_arr[i]]
            end
            storage[k] = eval_multivariate_function(
                data.operators,
                data.operators.multivariate_operators[node.index],
                f_input,
            )
        elseif node.type == NODE_CALL_UNIVARIATE
            child_idx = children_arr[adj.colptr[k]]
            storage[k] = eval_univariate_function(
                data.operators,
                data.operators.univariate_operators[node.index],
                storage[child_idx],
            )
        elseif node.type == NODE_COMPARISON
            children_idx = SparseArrays.nzrange(adj, k)
            result = true
            for r in 2:length(children_idx)
                lhs = children_arr[children_idx[r-1]]
                rhs = children_arr[children_idx[r]]
                result &= eval_comparison_function(
                    data.operators,
                    data.operators.comparison_operators[node.index],
                    storage[lhs],
                    storage[rhs],
                )
            end
            storage[k] = result
        else
            @assert node.type == NODE_LOGIC
            children_idx = SparseArrays.nzrange(adj, k)
            lhs = children_arr[children_idx[1]]
            rhs = children_arr[children_idx[2]]
            storage[k] = eval_logic_function(
                data.operators,
                data.operators.logic_operators[node.index],
                storage[lhs] == 1,
                storage[rhs] == 1,
            )
        end
    end
    return storage[1]
end

end  # module
