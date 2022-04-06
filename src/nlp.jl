#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function _init_NLP(model::Model)
    if model.nlp_data === nothing
        model.nlp_data = Nonlinear.NonlinearData()
    end
    return
end

function _init_NLP(m::AbstractModel)
    return error(
        "Encountered an error parsing nonlinear expression: we don't support " *
        "models of type $(typeof(m)). In general, JuMP's nonlinear features " *
        "don't work with JuMP-extensions.",
    )
end

function Nonlinear.check_return_type(
    ::Type{T},
    ret::AbstractJuMPScalar,
) where {T}
    return error(
        "Expected return type of $T from a user-defined function, but got " *
        "$(typeof(ret)). Make sure your user-defined function only depends " *
        "on variables passed as arguments.",
    )
end

function Nonlinear._parse_expression(
    data::Nonlinear.NonlinearData,
    expr::Nonlinear.NonlinearExpression,
    x::VariableRef,
    parent::Int,
)
    Nonlinear._parse_expression(data, expr, index(x), parent)
    return
end

function Nonlinear._parse_expression(
    data::Nonlinear.NonlinearData,
    expr::Nonlinear.NonlinearExpression,
    x::GenericAffExpr,
    parent::Int,
)
    push!(
        expr.nodes,
        Nonlinear.Node(
            Nonlinear.NODE_CALL_MULTIVARIATE,
            data.operators.multivariate_operator_to_id[:+],
            parent,
        ),
    )
    sum_parent = length(expr.nodes)
    if !iszero(x.constant)
        Nonlinear._parse_expression(data, expr, x.constant, sum_parent)
    end
    for (v, c) in x.terms
        if isone(c)
            Nonlinear._parse_expression(data, expr, v, sum_parent)
        else
            push!(
                expr.nodes,
                Nonlinear.Node(
                    Nonlinear.NODE_CALL_MULTIVARIATE,
                    data.operators.multivariate_operator_to_id[:*],
                    sum_parent,
                ),
            )
            mult_parent = length(expr.nodes)
            Nonlinear._parse_expression(data, expr, c, mult_parent)
            Nonlinear._parse_expression(data, expr, v, mult_parent)
        end
    end
    return
end

function Nonlinear._parse_expression(
    data::Nonlinear.NonlinearData,
    expr::Nonlinear.NonlinearExpression,
    x::QuadExpr,
    parent::Int,
)
    sum = data.operators.multivariate_operator_to_id[:+]
    prod = data.operators.multivariate_operator_to_id[:*]
    push!(
        expr.nodes,
        Nonlinear.Node(Nonlinear.NODE_CALL_MULTIVARIATE, sum, parent),
    )
    sum_parent = length(expr.nodes)
    Nonlinear._parse_expression(data, expr, x.aff, sum_parent)
    for (xy, c) in x.terms
        push!(
            expr.nodes,
            Nonlinear.Node(Nonlinear.NODE_CALL_MULTIVARIATE, prod, sum_parent),
        )
        mult_parent = length(expr.nodes)
        Nonlinear._parse_expression(data, expr, xy.a, mult_parent)
        Nonlinear._parse_expression(data, expr, xy.b, mult_parent)
        if !isone(c)
            Nonlinear._parse_expression(data, expr, c, mult_parent)
        end
    end
    return
end

###
### Nonlinear objectives
###

"""
    set_nonlinear_objective(
        model::Model,
        sense::MOI.OptimizationSense,
        expr::Expr,
    )

Set the nonlinear objective of `model` to the expression `expr`, with the
optimization sense `sense`.

This function is most useful if the expression `expr` is generated
programmatically, and you cannot use [`@NLobjective`](@ref).

## Notes

 * You must interpolate the variables directly into the expression `expr`.
 * You must use `MIN_SENSE` or `MAX_SENSE` instead of `Min` and `Max`.

## Examples

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> set_nonlinear_objective(model, MIN_SENSE, :(\$(x) + \$(x)^2))
```
"""
function set_nonlinear_objective(model::Model, sense::MOI.OptimizationSense, x)
    _init_NLP(model)
    set_objective_sense(model, sense)
    Nonlinear.set_objective(model.nlp_data, x)
    return
end

"""
    _nlp_objective_function(model::Model)

Returns the nonlinear objective function or `nothing` if no nonlinear objective
function is set.
"""
function _nlp_objective_function(model::Model)
    if model.nlp_data === nothing
        return nothing
    end
    return model.nlp_data.objective
end

###
### Nonlinear parameters
###

"""
    NonlinearParameter <: AbstractJuMPScalar

A struct to represent a nonlinear parameter.

Create a parameter using [`@NLparameter`](@ref).
"""
struct NonlinearParameter <: AbstractJuMPScalar
    model::Model
    index::Int
end

function Nonlinear._parse_expression(
    data::Nonlinear.NonlinearData,
    expr::Nonlinear.NonlinearExpression,
    x::NonlinearParameter,
    parent::Int,
)
    index = Nonlinear.ParameterIndex(x.index)
    return Nonlinear._parse_expression(data, expr, index, parent)
end

function add_nonlinear_parameter(model::Model, value::Real)
    _init_NLP(model)
    p = Nonlinear.add_parameter(model.nlp_data, Float64(value))
    return NonlinearParameter(model, p.value)
end

"""
    value(p::NonlinearParameter)

Return the current value stored in the nonlinear parameter `p`.

# Example

```jldoctest; setup=:(using JuMP)
model = Model()
@NLparameter(model, p == 10)
value(p)

# output
10.0
```
"""
function value(p::NonlinearParameter)
    return p.model.nlp_data[Nonlinear.ParameterIndex(p.index)]
end

"""
    set_value(p::NonlinearParameter, v::Number)

Store the value `v` in the nonlinear parameter `p`.

# Example
```jldoctest; setup=:(using JuMP)
model = Model()
@NLparameter(model, p == 0)
set_value(p, 5)
value(p)

# output
5.0
```
"""
function set_value(p::NonlinearParameter, value::Number)
    p.model.nlp_data[Nonlinear.ParameterIndex(p.index)] = value
    return value
end

###
### Nonlinear expressions
###

"""
    NonlinearExpression <: AbstractJuMPScalar

A struct to represent a nonlinear expression.

Create an expression using [`@NLexpression`](@ref).
"""
struct NonlinearExpression <: AbstractJuMPScalar
    model::Model
    index::Int
end

function Nonlinear._parse_expression(
    data::Nonlinear.NonlinearData,
    expr::Nonlinear.NonlinearExpression,
    x::NonlinearExpression,
    parent::Int,
)
    index = Nonlinear.ExpressionIndex(x.index)
    return Nonlinear._parse_expression(data, expr, index, parent)
end

"""
    add_nonlinear_expression(model::Model, expr::Expr)

Add a nonlinear expression `expr` to `model`.

This function is most useful if the expression `expr` is generated
programmatically, and you cannot use [`@NLexpression`](@ref).

## Notes

 * You must interpolate the variables directly into the expression `expr`.

## Examples

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> add_nonlinear_expression(model, :(\$(x) + \$(x)^2))
subexpression[1]: x + x ^ 2.0
```
"""
function add_nonlinear_expression(model::Model, ex)
    _init_NLP(model)
    index = Nonlinear.add_expression(model.nlp_data, ex)
    return NonlinearExpression(model, index.value)
end

"""
    value(ex::NonlinearExpression; result::Int = 1)

Return the value of the `NonlinearExpression` `ex` associated with result index
`result` of the most-recent solution returned by the solver.

Replaces `getvalue` for most use cases.

See also: [`result_count`](@ref).
"""
function value(ex::NonlinearExpression; result::Int = 1)
    return value(ex) do x
        return value(x; result = result)
    end
end

"""
    _VariableValueMap{T,F}

A lazy cache used for computing the primal variable solution in `value`.

This avoids the need to rewrite the nonlinear expressions from MOI_VARIABLE to
VARIABLE, as well as eagerly computing the `var_value` for every variable. We
use a `cache` so we don't have to recompute variables we have already seen.
"""
struct _VariableValueMap{F} <: AbstractDict{MOI.VariableIndex,Float64}
    model::Model
    value::F
    cache::Dict{MOI.VariableIndex,Float64}
    function _VariableValueMap(model, value::F) where {F}
        return new{F}(model, value, Dict{MOI.VariableIndex,Float64}())
    end
end

function Base.getindex(map::_VariableValueMap, index::MOI.VariableIndex)
    return get!(map.cache, index) do
        return map.value(VariableRef(map.model, index))
    end
end

"""
    value(var_value::Function, ex::NonlinearExpression)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(var_value::Function, ex::NonlinearExpression)
    return Nonlinear.evaluate(
        _VariableValueMap(ex.model, var_value),
        ex.model.nlp_data,
        Nonlinear.ExpressionIndex(ex.index),
    )
end

###
### Nonlinear constraints
###

"""
    NonlinearConstraintIndex(index::Int64)

A struct to refer to the 1-indexed nonlinear constraint `index`.
"""
struct NonlinearConstraintIndex
    value::Int64
end

const NonlinearConstraintRef = ConstraintRef{Model,NonlinearConstraintIndex}

"""
    add_nonlinear_constraint(model::Model, expr::Expr)

Add a nonlinear constraint described by the Julia expression `ex` to `model`.

This function is most useful if the expression `ex` is generated
programmatically, and you cannot use [`@NLconstraint`](@ref).

## Notes

 * You must interpolate the variables directly into the expression `expr`.

## Examples

```jldoctest; setup=:(using JuMP; model = Model(); @variable(model, x))
julia> add_nonlinear_constraint(model, :(\$(x) + \$(x)^2 <= 1))
(x + x ^ 2.0) - 1.0 ≤ 0
```
"""
function add_nonlinear_constraint(model::Model, ex::Expr)
    _init_NLP(model)
    c = Nonlinear.add_constraint(model.nlp_data, ex)
    index = NonlinearConstraintIndex(c.value)
    return ConstraintRef(model, index, ScalarShape())
end

"""
    is_valid(model::Model, c::NonlinearConstraintRef)

Return `true` if `c` refers to a valid nonlinear constraint in `model`.
"""
function is_valid(model::Model, c::NonlinearConstraintRef)
    if model !== c.model
        return false
    end
    _init_NLP(model)
    index = Nonlinear.ConstraintIndex(c.index.value)
    return MOI.is_valid(model.nlp_data, index)
end

"""
    num_nonlinear_constraints(model::Model)

Returns the number of nonlinear constraints associated with the `model`.
"""
function num_nonlinear_constraints(model::Model)
    return model.nlp_data !== nothing ? length(model.nlp_data.constraints) : 0
end

"""
    all_nonlinear_constraints(model::Model)

Return a vector of all nonlinear constraint references in the model in the
order they were added to the model.
"""
function all_nonlinear_constraints(model::Model)
    return map(1:num_nonlinear_constraints(model)) do i
        return ConstraintRef(model, NonlinearConstraintIndex(i), ScalarShape())
    end
end

###
### Nonlinear dual solutions
###

"""
    dual(c::NonlinearConstraintRef)

Return the dual of the nonlinear constraint `c`.
"""
function dual(c::NonlinearConstraintRef)
    _init_NLP(c.model)
    dual = c.model.nlp_data.constraint_dual
    if length(dual) == 0
        append!(dual, MOI.get(c.model, MOI.NLPBlockDual()))
    end
    index = Nonlinear.ConstraintIndex(c.index.value)
    return dual[Nonlinear.row(c.model.nlp_data, index)]
end

"""
    nonlinear_dual_start_value(model::Model)

Return the current value of the MOI attribute [`MOI.NLPBlockDualStart`](@ref).
"""
function nonlinear_dual_start_value(model::Model)
    return MOI.get(model, MOI.NLPBlockDualStart())
end

"""
    set_nonlinear_dual_start_value(
        model::Model,
        start::Union{Nothing,Vector{Float64}},
    )

Set the value of the MOI attribute [`MOI.NLPBlockDualStart`](@ref).

The start vector corresponds to the Lagrangian duals of the nonlinear
constraints, in the order given by [`all_nonlinear_constraints`](@ref). That is, you
must pass a single start vector corresponding to all of the nonlinear
constraints in a single function call; you cannot set the dual start value of
nonlinear constraints one-by-one. The example below demonstrates how to use
[`all_nonlinear_constraints`](@ref) to create a mapping between the nonlinear
constraint references and the start vector.

Pass `nothing` to unset a previous start.

## Examples

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2]);

julia> nl1 = @NLconstraint(model, x[1] <= sqrt(x[2]));

julia> nl2 = @NLconstraint(model, x[1] >= exp(x[2]));

julia> start = Dict(nl1 => -1.0, nl2 => 1.0);

julia> start_vector = [start[con] for con in all_nonlinear_constraints(model)]
2-element Vector{Float64}:
 -1.0
  1.0

julia> set_nonlinear_dual_start_value(model, start_vector)

julia> nonlinear_dual_start_value(model)
2-element Vector{Float64}:
 -1.0
  1.0
```
"""
function set_nonlinear_dual_start_value(model::Model, start::Vector{Float64})
    _init_NLP(model)
    N = num_nonlinear_constraints(model)
    if length(start) != N
        throw(
            ArgumentError(
                "length start vector ($(length(start))) does not match the " *
                "number of nonlinear constraints ($N).",
            ),
        )
    end
    MOI.set(model, MOI.NLPBlockDualStart(), start)
    return
end

function set_nonlinear_dual_start_value(model::Model, start::Nothing)
    MOI.set(model, MOI.NLPBlockDualStart(), start)
    return
end

"""
    register(
        model::Model,
        op::Symbol,
        dimension::Integer,
        f::Function;
        autodiff:Bool = false,
    )

Register the user-defined function `f` that takes `dimension` arguments in
`model` as the symbol `op`.

The function `f` must support all subtypes of `Real` as arguments. Do not assume
that the inputs are `Float64`.

## Notes

 * For this method, you must explicitly set `autodiff = true`, because no
   user-provided gradient function `∇f` is given.
 * Second-derivative information is only computed if `dimension == 1`.
 * `op` does not have to be the same symbol as `f`, but it is generally more
   readable if it is.

## Examples

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x)
f(x::T) where {T<:Real} = x^2
register(model, :foo, 1, f; autodiff = true)
@NLobjective(model, Min, foo(x))
```

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x[1:2])
g(x::T, y::T) where {T<:Real} = x * y
register(model, :g, 2, g; autodiff = true)
@NLobjective(model, Min, g(x[1], x[2]))
```
"""
function register(
    model::Model,
    op::Symbol,
    dimension::Integer,
    f::Function;
    autodiff::Bool = false,
)
    if autodiff == false
        error("If only the function is provided, must set autodiff=true")
    end
    _init_NLP(model)
    Nonlinear.register_operator(model.nlp_data, op, dimension, f)
    return
end

"""
    register(
        model::Model,
        s::Symbol,
        dimension::Integer,
        f::Function,
        ∇f::Function;
        autodiff:Bool = false,
    )

Register the user-defined function `f` that takes `dimension` arguments in
`model` as the symbol `s`. In addition, provide a gradient function `∇f`.

The functions `f`and `∇f` must support all subtypes of `Real` as arguments. Do
not assume that the inputs are `Float64`.

## Notes

 * If the function `f` is univariate (i.e., `dimension == 1`), `∇f` must return
   a number which represents the first-order derivative of the function `f`.
 * If the function `f` is multi-variate, `∇f` must have a signature matching
   `∇f(g::AbstractVector{T}, args::T...) where {T<:Real}`, where the first
   argument is a vector `g` that is modified in-place with the gradient.
 * If `autodiff = true` and `dimension == 1`, use automatic differentiation to
   compute the second-order derivative information. If `autodiff = false`, only
   first-order derivative information will be used.
 * `s` does not have to be the same symbol as `f`, but it is generally more
   readable if it is.

## Examples

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x)
f(x::T) where {T<:Real} = x^2
∇f(x::T) where {T<:Real} = 2 * x
register(model, :foo, 1, f, ∇f; autodiff = true)
@NLobjective(model, Min, foo(x))
```

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x[1:2])
g(x::T, y::T) where {T<:Real} = x * y
function ∇g(g::AbstractVector{T}, x::T, y::T) where {T<:Real}
    g[1] = y
    g[2] = x
    return
end
register(model, :g, 2, g, ∇g; autodiff = true)
@NLobjective(model, Min, g(x[1], x[2]))
```
"""
function register(
    model::Model,
    op::Symbol,
    dimension::Integer,
    f::Function,
    ∇f::Function;
    autodiff::Bool = false,
)
    _init_NLP(model)
    if dimension == 1
        if autodiff == false
            error(
                "Currently must provide 2nd order derivatives of univariate functions. Try setting autodiff=true.",
            )
        end
        Nonlinear.register_operator(model.nlp_data, op, dimension, f, ∇f)
    else
        if autodiff == true
            @warn("autodiff = true ignored since gradient is already provided.")
        end
        Nonlinear.register_operator(model.nlp_data, op, dimension, f, ∇f)
    end
    return
end

"""
    register(
        model::Model,
        s::Symbol,
        dimension::Integer,
        f::Function,
        ∇f::Function,
        ∇²f::Function,
    )

Register the user-defined function `f` that takes `dimension` arguments in
`model` as the symbol `s`. In addition, provide a gradient function `∇f` and a
hessian function `∇²f`.

`∇f` and `∇²f` must return numbers corresponding to the first- and second-order
derivatives of the function `f` respectively.

## Notes

 * Because automatic differentiation is not used, you can assume the inputs are
   all `Float64`.
 * This method will throw an error if `dimension > 1`.
 * `s` does not have to be the same symbol as `f`, but it is generally more
   readable if it is.

## Examples

```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x)
f(x::Float64) = x^2
∇f(x::Float64) = 2 * x
∇²f(x::Float64) = 2.0
register(model, :foo, 1, f, ∇f, ∇²f)
@NLobjective(model, Min, foo(x))
```
"""
function register(
    model::Model,
    op::Symbol,
    dimension::Integer,
    f::Function,
    ∇f::Function,
    ∇²f::Function,
)
    if dimension > 1
        error(
            "Providing hessians for multivariate functions is not yet supported",
        )
    end
    _init_NLP(model)
    Nonlinear.register_operator(model.nlp_data, op, dimension, f, ∇f, ∇²f)
    return
end

"""
    NLPEvaluator(model)

Return an [`MOI.AbstractNLPEvaluator`](@ref) constructed from `model`.

Before using, you must initialize the evaluator using [`MOI.initialize`](@ref).
"""
function NLPEvaluator(model::Model)
    _init_NLP(model)
    variables = all_variables(model)
    Nonlinear.set_differentiation_backend(
        model.nlp_data,
        Nonlinear.SparseReverseMode(),
        index.(variables),
    )
    return model.nlp_data
end
