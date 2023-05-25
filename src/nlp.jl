#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    nonlinear_model(
        model::Model;
        force::Bool = false,
    )::Union{MOI.Nonlinear.Model,Nothing}

If `model` has nonlinear components, return a [`MOI.Nonlinear.Model`](@ref),
otherwise return `nothing`.

If `force`, always return a [`MOI.Nonlinear.Model`](@ref), and if one does not
exist for the model, create an empty one.
"""
function nonlinear_model(model::Model; force::Bool = false)
    if force
        _init_NLP(model)
    end
    return model.nlp_model
end

function _init_NLP(model::Model)
    if model.nlp_model === nothing
        model.nlp_model = MOI.Nonlinear.Model()
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

function MOI.Nonlinear.check_return_type(
    ::Type{T},
    ret::AbstractJuMPScalar,
) where {T}
    return error(
        "Expected return type of $T from a user-defined function, but got " *
        "$(typeof(ret)). Make sure your user-defined function only depends " *
        "on variables passed as arguments.",
    )
end

function MOI.Nonlinear.parse_expression(
    model::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::VariableRef,
    parent::Int,
)
    MOI.Nonlinear.parse_expression(model, expr, index(x), parent)
    return
end

function MOI.Nonlinear.parse_expression(
    model::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::GenericAffExpr,
    parent::Int,
)
    sum_id = model.operators.multivariate_operator_to_id[:+]
    prod_id = model.operators.multivariate_operator_to_id[:*]
    call_multi = MOI.Nonlinear.NODE_CALL_MULTIVARIATE
    n_terms = length(x.terms) + !iszero(x.constant)
    if n_terms == 0
        # If `x` is empty, then substitute `0.0` as the expression.
        MOI.Nonlinear.parse_expression(model, expr, 0.0, parent)
        return
    elseif n_terms == 1
        # If `x` contains 1 term, then we don't need the leading `+` operator.
    else
        push!(expr.nodes, MOI.Nonlinear.Node(call_multi, sum_id, parent))
        parent = length(expr.nodes)
    end
    if !iszero(x.constant)
        MOI.Nonlinear.parse_expression(model, expr, x.constant, parent)
    end
    for (v, c) in x.terms
        if isone(c)
            MOI.Nonlinear.parse_expression(model, expr, v, parent)
        else
            push!(expr.nodes, MOI.Nonlinear.Node(call_multi, prod_id, parent))
            mult_parent = length(expr.nodes)
            MOI.Nonlinear.parse_expression(model, expr, c, mult_parent)
            MOI.Nonlinear.parse_expression(model, expr, v, mult_parent)
        end
    end
    return
end

function MOI.Nonlinear.parse_expression(
    model::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::QuadExpr,
    parent::Int,
)
    sum_id = model.operators.multivariate_operator_to_id[:+]
    prod_id = model.operators.multivariate_operator_to_id[:*]
    call_multi = MOI.Nonlinear.NODE_CALL_MULTIVARIATE
    n_terms = length(x.terms) + !iszero(x.aff)
    if n_terms == 0
        # If `x` is empty, then substitute `0.0` as the expression.
        MOI.Nonlinear.parse_expression(model, expr, 0.0, parent)
        return
    elseif n_terms == 1
        # If `x` contains 1 term, then we don't need the leading `+` operator.
    else
        push!(expr.nodes, MOI.Nonlinear.Node(call_multi, sum_id, parent))
        parent = length(expr.nodes)
    end
    if !iszero(x.aff)
        MOI.Nonlinear.parse_expression(model, expr, x.aff, parent)
    end
    for (xy, c) in x.terms
        push!(expr.nodes, MOI.Nonlinear.Node(call_multi, prod_id, parent))
        mult_parent = length(expr.nodes)
        MOI.Nonlinear.parse_expression(model, expr, xy.a, mult_parent)
        MOI.Nonlinear.parse_expression(model, expr, xy.b, mult_parent)
        if !isone(c)
            MOI.Nonlinear.parse_expression(model, expr, c, mult_parent)
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

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> set_nonlinear_objective(model, MIN_SENSE, :(\$(x) + \$(x)^2))
```
"""
function set_nonlinear_objective(model::Model, sense::MOI.OptimizationSense, x)
    _init_NLP(model)
    set_objective_sense(model, sense)
    MOI.Nonlinear.set_objective(model.nlp_model, x)
    return
end

"""
    _nlp_objective_function(model::Model)

Returns the nonlinear objective function or `nothing` if no nonlinear objective
function is set.
"""
function _nlp_objective_function(model::Model)
    if model.nlp_model === nothing
        return nothing
    end
    return model.nlp_model.objective
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

function MOI.Nonlinear.parse_expression(
    model::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::NonlinearParameter,
    parent::Int,
)
    index = MOI.Nonlinear.ParameterIndex(x.index)
    return MOI.Nonlinear.parse_expression(model, expr, index, parent)
end

"""
    add_nonlinear_parameter(model::Model, value::Real)

Add an anonymous parameter to the model.
"""
function add_nonlinear_parameter(model::Model, value::Real)
    _init_NLP(model)
    p = MOI.Nonlinear.add_parameter(model.nlp_model, Float64(value))
    return NonlinearParameter(model, p.value)
end

"""
    index(p::NonlinearParameter)::MOI.Nonlinear.ParameterIndex

Return the index of the nonlinear parameter associated with `p`.
"""
index(p::NonlinearParameter) = MOI.Nonlinear.ParameterIndex(p.index)

"""
    value(p::NonlinearParameter)

Return the current value stored in the nonlinear parameter `p`.

## Example

```jldoctest
julia> model = Model();

julia> @NLparameter(model, p == 10)
p == 10.0

julia> value(p)
10.0
```
"""
function value(p::NonlinearParameter)
    return p.model.nlp_model[MOI.Nonlinear.ParameterIndex(p.index)]
end

"""
    set_value(p::NonlinearParameter, v::Number)

Store the value `v` in the nonlinear parameter `p`.

## Example

```jldoctest
julia> model = Model();

julia> @NLparameter(model, p == 0)
p == 0.0

julia> set_value(p, 5)
5

julia> value(p)
5.0
```
"""
function set_value(p::NonlinearParameter, value::Number)
    p.model.nlp_model[MOI.Nonlinear.ParameterIndex(p.index)] = value
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

function MOI.Nonlinear.parse_expression(
    model::MOI.Nonlinear.Model,
    expr::MOI.Nonlinear.Expression,
    x::NonlinearExpression,
    parent::Int,
)
    index = MOI.Nonlinear.ExpressionIndex(x.index)
    return MOI.Nonlinear.parse_expression(model, expr, index, parent)
end

"""
    index(ex::NonlinearExpression)::MOI.Nonlinear.ExpressionIndex

Return the index of the nonlinear expression associated with `ex`.
"""
index(ex::NonlinearExpression) = MOI.Nonlinear.ExpressionIndex(ex.index)

"""
    add_nonlinear_expression(model::Model, expr::Expr)

Add a nonlinear expression `expr` to `model`.

This function is most useful if the expression `expr` is generated
programmatically, and you cannot use [`@NLexpression`](@ref).

## Notes

 * You must interpolate the variables directly into the expression `expr`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> add_nonlinear_expression(model, :(\$(x) + \$(x)^2))
subexpression[1]: x + x ^ 2.0
```
"""
function add_nonlinear_expression(model::Model, ex)
    _init_NLP(model)
    index = MOI.Nonlinear.add_expression(model.nlp_model, ex)
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
    return MOI.Nonlinear.evaluate(
        _VariableValueMap(ex.model, var_value),
        ex.model.nlp_model,
        MOI.Nonlinear.ExpressionIndex(ex.index),
    )
end

###
### Nonlinear constraints
###

# For backwards compatibility, we need to retain the previously exported
# `NonlinearConstraintIndex`. They both have a `.value` field, but the new
# `Nonlinear.ConstraintIndex` is not guaranteed to be 1-based because some
# constraints may have been deleted.
const NonlinearConstraintIndex = MOI.Nonlinear.ConstraintIndex

"""
    NonlinearConstraintRef
"""
const NonlinearConstraintRef =
    ConstraintRef{Model,MOI.Nonlinear.ConstraintIndex}

function _normalize_constraint_expr(lhs::Real, body, rhs::Real)
    return Float64(lhs), body, Float64(rhs)
end

function _normalize_constraint_expr(lhs, body, rhs)
    return error(
        "Interval constraint contains non-constant left- or right-hand " *
        "sides. Reformulate as two separate constraints, or move all " *
        "variables into the central term.",
    )
end

_normalize_constraint_expr(lhs, rhs::Real) = lhs, Float64(rhs)

_normalize_constraint_expr(lhs, rhs) = Expr(:call, :-, lhs, rhs), 0.0

function _expr_to_constraint(expr::Expr)
    if isexpr(expr, :comparison)
        @assert expr.args[2] == expr.args[4]
        @assert expr.args[2] in (:<=, :>=)
        lhs, body, rhs =
            _normalize_constraint_expr(expr.args[1], expr.args[3], expr.args[5])
        return body, MOI.Interval(lhs, rhs)
    end
    lhs, rhs = _normalize_constraint_expr(expr.args[2], expr.args[3])
    if expr.args[1] == :<=
        return :($lhs - $rhs), MOI.LessThan(0.0)
    elseif expr.args[1] == :>=
        return :($lhs - $rhs), MOI.GreaterThan(0.0)
    else
        @assert expr.args[1] == :(==)
        return :($lhs - $rhs), MOI.EqualTo(0.0)
    end
end

"""
    add_nonlinear_constraint(model::Model, expr::Expr)

Add a nonlinear constraint described by the Julia expression `ex` to `model`.

This function is most useful if the expression `ex` is generated
programmatically, and you cannot use [`@NLconstraint`](@ref).

## Notes

 * You must interpolate the variables directly into the expression `expr`.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> add_nonlinear_constraint(model, :(\$(x) + \$(x)^2 <= 1))
(x + x ^ 2.0) - 1.0 ≤ 0
```
"""
function add_nonlinear_constraint(model::Model, ex::Expr)
    _init_NLP(model)
    f, set = _expr_to_constraint(ex)
    c = MOI.Nonlinear.add_constraint(model.nlp_model, f, set)
    return ConstraintRef(model, c, ScalarShape())
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
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    return MOI.is_valid(model.nlp_model, index)
end

"""
    delete(model::Model, c::NonlinearConstraintRef)

Delete the nonlinear constraint `c` from `model`.
"""
function delete(model::Model, c::NonlinearConstraintRef)
    _init_NLP(model)
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    MOI.Nonlinear.delete(model.nlp_model, index)
    return
end

"""
    num_nonlinear_constraints(model::Model)

Returns the number of nonlinear constraints associated with the `model`.
"""
function num_nonlinear_constraints(model::Model)
    nlp_model = nonlinear_model(model)
    if nlp_model === nothing
        return 0
    end
    return length(nlp_model.constraints)
end

"""
    all_nonlinear_constraints(model::Model)

Return a vector of all nonlinear constraint references in the model in the
order they were added to the model.
"""
function all_nonlinear_constraints(model::Model)
    nlp_model = nonlinear_model(model)
    if nlp_model === nothing
        return NonlinearConstraintRef[]
    end
    return NonlinearConstraintRef[
        ConstraintRef(model, index, ScalarShape()) for
        (index, _) in nlp_model.constraints
    ]
end

"""
    value(c::NonlinearConstraintRef; result::Int = 1)

Return the value of the `NonlinearConstraintRef` `c` associated with result
index `result` of the most-recent solution returned by the solver.

See also: [`result_count`](@ref).
"""
function value(c::NonlinearConstraintRef; result::Int = 1)
    return value(c) do x
        return value(x; result = result)
    end
end

"""
    value(var_value::Function, c::NonlinearConstraintRef)

Evaluate `c` using `var_value(v)` as the value for each variable `v`.
"""
function value(var_value::Function, c::NonlinearConstraintRef)
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    return MOI.Nonlinear.evaluate(
        _VariableValueMap(c.model, var_value),
        c.model.nlp_model,
        c.model.nlp_model[index].expression,
    )
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
    evaluator =
        MOI.get(c.model, MOI.NLPBlock()).evaluator::MOI.Nonlinear.Evaluator
    if length(evaluator.constraint_dual) == 0
        append!(evaluator.constraint_dual, MOI.get(c.model, MOI.NLPBlockDual()))
    end
    index = MOI.Nonlinear.ConstraintIndex(c.index.value)
    return evaluator.constraint_dual[MOI.Nonlinear.ordinal_index(
        evaluator,
        index,
    )]
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

## Example

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

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::T) where {T<:Real} = x^2
f (generic function with 1 method)

julia> register(model, :foo, 1, f; autodiff = true)

julia> @NLobjective(model, Min, foo(x))
```

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> g(x::T, y::T) where {T<:Real} = x * y
g (generic function with 1 method)

julia> register(model, :g, 2, g; autodiff = true)

julia> @NLobjective(model, Min, g(x[1], x[2]))
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
    MOI.Nonlinear.register_operator(model.nlp_model, op, dimension, f)
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

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::T) where {T<:Real} = x^2
f (generic function with 1 method)

julia> ∇f(x::T) where {T<:Real} = 2 * x
∇f (generic function with 1 method)

julia> register(model, :foo, 1, f, ∇f; autodiff = true)

julia> @NLobjective(model, Min, foo(x))
```

```jldoctest
julia> model = Model();

julia> @variable(model, x[1:2])
2-element Vector{VariableRef}:
 x[1]
 x[2]

julia> g(x::T, y::T) where {T<:Real} = x * y
g (generic function with 1 method)

julia> function ∇g(g::AbstractVector{T}, x::T, y::T) where {T<:Real}
           g[1] = y
           g[2] = x
           return
       end
∇g (generic function with 1 method)

julia> register(model, :g, 2, g, ∇g)

julia> @NLobjective(model, Min, g(x[1], x[2]))
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
        MOI.Nonlinear.register_operator(model.nlp_model, op, dimension, f, ∇f)
    else
        if autodiff == true
            @warn("autodiff = true ignored since gradient is already provided.")
        end
        MOI.Nonlinear.register_operator(model.nlp_model, op, dimension, f, ∇f)
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

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> f(x::Float64) = x^2
f (generic function with 1 method)

julia> ∇f(x::Float64) = 2 * x
∇f (generic function with 1 method)

julia> ∇²f(x::Float64) = 2.0
∇²f (generic function with 1 method)

julia> register(model, :foo, 1, f, ∇f, ∇²f)

julia> @NLobjective(model, Min, foo(x))

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
    _init_NLP(model)
    MOI.Nonlinear.register_operator(model.nlp_model, op, dimension, f, ∇f, ∇²f)
    return
end

"""
    NLPEvaluator(
        model::Model,
        _differentiation_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation =
            MOI.Nonlinear.SparseReverseMode(),
    )

Return an [`MOI.AbstractNLPEvaluator`](@ref) constructed from `model`

!!! warning
    Before using, you must initialize the evaluator using
    [`MOI.initialize`](@ref).

## Experimental

These features may change or be removed in any future version of JuMP.

Pass `_differentiation_backend` to specify the differentiation backend used to
compute derivatives.
"""
function NLPEvaluator(
    model::Model;
    _differentiation_backend::MOI.Nonlinear.AbstractAutomaticDifferentiation = MOI.Nonlinear.SparseReverseMode(),
)
    _init_NLP(model)
    return MOI.Nonlinear.Evaluator(
        model.nlp_model,
        _differentiation_backend,
        index.(all_variables(model)),
    )
end
