#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

mutable struct _NonlinearExprData
    nd::Vector{NodeData}
    const_values::Vector{Float64}
end

function set_objective(
    m::Model,
    sense::MOI.OptimizationSense,
    ex::_NonlinearExprData,
)
    _init_NLP(m)
    set_objective_sense(m, sense)
    m.nlp_data.nlobj = ex
    return
end

include("parse_nlp.jl")

# GenericRangeConstraint
# l ≤ ∑ aᵢ xᵢ ≤ u
# The constant part of the internal expression is assumed to be zero
mutable struct _NonlinearConstraint <: AbstractConstraint
    terms::_NonlinearExprData
    lb::Float64
    ub::Float64
end

"""
    NonlinearConstraintIndex(index::Int64)

A struct to refer to the 1-indexed nonlinear constraint `index`.
"""
struct NonlinearConstraintIndex
    value::Int64
end

#  b ≤ expr ≤ b   →   ==
# -∞ ≤ expr ≤ u   →   <=
#  l ≤ expr ≤ ∞   →   >=
#  l ≤ expr ≤ u   →   range
function _sense(c::_NonlinearConstraint)
    if c.lb != -Inf
        if c.ub != Inf
            if c.ub == c.lb
                return :(==)
            else
                return :range
            end
        else
            return :(>=)
        end
    else #if c.lb == -Inf
        c.ub == Inf && error("'Free' constraint sense not supported")
        return :(<=)
    end
end

function _rhs(c::_NonlinearConstraint)
    s = _sense(c)
    s == :range && error("Range constraints do not have a well-defined RHS")
    return s == :(<=) ? c.ub : c.lb
end

mutable struct _NLPData
    nlobj::Union{Nothing,_NonlinearExprData}
    nlconstr::Vector{_NonlinearConstraint}
    nlexpr::Vector{_NonlinearExprData}
    nlconstr_duals::Vector{Float64}
    nlparamvalues::Vector{Float64}
    user_operators::_Derivatives.UserOperatorRegistry
    largest_user_input_dimension::Int
    evaluator::Any
end

"""
    _nlp_objective_function(model::Model)

Returns the nonlinear objective function or `nothing` if no nonlinear objective
function is set.
"""
function _nlp_objective_function(model::Model)
    if model.nlp_data === nothing
        return nothing
    else
        return model.nlp_data.nlobj
    end
end

function _create_nlp_block_data(m::Model)
    @assert m.nlp_data !== nothing
    bounds = MOI.NLPBoundsPair[]
    for constr in m.nlp_data.nlconstr
        push!(bounds, MOI.NLPBoundsPair(constr.lb, constr.ub))
    end
    return MOI.NLPBlockData(
        bounds,
        NLPEvaluator(m),
        isa(m.nlp_data.nlobj, _NonlinearExprData),
    )
end

function NonlinearExpression(m::Model, ex::_NonlinearExprData)
    _init_NLP(m)
    nldata::_NLPData = m.nlp_data
    push!(nldata.nlexpr, ex)
    return NonlinearExpression(m, length(nldata.nlexpr))
end

function _new_parameter(m::Model, value::Number)
    _init_NLP(m)
    nldata::_NLPData = m.nlp_data
    push!(nldata.nlparamvalues, value)
    return NonlinearParameter(m, length(nldata.nlparamvalues))
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
value(p::NonlinearParameter) = p.model.nlp_data.nlparamvalues[p.index]::Float64

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
function set_value(p::NonlinearParameter, v::Number)
    return p.model.nlp_data.nlparamvalues[p.index] = v
end

function _NLPData()
    return _NLPData(
        nothing,
        _NonlinearConstraint[],
        _NonlinearExprData[],
        Float64[],
        Float64[],
        _Derivatives.UserOperatorRegistry(),
        0,
        nothing,
    )
end

Base.copy(::_NLPData) = error("Copying nonlinear problems not yet implemented")

function _init_NLP(m::Model)
    if m.nlp_data === nothing
        m.nlp_data = _NLPData()
    end
end

const NonlinearConstraintRef = ConstraintRef{Model,NonlinearConstraintIndex}

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

function all_nl_constraints(::Model)
    return error(
        "`all_nl_constraints`` has been renamed `all_nonlinear_constraints`",
    )
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
    return 1 <= c.index.value <= num_nonlinear_constraints(model)
end

"""
    dual(c::NonlinearConstraintRef)

Return the dual of the nonlinear constraint `c`.
"""
function dual(c::NonlinearConstraintRef)
    _init_NLP(c.model)
    nldata::_NLPData = c.model.nlp_data
    # The array is cleared on every solve.
    if length(nldata.nlconstr_duals) != length(nldata.nlconstr)
        nldata.nlconstr_duals = MOI.get(c.model, MOI.NLPBlockDual())
    end
    return nldata.nlconstr_duals[c.index.value]
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
    nldata::_NLPData = model.nlp_data
    if length(nldata.nlconstr) != length(start)
        throw(
            ArgumentError(
                "length start vector ($(length(start))) does not match the " *
                "number of nonlinear constraints ($(length(nldata.nlconstr))).",
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

mutable struct _FunctionStorage
    nd::Vector{NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    grad_sparsity::Vector{Int}
    hess_I::Vector{Int} # nonzero pattern of hessian
    hess_J::Vector{Int}
    rinfo::Coloring.RecoveryInfo # coloring info for hessians
    seed_matrix::Matrix{Float64}
    linearity::Linearity
    dependent_subexpressions::Vector{Int} # subexpressions which this function depends on, ordered for forward pass
end

mutable struct _SubexpressionStorage
    nd::Vector{NodeData}
    adj::SparseMatrixCSC{Bool,Int}
    const_values::Vector{Float64}
    forward_storage::Vector{Float64}
    partials_storage::Vector{Float64}
    reverse_storage::Vector{Float64}
    forward_storage_ϵ::Vector{Float64}
    partials_storage_ϵ::Vector{Float64}
    reverse_storage_ϵ::Vector{Float64}
    linearity::Linearity
end

"""
    NLPEvaluator(m::Model)

Return an `MOI.AbstractNLPEvaluator` constructed from the model `model`.

Before using, you must initialize the evaluator using `MOI.initialize`.
"""
mutable struct NLPEvaluator <: MOI.AbstractNLPEvaluator
    model::Model
    parameter_values::Vector{Float64}
    has_nlobj::Bool
    objective::_FunctionStorage
    constraints::Vector{_FunctionStorage}
    subexpressions::Vector{_SubexpressionStorage}
    subexpression_order::Vector{Int}
    subexpression_forward_values::Vector{Float64}
    subexpression_reverse_values::Vector{Float64}
    subexpression_linearity::Vector{_Derivatives.Linearity}
    subexpressions_as_julia_expressions::Vector{Any}
    last_x::Vector{Float64}
    jac_storage::Vector{Float64} # temporary storage for computing jacobians
    user_output_buffer::Vector{Float64} # temporary storage for user-defined functions
    # storage for computing hessians
    # these Float64 vectors are reinterpreted to hold multiple epsilon components
    # so the length should be multiplied by the maximum number of epsilon components
    disable_2ndorder::Bool # don't offer Hess or HessVec
    want_hess::Bool
    forward_storage_ϵ::Vector{Float64} # (longest expression)
    partials_storage_ϵ::Vector{Float64} # (longest expression)
    reverse_storage_ϵ::Vector{Float64} # (longest expression)
    input_ϵ::Vector{Float64} # (number of variables)
    output_ϵ::Vector{Float64}# (number of variables)
    subexpression_forward_values_ϵ::Vector{Float64} # (number of subexpressions)
    subexpression_reverse_values_ϵ::Vector{Float64} # (number of subexpressions)
    hessian_sparsity::Vector{Tuple{Int64,Int64}}
    max_chunk::Int # chunk size for which we've allocated storage
    # timers
    eval_objective_timer::Float64
    eval_constraint_timer::Float64
    eval_objective_gradient_timer::Float64
    eval_constraint_jacobian_timer::Float64
    eval_hessian_lagrangian_timer::Float64

    function NLPEvaluator(m::Model)
        d = new(m)
        d.eval_objective_timer = 0
        d.eval_constraint_timer = 0
        d.eval_objective_gradient_timer = 0
        d.eval_constraint_jacobian_timer = 0
        d.eval_hessian_lagrangian_timer = 0
        return d
    end
end

function _replace_moi_variables(
    nd::Vector{NodeData},
    moi_index_to_consecutive_index,
)
    new_nd = Vector{NodeData}(undef, length(nd))
    for i in 1:length(nd)
        node = nd[i]
        if node.nodetype == MOIVARIABLE
            new_nd[i] = NodeData(
                VARIABLE,
                moi_index_to_consecutive_index[_MOIVAR(node.index)],
                node.parent,
            )
        else
            new_nd[i] = node
        end
    end
    return new_nd
end

function _FunctionStorage(
    nd::Vector{NodeData},
    const_values,
    num_variables,
    coloring_storage,
    want_hess::Bool,
    subexpressions::Vector{_SubexpressionStorage},
    dependent_subexpressions,
    subexpression_linearity,
    subexpression_edgelist,
    subexpression_variables,
    moi_index_to_consecutive_index,
)
    nd = _replace_moi_variables(nd, moi_index_to_consecutive_index)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    empty!(coloring_storage)
    compute_gradient_sparsity!(coloring_storage, nd)

    for k in dependent_subexpressions
        compute_gradient_sparsity!(coloring_storage, subexpressions[k].nd)
    end
    grad_sparsity = sort!(collect(coloring_storage))
    empty!(coloring_storage)

    if want_hess
        # compute hessian sparsity
        linearity = classify_linearity(nd, adj, subexpression_linearity)
        edgelist = compute_hessian_sparsity(
            nd,
            adj,
            linearity,
            coloring_storage,
            subexpression_edgelist,
            subexpression_variables,
        )
        hess_I, hess_J, rinfo = Coloring.hessian_color_preprocess(
            edgelist,
            num_variables,
            coloring_storage,
        )
        seed_matrix = Coloring.seed_matrix(rinfo)
    else
        hess_I = hess_J = Int[]
        rinfo = Coloring.RecoveryInfo()
        seed_matrix = Array{Float64}(undef, 0, 0)
        linearity = [NONLINEAR]
    end

    return _FunctionStorage(
        nd,
        adj,
        const_values,
        forward_storage,
        partials_storage,
        reverse_storage,
        grad_sparsity,
        hess_I,
        hess_J,
        rinfo,
        seed_matrix,
        linearity[1],
        dependent_subexpressions,
    )
end

function _SubexpressionStorage(
    nd::Vector{NodeData},
    const_values,
    num_variables,
    subexpression_linearity,
    moi_index_to_consecutive_index,
)
    nd = _replace_moi_variables(nd, moi_index_to_consecutive_index)
    adj = adjmat(nd)
    forward_storage = zeros(length(nd))
    partials_storage = zeros(length(nd))
    reverse_storage = zeros(length(nd))
    linearity = classify_linearity(nd, adj, subexpression_linearity)

    empty_arr = Array{Float64}(undef, 0)

    return _SubexpressionStorage(
        nd,
        adj,
        const_values,
        forward_storage,
        partials_storage,
        reverse_storage,
        empty_arr,
        empty_arr,
        empty_arr,
        linearity[1],
    )
end

function MOI.initialize(d::NLPEvaluator, requested_features::Vector{Symbol})
    nldata::_NLPData = d.model.nlp_data

    for feat in requested_features
        if !(feat in MOI.features_available(d))
            error("Unsupported feature $feat")
        end
    end
    if d.eval_objective_timer != 0
        # we've already been initialized
        # assume no new features are being requested.
        return
    end

    num_variables_ = num_variables(d.model)

    moi_index_to_consecutive_index = Dict(
        moi_index => consecutive_index for (consecutive_index, moi_index) in
        enumerate(MOI.get(d.model, MOI.ListOfVariableIndices()))
    )

    d.user_output_buffer =
        Array{Float64}(undef, d.model.nlp_data.largest_user_input_dimension)
    d.jac_storage = Array{Float64}(
        undef,
        max(num_variables_, d.model.nlp_data.largest_user_input_dimension),
    )

    d.constraints = _FunctionStorage[]
    d.last_x = fill(NaN, num_variables_)

    d.parameter_values = nldata.nlparamvalues

    d.want_hess = (:Hess in requested_features)
    want_hess_storage = (:HessVec in requested_features) || d.want_hess
    coloring_storage = _Derivatives.Coloring.IndexedSet(num_variables_)

    d.has_nlobj = nldata.nlobj !== nothing
    max_expr_length = 0
    main_expressions = Array{Vector{NodeData}}(undef, 0)
    subexpr = Array{Vector{NodeData}}(undef, 0)
    for nlexpr in nldata.nlexpr
        push!(subexpr, nlexpr.nd)
    end
    if d.has_nlobj
        push!(main_expressions, nldata.nlobj.nd)
    end
    for nlconstr in nldata.nlconstr
        push!(main_expressions, nlconstr.terms.nd)
    end
    d.subexpression_order, individual_order =
        order_subexpressions(main_expressions, subexpr)

    d.subexpression_linearity = Array{Linearity}(undef, length(nldata.nlexpr))
    subexpression_variables = Array{Vector{Int}}(undef, length(nldata.nlexpr))
    subexpression_edgelist =
        Array{Set{Tuple{Int,Int}}}(undef, length(nldata.nlexpr))
    d.subexpressions =
        Array{_SubexpressionStorage}(undef, length(nldata.nlexpr))
    d.subexpression_forward_values =
        Array{Float64}(undef, length(d.subexpressions))
    d.subexpression_reverse_values =
        Array{Float64}(undef, length(d.subexpressions))

    empty_edgelist = Set{Tuple{Int,Int}}()
    for k in d.subexpression_order # only load expressions which actually are used
        d.subexpression_forward_values[k] = NaN
        d.subexpressions[k] = _SubexpressionStorage(
            nldata.nlexpr[k].nd,
            nldata.nlexpr[k].const_values,
            num_variables_,
            d.subexpression_linearity,
            moi_index_to_consecutive_index,
        )
        subex = d.subexpressions[k]
        d.subexpression_linearity[k] = subex.linearity
        if d.want_hess
            empty!(coloring_storage)
            compute_gradient_sparsity!(coloring_storage, subex.nd)
            # union with all dependent expressions
            for idx in list_subexpressions(subex.nd)
                union!(coloring_storage, subexpression_variables[idx])
            end
            subexpression_variables[k] = collect(coloring_storage)
            empty!(coloring_storage)
            linearity = classify_linearity(
                subex.nd,
                subex.adj,
                d.subexpression_linearity,
            )
            edgelist = compute_hessian_sparsity(
                subex.nd,
                subex.adj,
                linearity,
                coloring_storage,
                subexpression_edgelist,
                subexpression_variables,
            )
            subexpression_edgelist[k] = edgelist
        end
    end

    if :ExprGraph in requested_features
        d.subexpressions_as_julia_expressions =
            Array{Any}(undef, length(subexpr))
        for k in d.subexpression_order
            ex = d.subexpressions[k]
            d.subexpressions_as_julia_expressions[k] = _tape_to_expr(
                d.model,
                1,
                nldata.nlexpr[k].nd,
                ex.adj,
                ex.const_values,
                d.parameter_values,
                d.subexpressions_as_julia_expressions,
                nldata.user_operators,
                true,
                true,
            )
        end
    end

    max_chunk = 1

    if d.has_nlobj
        nd = main_expressions[1]
        d.objective = _FunctionStorage(
            nd,
            nldata.nlobj.const_values,
            num_variables_,
            coloring_storage,
            d.want_hess,
            d.subexpressions,
            individual_order[1],
            d.subexpression_linearity,
            subexpression_edgelist,
            subexpression_variables,
            moi_index_to_consecutive_index,
        )
        max_expr_length = max(max_expr_length, length(d.objective.nd))
        max_chunk = max(max_chunk, size(d.objective.seed_matrix, 2))
    end

    for k in 1:length(nldata.nlconstr)
        nlconstr = nldata.nlconstr[k]
        idx = (d.has_nlobj) ? k + 1 : k
        nd = main_expressions[idx]
        push!(
            d.constraints,
            _FunctionStorage(
                nd,
                nlconstr.terms.const_values,
                num_variables_,
                coloring_storage,
                d.want_hess,
                d.subexpressions,
                individual_order[idx],
                d.subexpression_linearity,
                subexpression_edgelist,
                subexpression_variables,
                moi_index_to_consecutive_index,
            ),
        )
        max_expr_length = max(max_expr_length, length(d.constraints[end].nd))
        max_chunk = max(max_chunk, size(d.constraints[end].seed_matrix, 2))
    end

    max_chunk = min(max_chunk, 10) # 10 is hardcoded upper bound to avoid excess memory allocation

    if d.want_hess || want_hess_storage # storage for Hess or HessVec
        d.input_ϵ = Array{Float64}(undef, max_chunk * num_variables_)
        d.output_ϵ = Array{Float64}(undef, max_chunk * num_variables_)
        d.forward_storage_ϵ = Array{Float64}(undef, max_chunk * max_expr_length)
        d.partials_storage_ϵ =
            Array{Float64}(undef, max_chunk * max_expr_length)
        d.reverse_storage_ϵ = Array{Float64}(undef, max_chunk * max_expr_length)
        d.subexpression_forward_values_ϵ =
            Array{Float64}(undef, max_chunk * length(d.subexpressions))
        d.subexpression_reverse_values_ϵ =
            Array{Float64}(undef, max_chunk * length(d.subexpressions))
        for k in d.subexpression_order
            subex = d.subexpressions[k]
            subex.forward_storage_ϵ =
                zeros(Float64, max_chunk * length(subex.nd))
            subex.partials_storage_ϵ =
                zeros(Float64, max_chunk * length(subex.nd))
            subex.reverse_storage_ϵ =
                zeros(Float64, max_chunk * length(subex.nd))
        end
        d.max_chunk = max_chunk
        if d.want_hess
            d.hessian_sparsity = _hessian_lagrangian_structure(d)
            # JIT warm-up
            # TODO: rewrite without MPB
            #MathProgBase.eval_hessian_lagrangian(d, Array{Float64}(undef,length(d.hess_I)), d.model.colVal, 1.0, ones(MathProgBase.numconstr(d.model)))
        end
    end

    # JIT warm-up
    # TODO: rewrite without MPB
    # if :Grad in requested_features
    #     MOI.eval_objective_gradient(d, zeros(numVar), d.model.colVal)
    #     MOI.eval_constraint(d, zeros(MathProgBase.numconstr(d.model)), d.model.colVal)
    # end

    # reset timers
    d.eval_objective_timer = 0
    d.eval_objective_gradient_timer = 0
    d.eval_constraint_timer = 0
    d.eval_constraint_jacobian_timer = 0
    d.eval_hessian_lagrangian_timer = 0

    return nothing
end

function _recompute_disable_2ndorder(evaluator::NLPEvaluator)
    # Check if we have any user-defined operators, in which case we need to
    # disable hessians. The result of features_available depends on this.
    nldata::_NLPData = evaluator.model.nlp_data
    has_nlobj = nldata.nlobj !== nothing
    has_user_mv_operator = false
    for nlexpr in nldata.nlexpr
        has_user_mv_operator |=
            _Derivatives.has_user_multivariate_operators(nlexpr.nd)
    end
    if has_nlobj
        has_user_mv_operator |=
            _Derivatives.has_user_multivariate_operators(nldata.nlobj.nd)
    end
    for nlconstr in nldata.nlconstr
        has_user_mv_operator |=
            _Derivatives.has_user_multivariate_operators(nlconstr.terms.nd)
    end
    evaluator.disable_2ndorder = has_user_mv_operator
    return
end

function MOI.features_available(d::NLPEvaluator)
    _recompute_disable_2ndorder(d)
    features = [:Grad, :Jac, :JacVec, :ExprGraph]
    if !d.disable_2ndorder
        push!(features, :Hess)
        push!(features, :HessVec)
    end
    return features
end

function _forward_eval_all(d::NLPEvaluator, x)
    # do a forward pass on all expressions at x
    subexpr_values = d.subexpression_forward_values
    user_operators =
        d.model.nlp_data.user_operators::_Derivatives.UserOperatorRegistry
    user_input_buffer = d.jac_storage
    user_output_buffer = d.user_output_buffer
    for k in d.subexpression_order
        ex = d.subexpressions[k]
        subexpr_values[k] = forward_eval(
            ex.forward_storage,
            ex.partials_storage,
            ex.nd,
            ex.adj,
            ex.const_values,
            d.parameter_values,
            x,
            subexpr_values,
            user_input_buffer,
            user_output_buffer,
            user_operators,
        )
    end
    if d.has_nlobj
        obj = d.objective
        forward_eval(
            obj.forward_storage,
            obj.partials_storage,
            obj.nd,
            obj.adj,
            obj.const_values,
            d.parameter_values,
            x,
            subexpr_values,
            user_input_buffer,
            user_output_buffer,
            user_operators,
        )
    end
    for con in d.constraints
        forward_eval(
            con.forward_storage,
            con.partials_storage,
            con.nd,
            con.adj,
            con.const_values,
            d.parameter_values,
            x,
            subexpr_values,
            user_input_buffer,
            user_output_buffer,
            user_operators,
        )
    end
end

function _reverse_eval_all(d::NLPEvaluator, x)
    # do a reverse pass on all expressions at x
    subexpr_reverse_values = d.subexpression_reverse_values
    subexpr_values = d.subexpression_forward_values
    grad_storage = d.jac_storage
    for k in d.subexpression_order
        ex = d.subexpressions[k]
        reverse_eval(ex.reverse_storage, ex.partials_storage, ex.nd, ex.adj)
    end
    if d.has_nlobj
        obj = d.objective
        reverse_eval(obj.reverse_storage, obj.partials_storage, obj.nd, obj.adj)
    end
    for con in d.constraints
        reverse_eval(con.reverse_storage, con.partials_storage, con.nd, con.adj)
    end
    return copyto!(d.last_x, x)
end

function MOI.eval_objective(d::NLPEvaluator, x)
    d.eval_objective_timer += @elapsed begin
        if d.last_x != x
            _forward_eval_all(d, x)
            _reverse_eval_all(d, x)
        end
        val = zero(eltype(x))
        if d.has_nlobj
            val = d.objective.forward_storage[1]
        else
            error("No nonlinear objective.")
        end
    end
    return val
end

function MOI.eval_objective_gradient(d::NLPEvaluator, g, x)
    d.eval_objective_gradient_timer += @elapsed begin
        if d.last_x != x
            _forward_eval_all(d, x)
            _reverse_eval_all(d, x)
        end
        if d.has_nlobj
            fill!(g, 0.0)
            ex = d.objective
            subexpr_reverse_values = d.subexpression_reverse_values
            subexpr_reverse_values[ex.dependent_subexpressions] .= 0.0
            reverse_extract(
                g,
                ex.reverse_storage,
                ex.nd,
                ex.adj,
                subexpr_reverse_values,
                1.0,
            )
            for i in length(ex.dependent_subexpressions):-1:1
                k = ex.dependent_subexpressions[i]
                subexpr = d.subexpressions[k]
                reverse_extract(
                    g,
                    subexpr.reverse_storage,
                    subexpr.nd,
                    subexpr.adj,
                    subexpr_reverse_values,
                    subexpr_reverse_values[k],
                )
            end
        else
            error("No nonlinear objective.")
        end
    end
    return
end

function MOI.eval_constraint(d::NLPEvaluator, g, x)
    d.eval_constraint_timer += @elapsed begin
        if d.last_x != x
            _forward_eval_all(d, x)
            _reverse_eval_all(d, x)
        end

        for i in 1:length(d.constraints)
            g[i] = d.constraints[i].forward_storage[1]
        end
    end
    return
end

function MOI.eval_constraint_jacobian(d::NLPEvaluator, J, x)
    d.eval_constraint_jacobian_timer += @elapsed begin
        if d.last_x != x
            _forward_eval_all(d, x)
            _reverse_eval_all(d, x)
        end
        fill!(J, 0.0)
        grad_storage = d.jac_storage
        subexpr_reverse_values = d.subexpression_reverse_values
        idx = 0
        for ex in d.constraints
            nzidx = ex.grad_sparsity
            for i in nzidx
                @inbounds grad_storage[i] = 0.0
            end
            for i in ex.dependent_subexpressions
                @inbounds subexpr_reverse_values[i] = 0.0
            end

            reverse_extract(
                grad_storage,
                ex.reverse_storage,
                ex.nd,
                ex.adj,
                subexpr_reverse_values,
                1.0,
            )
            for i in length(ex.dependent_subexpressions):-1:1
                k = ex.dependent_subexpressions[i]
                subexpr = d.subexpressions[k]
                reverse_extract(
                    grad_storage,
                    subexpr.reverse_storage,
                    subexpr.nd,
                    subexpr.adj,
                    subexpr_reverse_values,
                    subexpr_reverse_values[k],
                )
            end

            for k in 1:length(nzidx)
                J[idx+k] = grad_storage[nzidx[k]]
            end
            idx += length(nzidx)
        end
    end
    return
end

function MOI.eval_constraint_jacobian_product(
    d::NLPEvaluator,
    y::AbstractVector{Float64},
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
)
    fill!(y, 0.0)
    jac_struct = MOI.jacobian_structure(d)
    nnz_jac = length(jac_struct)
    J = zeros(Float64, nnz_jac)
    MOI.eval_constraint_jacobian(d, J, x)
    for k in 1:nnz_jac
        i, j = jac_struct[k]
        y[i] += J[k] * w[j]
    end
    return y
end

function MOI.eval_constraint_jacobian_transpose_product(
    d::NLPEvaluator,
    y::AbstractVector{Float64},
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
)
    fill!(y, 0.0)
    jac_struct = MOI.jacobian_structure(d)
    nnz_jac = length(jac_struct)
    J = zeros(Float64, nnz_jac)
    MOI.eval_constraint_jacobian(d, J, x)
    for k in 1:nnz_jac
        i, j = jac_struct[k]
        y[j] += J[k] * w[i]
    end
    return y
end

function MOI.eval_hessian_lagrangian_product(
    d::NLPEvaluator,
    h::AbstractVector{Float64}, # output vector
    x::AbstractVector{Float64}, # current solution
    v::AbstractVector{Float64}, # rhs vector
    σ::Float64,                 # multiplier for objective
    μ::AbstractVector{Float64}, # multipliers for each constraint
)
    nldata = d.model.nlp_data::_NLPData

    if d.last_x != x
        _forward_eval_all(d, x)
        _reverse_eval_all(d, x)
    end

    fill!(h, 0.0)

    input_ϵ = reinterpret(ForwardDiff.Partials{1,Float64}, d.input_ϵ)
    output_ϵ = reinterpret(ForwardDiff.Partials{1,Float64}, d.output_ϵ)
    for i in 1:length(x)
        input_ϵ[i] = ForwardDiff.Partials((v[i],))
    end

    # forward evaluate all subexpressions once
    subexpr_forward_values_ϵ = reinterpret(
        ForwardDiff.Partials{1,Float64},
        d.subexpression_forward_values_ϵ,
    )
    subexpr_reverse_values_ϵ = reinterpret(
        ForwardDiff.Partials{1,Float64},
        d.subexpression_reverse_values_ϵ,
    )
    forward_storage_ϵ =
        reinterpret(ForwardDiff.Partials{1,Float64}, d.forward_storage_ϵ)
    reverse_storage_ϵ =
        reinterpret(ForwardDiff.Partials{1,Float64}, d.reverse_storage_ϵ)
    partials_storage_ϵ =
        reinterpret(ForwardDiff.Partials{1,Float64}, d.partials_storage_ϵ)
    for expridx in d.subexpression_order
        subexpr = d.subexpressions[expridx]
        sub_forward_storage_ϵ = reinterpret(
            ForwardDiff.Partials{1,Float64},
            subexpr.forward_storage_ϵ,
        )
        sub_partials_storage_ϵ = reinterpret(
            ForwardDiff.Partials{1,Float64},
            subexpr.partials_storage_ϵ,
        )
        subexpr_forward_values_ϵ[expridx] = forward_eval_ϵ(
            subexpr.forward_storage,
            sub_forward_storage_ϵ,
            subexpr.partials_storage,
            sub_partials_storage_ϵ,
            subexpr.nd,
            subexpr.adj,
            input_ϵ,
            subexpr_forward_values_ϵ,
            nldata.user_operators,
        )
    end
    # we only need to do one reverse pass through the subexpressions as well
    zero_ϵ = zero(ForwardDiff.Partials{1,Float64})
    fill!(subexpr_reverse_values_ϵ, zero_ϵ)
    fill!(d.subexpression_reverse_values, 0.0)
    fill!(reverse_storage_ϵ, zero_ϵ)
    fill!(output_ϵ, zero_ϵ)
    if d.has_nlobj
        ex = d.objective
        forward_eval_ϵ(
            ex.forward_storage,
            forward_storage_ϵ,
            ex.partials_storage,
            partials_storage_ϵ,
            ex.nd,
            ex.adj,
            input_ϵ,
            subexpr_forward_values_ϵ,
            nldata.user_operators,
        )
        reverse_eval_ϵ(
            output_ϵ,
            ex.reverse_storage,
            reverse_storage_ϵ,
            ex.partials_storage,
            partials_storage_ϵ,
            ex.nd,
            ex.adj,
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            σ,
            zero_ϵ,
        )
    end

    for i in 1:length(d.constraints)
        ex = d.constraints[i]
        l = μ[i]
        forward_eval_ϵ(
            ex.forward_storage,
            forward_storage_ϵ,
            ex.partials_storage,
            partials_storage_ϵ,
            ex.nd,
            ex.adj,
            input_ϵ,
            subexpr_forward_values_ϵ,
            nldata.user_operators,
        )
        reverse_eval_ϵ(
            output_ϵ,
            ex.reverse_storage,
            reverse_storage_ϵ,
            ex.partials_storage,
            partials_storage_ϵ,
            ex.nd,
            ex.adj,
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            l,
            zero_ϵ,
        )
    end

    for i in length(d.subexpression_order):-1:1
        expridx = d.subexpression_order[i]
        subexpr = d.subexpressions[expridx]
        sub_reverse_storage_ϵ = reinterpret(
            ForwardDiff.Partials{1,Float64},
            subexpr.reverse_storage_ϵ,
        )
        sub_partials_storage_ϵ = reinterpret(
            ForwardDiff.Partials{1,Float64},
            subexpr.partials_storage_ϵ,
        )
        reverse_eval_ϵ(
            output_ϵ,
            subexpr.reverse_storage,
            sub_reverse_storage_ϵ,
            subexpr.partials_storage,
            sub_partials_storage_ϵ,
            subexpr.nd,
            subexpr.adj,
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            d.subexpression_reverse_values[expridx],
            subexpr_reverse_values_ϵ[expridx],
        )
    end

    for i in 1:length(x)
        h[i] += output_ϵ[i].values[1]
    end
end

function MOI.eval_hessian_lagrangian(
    d::NLPEvaluator,
    H::AbstractVector{Float64},      # Sparse hessian entry vector
    x::AbstractVector{Float64},      # Current solution
    obj_factor::Float64,             # Lagrangian multiplier for objective
    lambda::AbstractVector{Float64}, # Multipliers for each constraint
)
    nldata = d.model.nlp_data::_NLPData

    d.want_hess || error(
        "Hessian computations were not requested on the call to initialize!.",
    )

    if d.last_x != x
        _forward_eval_all(d, x)
        _reverse_eval_all(d, x)
    end

    d.eval_hessian_lagrangian_timer += @elapsed begin
        fill!(d.input_ϵ, 0.0)
        recovery_tmp_storage = d.output_ϵ
        nzcount = 0

        if d.has_nlobj
            ex = d.objective
            chunk = min(size(ex.seed_matrix, 2), d.max_chunk)
            if chunk == 1
                # skip dynamic dispatch
                nzthis = _hessian_slice(
                    d,
                    ex,
                    x,
                    H,
                    obj_factor,
                    nzcount,
                    recovery_tmp_storage,
                    Val{1},
                )::Int
            else
                nzthis = _hessian_slice(
                    d,
                    ex,
                    x,
                    H,
                    obj_factor,
                    nzcount,
                    recovery_tmp_storage,
                    Val{chunk},
                )::Int
            end
            nzcount += nzthis
        end # else, obj_factor is ignored.

        for i in 1:length(d.constraints)
            ex = d.constraints[i]
            chunk = min(size(ex.seed_matrix, 2), d.max_chunk)
            if chunk == 1
                nzthis = _hessian_slice(
                    d,
                    ex,
                    x,
                    H,
                    lambda[i],
                    nzcount,
                    recovery_tmp_storage,
                    Val{1},
                )::Int
            else
                nzthis = _hessian_slice(
                    d,
                    ex,
                    x,
                    H,
                    lambda[i],
                    nzcount,
                    recovery_tmp_storage,
                    Val{chunk},
                )::Int
            end
            nzcount += nzthis
        end
    end

    return
end

function _hessian_slice_inner(
    d,
    ex,
    input_ϵ,
    output_ϵ,
    ::Type{Val{CHUNK}},
) where {CHUNK}
    subexpr_forward_values_ϵ = _reinterpret_unsafe(
        ForwardDiff.Partials{CHUNK,Float64},
        d.subexpression_forward_values_ϵ,
    )
    subexpr_reverse_values_ϵ = _reinterpret_unsafe(
        ForwardDiff.Partials{CHUNK,Float64},
        d.subexpression_reverse_values_ϵ,
    )
    forward_storage_ϵ = _reinterpret_unsafe(
        ForwardDiff.Partials{CHUNK,Float64},
        d.forward_storage_ϵ,
    )
    reverse_storage_ϵ = _reinterpret_unsafe(
        ForwardDiff.Partials{CHUNK,Float64},
        d.reverse_storage_ϵ,
    )
    partials_storage_ϵ = _reinterpret_unsafe(
        ForwardDiff.Partials{CHUNK,Float64},
        d.partials_storage_ϵ,
    )
    zero_ϵ = zero(ForwardDiff.Partials{CHUNK,Float64})

    user_operators =
        d.model.nlp_data.user_operators::_Derivatives.UserOperatorRegistry
    # do a forward pass
    for expridx in ex.dependent_subexpressions
        subexpr = d.subexpressions[expridx]
        sub_forward_storage_ϵ = _reinterpret_unsafe(
            ForwardDiff.Partials{CHUNK,Float64},
            subexpr.forward_storage_ϵ,
        )
        sub_partials_storage_ϵ = _reinterpret_unsafe(
            ForwardDiff.Partials{CHUNK,Float64},
            subexpr.partials_storage_ϵ,
        )
        subexpr_forward_values_ϵ[expridx] = forward_eval_ϵ(
            subexpr.forward_storage,
            sub_forward_storage_ϵ,
            subexpr.partials_storage,
            sub_partials_storage_ϵ,
            subexpr.nd,
            subexpr.adj,
            input_ϵ,
            subexpr_forward_values_ϵ,
            user_operators,
        )
    end
    forward_eval_ϵ(
        ex.forward_storage,
        forward_storage_ϵ,
        ex.partials_storage,
        partials_storage_ϵ,
        ex.nd,
        ex.adj,
        input_ϵ,
        subexpr_forward_values_ϵ,
        user_operators,
    )

    # do a reverse pass
    @inbounds for idx in ex.dependent_subexpressions
        subexpr_reverse_values_ϵ[idx] = zero_ϵ
        d.subexpression_reverse_values[idx] = 0.0
    end

    reverse_eval_ϵ(
        output_ϵ,
        ex.reverse_storage,
        reverse_storage_ϵ,
        ex.partials_storage,
        partials_storage_ϵ,
        ex.nd,
        ex.adj,
        d.subexpression_reverse_values,
        subexpr_reverse_values_ϵ,
        1.0,
        zero_ϵ,
    )
    for i in length(ex.dependent_subexpressions):-1:1
        expridx = ex.dependent_subexpressions[i]
        subexpr = d.subexpressions[expridx]
        sub_reverse_storage_ϵ = _reinterpret_unsafe(
            ForwardDiff.Partials{CHUNK,Float64},
            subexpr.reverse_storage_ϵ,
        )
        sub_partials_storage_ϵ = _reinterpret_unsafe(
            ForwardDiff.Partials{CHUNK,Float64},
            subexpr.partials_storage_ϵ,
        )
        reverse_eval_ϵ(
            output_ϵ,
            subexpr.reverse_storage,
            sub_reverse_storage_ϵ,
            subexpr.partials_storage,
            sub_partials_storage_ϵ,
            subexpr.nd,
            subexpr.adj,
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            d.subexpression_reverse_values[expridx],
            subexpr_reverse_values_ϵ[expridx],
        )
    end
end

function _hessian_slice(
    d,
    ex,
    x,
    H,
    scale,
    nzcount,
    recovery_tmp_storage,
    ::Type{Val{CHUNK}},
) where {CHUNK}
    nzthis = length(ex.hess_I)
    if ex.linearity == LINEAR
        @assert nzthis == 0
        return 0
    end
    R = ex.seed_matrix
    Coloring.prepare_seed_matrix!(R, ex.rinfo)
    local_to_global_idx = ex.rinfo.local_indices

    zero_ϵ = zero(ForwardDiff.Partials{CHUNK,Float64})

    input_ϵ_raw = d.input_ϵ
    output_ϵ_raw = d.output_ϵ
    input_ϵ =
        _reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64}, input_ϵ_raw)
    output_ϵ =
        _reinterpret_unsafe(ForwardDiff.Partials{CHUNK,Float64}, output_ϵ_raw)

    # compute hessian-vector products
    num_products = size(R, 2) # number of hessian-vector products
    num_chunks = div(num_products, CHUNK)
    @assert size(R, 1) == length(local_to_global_idx)

    for k in 1:CHUNK:CHUNK*num_chunks
        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up R[r,k,k+1,...,k+CHUNK-1] into input_ϵ
            for s in 1:CHUNK
                input_ϵ_raw[(idx-1)*CHUNK+s] = R[r, k+s-1]
            end
            @inbounds output_ϵ[idx] = zero_ϵ
        end

        _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, Val{CHUNK})

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into R[r,k,k+1,...,k+CHUNK-1]
            for s in 1:CHUNK
                R[r, k+s-1] = output_ϵ_raw[(idx-1)*CHUNK+s]
            end
            @inbounds input_ϵ[idx] = zero_ϵ
        end
    end

    # leftover chunk
    remaining = num_products - CHUNK * num_chunks
    if remaining > 0
        k = CHUNK * num_chunks + 1
        for r in 1:length(local_to_global_idx)
            # set up directional derivatives
            @inbounds idx = local_to_global_idx[r]
            # load up R[r,k,k+1,...,k+remaining-1] into input_ϵ
            for s in 1:remaining
                # leave junk in the unused components
                input_ϵ_raw[(idx-1)*CHUNK+s] = R[r, k+s-1]
            end
            @inbounds output_ϵ[idx] = zero_ϵ
        end

        _hessian_slice_inner(d, ex, input_ϵ, output_ϵ, Val{CHUNK})

        # collect directional derivatives
        for r in 1:length(local_to_global_idx)
            idx = local_to_global_idx[r]
            # load output_ϵ into R[r,k,k+1,...,k+remaining-1]
            for s in 1:remaining
                R[r, k+s-1] = output_ϵ_raw[(idx-1)*CHUNK+s]
            end
            @inbounds input_ϵ[idx] = zero_ϵ
        end
    end

    # Output is in R, now recover

    #output_slice = view(H, (nzcount+1):(nzcount+nzthis))
    output_slice = _VectorView(nzcount, nzthis, pointer(H))
    Coloring.recover_from_matmat!(
        output_slice,
        R,
        ex.rinfo,
        recovery_tmp_storage,
    )
    _rmul!(output_slice, scale)
    return nzthis
end

function MOI.jacobian_structure(d::NLPEvaluator)
    jacobian_sparsity = Tuple{Int64,Int64}[]
    for row in 1:length(d.constraints)
        row_sparsity = d.constraints[row].grad_sparsity
        for idx in row_sparsity
            push!(jacobian_sparsity, (row, idx))
        end
    end
    return jacobian_sparsity
end
function MOI.hessian_lagrangian_structure(d::NLPEvaluator)
    d.want_hess || error(
        "Hessian computations were not requested on the call to initialize!.",
    )
    return d.hessian_sparsity
end
function _hessian_lagrangian_structure(d::NLPEvaluator)
    hessian_sparsity = Tuple{Int64,Int64}[]
    if d.has_nlobj
        for idx in 1:length(d.objective.hess_I)
            push!(
                hessian_sparsity,
                (d.objective.hess_I[idx], d.objective.hess_J[idx]),
            )
        end
    end
    for ex in d.constraints
        for idx in 1:length(ex.hess_I)
            push!(hessian_sparsity, (ex.hess_I[idx], ex.hess_J[idx]))
        end
    end
    return hessian_sparsity
end

mutable struct _VariablePrintWrapper
    v::VariableRef
    mode::Any
end

function Base.show(io::IO, v::_VariablePrintWrapper)
    return print(io, function_string(v.mode, v.v))
end

mutable struct _ParameterPrintWrapper
    model::Model
    idx::Int
    mode::Any
end

function Base.show(io::IO, wrapper::_ParameterPrintWrapper)
    p = NonlinearParameter(wrapper.model, wrapper.idx)
    for (k, v) in object_dictionary(p.model)
        if v == p
            print(io, k)
            return
        end
    end
    # No named parameter; use a generic name.
    if wrapper.mode == MIME("text/latex")
        print(io, "parameter_{$(p.index)}")
    else
        print(io, "parameter[$(p.index)]")
    end
    return
end

mutable struct _SubexpressionPrintWrapper
    idx::Int
    mode::Any
end

function Base.show(io::IO, s::_SubexpressionPrintWrapper)
    if s.mode == MIME("text/latex")
        print(io, "subexpression_{$(s.idx)}")
    else
        print(io, "subexpression[$(s.idx)]")
    end
end

# we splat in the subexpressions (for now)
function _tape_to_expr(
    m::Model,
    k,
    nd::Vector{NodeData},
    adj,
    const_values,
    parameter_values,
    subexpressions::Vector{Any},
    user_operators::_Derivatives.UserOperatorRegistry,
    generic_variable_names::Bool,
    splat_subexpressions::Bool,
    print_mode = MIME("text/plain"),
)
    children_arr = rowvals(adj)

    nod = nd[k]
    if nod.nodetype == MOIVARIABLE
        if generic_variable_names
            return Expr(:ref, :x, _MOIVAR(nod.index))
        else
            # mode only matters when generic_variable_names == false
            return _VariablePrintWrapper(
                VariableRef(m, _MOIVAR(nod.index)),
                print_mode,
            )
        end
    elseif nod.nodetype == VALUE
        return const_values[nod.index]
    elseif nod.nodetype == SUBEXPRESSION
        if splat_subexpressions
            return subexpressions[nod.index]
        else
            return _SubexpressionPrintWrapper(nod.index, print_mode)
        end
    elseif nod.nodetype == PARAMETER
        if splat_subexpressions
            return parameter_values[nod.index]
        else
            return _ParameterPrintWrapper(m, nod.index, print_mode)
        end
    elseif nod.nodetype == CALL
        op = nod.index
        opsymbol = :error
        if op < _Derivatives.USER_OPERATOR_ID_START
            opsymbol = operators[op]
        else
            for (key, value) in user_operators.multivariate_operator_to_id
                if value == op - _Derivatives.USER_OPERATOR_ID_START + 1
                    opsymbol = key
                end
            end
        end
        @assert opsymbol != :error
        children_idx = nzrange(adj, k)
        if opsymbol == :+ && length(children_idx) == 0
            return 0
        elseif opsymbol == :* && length(children_idx) == 0
            return 1
        end
        ex = Expr(:call, opsymbol)
        for cidx in children_idx
            push!(
                ex.args,
                _tape_to_expr(
                    m,
                    children_arr[cidx],
                    nd,
                    adj,
                    const_values,
                    parameter_values,
                    subexpressions,
                    user_operators,
                    generic_variable_names,
                    splat_subexpressions,
                    print_mode,
                ),
            )
        end
        return ex
    elseif nod.nodetype == CALLUNIVAR
        op = nod.index
        opsymbol = :error
        if op < _Derivatives.USER_UNIVAR_OPERATOR_ID_START
            opsymbol = univariate_operators[op]
        else
            for (key, value) in user_operators.univariate_operator_to_id
                if value == op - _Derivatives.USER_UNIVAR_OPERATOR_ID_START + 1
                    opsymbol = key
                end
            end
        end
        @assert opsymbol != :error
        cidx = first(nzrange(adj, k))
        return Expr(
            :call,
            opsymbol,
            _tape_to_expr(
                m,
                children_arr[cidx],
                nd,
                adj,
                const_values,
                parameter_values,
                subexpressions,
                user_operators,
                generic_variable_names,
                splat_subexpressions,
                print_mode,
            ),
        )
    elseif nod.nodetype == COMPARISON
        op = nod.index
        opsymbol = comparison_operators[op]
        children_idx = nzrange(adj, k)
        if length(children_idx) > 2
            ex = Expr(:comparison)
            for cidx in children_idx
                push!(
                    ex.args,
                    _tape_to_expr(
                        m,
                        children_arr[cidx],
                        nd,
                        adj,
                        const_values,
                        parameter_values,
                        subexpressions,
                        user_operators,
                        generic_variable_names,
                        splat_subexpressions,
                        print_mode,
                    ),
                )
                push!(ex.args, opsymbol)
            end
            pop!(ex.args)
        else
            ex = Expr(:call, opsymbol)
            push!(
                ex.args,
                _tape_to_expr(
                    m,
                    children_arr[children_idx[1]],
                    nd,
                    adj,
                    const_values,
                    parameter_values,
                    subexpressions,
                    user_operators,
                    generic_variable_names,
                    splat_subexpressions,
                    print_mode,
                ),
            )
            push!(
                ex.args,
                _tape_to_expr(
                    m,
                    children_arr[children_idx[2]],
                    nd,
                    adj,
                    const_values,
                    parameter_values,
                    subexpressions,
                    user_operators,
                    generic_variable_names,
                    splat_subexpressions,
                    print_mode,
                ),
            )
        end
        return ex
    elseif nod.nodetype == LOGIC
        op = nod.index
        opsymbol = logic_operators[op]
        children_idx = nzrange(adj, k)
        lhs = _tape_to_expr(
            m,
            children_arr[first(children_idx)],
            nd,
            adj,
            const_values,
            parameter_values,
            subexpressions,
            user_operators,
            generic_variable_names,
            splat_subexpressions,
            print_mode,
        )
        rhs = _tape_to_expr(
            m,
            children_arr[last(children_idx)],
            nd,
            adj,
            const_values,
            parameter_values,
            subexpressions,
            user_operators,
            generic_variable_names,
            splat_subexpressions,
            print_mode,
        )
        return Expr(opsymbol, lhs, rhs)
    end
    return error()
end

function MOI.objective_expr(d::NLPEvaluator)
    if d.has_nlobj
        ex = d.objective
        return _tape_to_expr(
            d.model,
            1,
            d.model.nlp_data.nlobj.nd,
            ex.adj,
            ex.const_values,
            d.parameter_values,
            d.subexpressions_as_julia_expressions,
            d.model.nlp_data.user_operators,
            true,
            true,
        )
    else
        error("No nonlinear objective present")
    end
end

function MOI.constraint_expr(d::NLPEvaluator, i::Integer)
    ex = d.constraints[i]
    constr = d.model.nlp_data.nlconstr[i]
    julia_expr = _tape_to_expr(
        d.model,
        1,
        constr.terms.nd,
        ex.adj,
        ex.const_values,
        d.parameter_values,
        d.subexpressions_as_julia_expressions,
        d.model.nlp_data.user_operators,
        true,
        true,
    )
    if _sense(constr) == :range
        return Expr(:comparison, constr.lb, :(<=), julia_expr, :(<=), constr.ub)
    else
        return Expr(:call, _sense(constr), julia_expr, _rhs(constr))
    end
end

"""
    _VarValueMap{T,F}

A lazy cache used for computing the primal variable solution in `forward_eval`.

This avoids the need to rewrite the nonlinear expressions from MOIVARIABLE to
VARIABLE, as well as eagerly computing the `var_value` for every variable. We
use a `cache` so we don't have to recompute variables we have already seen.
"""
struct _VarValueMap{T,F}
    model::Model
    var_value::F
    cache::Dict{Int64,T}
end
function Base.getindex(m::_VarValueMap{T}, moi_index::Int64) where {T}
    return get!(m.cache, moi_index) do
        return m.var_value(
            VariableRef(m.model, MOI.VariableIndex(moi_index)),
        )::T
    end
end

"""
    value(var_value::Function, ex::NonlinearExpression)

Evaluate `ex` using `var_value(v)` as the value for each variable `v`.
"""
function value(var_value::Function, ex::NonlinearExpression)
    model = ex.model
    nlp_data::_NLPData = model.nlp_data
    variable_values = _VarValueMap(model, var_value, Dict{Int64,Float64}())
    subexpressions = Vector{NodeData}[nl_expr.nd for nl_expr in nlp_data.nlexpr]
    original_ex = nlp_data.nlexpr[ex.index]
    subexpression_order = order_subexpressions(original_ex.nd, subexpressions)
    max_len = length(original_ex.nd)
    for k in subexpression_order
        max_len = max(max_len, length(subexpressions[k]))
    end
    subexpr_values = Vector{Float64}(undef, length(subexpressions))
    forward_storage = Vector{Float64}(undef, max_len)
    partials_storage = Vector{Float64}(undef, max_len)
    user_input_buffer = zeros(nlp_data.largest_user_input_dimension)
    user_output_buffer = zeros(nlp_data.largest_user_input_dimension)
    for k in subexpression_order # Compute value of dependent subexpressions.
        subexpr_nd = subexpressions[k]
        original_subexpr = nlp_data.nlexpr[k]
        subexpr_values[k] = forward_eval(
            forward_storage,
            partials_storage,
            subexpr_nd,
            adjmat(subexpr_nd),
            original_subexpr.const_values,
            nlp_data.nlparamvalues,
            variable_values,
            subexpr_values,
            user_input_buffer,
            user_output_buffer,
            nlp_data.user_operators,
        )
    end
    return forward_eval(
        forward_storage,
        partials_storage,
        original_ex.nd,
        adjmat(original_ex.nd),
        original_ex.const_values,
        nlp_data.nlparamvalues,
        variable_values,
        subexpr_values,
        user_input_buffer,
        user_output_buffer,
        nlp_data.user_operators,
    )
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

mutable struct _UserFunctionEvaluator <: MOI.AbstractNLPEvaluator
    f::Any
    ∇f::Any
    len::Int
end

# TODO: This is a slightly confusing use for AbstractNLPEvaluator. Switch to a
# better data structure.
function MOI.eval_objective(d::_UserFunctionEvaluator, x)
    @assert length(x) == d.len
    return d.f(x)::eltype(x)
end
function MOI.eval_objective_gradient(d::_UserFunctionEvaluator, grad, x)
    d.∇f(grad, x)
    return nothing
end

const _FORWARD_DIFF_METHOD_ERROR_HELPER = raw"""
Common reasons for this include:

 * the function assumes `Float64` will be passed as input, it must work for any
   generic `Real` type.
 * the function allocates temporary storage using `zeros(3)` or similar. This
   defaults to `Float64`, so use `zeros(T, 3)` instead.

As an example, instead of:
```julia
function my_function(x::Float64...)
    y = zeros(length(x))
    for i in 1:length(x)
        y[i] = x[i]^2
    end
    return sum(y)
end
```
use:
```julia
function my_function(x::T...) where {T<:Real}
    y = zeros(T, length(x))
    for i in 1:length(x)
        y[i] = x[i]^2
    end
    return sum(y)
end
```

Review the stacktrace below for more information, but it can often be hard to
understand why and where your function is failing. You can also debug this
outside of JuMP as follows:
```julia
import ForwardDiff

# If the input dimension is 1
x = 1.0
my_function(a) = a^2
ForwardDiff.derivative(my_function, x)

# If the input dimension is more than 1
x = [1.0, 2.0]
my_function(a, b) = a^2 + b^2
ForwardDiff.gradient(x -> my_function(x...), x)
```
"""

_intercept_ForwardDiff_MethodError(err, ::Any) = rethrow(err)

function _intercept_ForwardDiff_MethodError(::MethodError, s)
    return error(
        "JuMP's autodiff of the user-defined function $(s) failed with a " *
        "MethodError.\n\n$(_FORWARD_DIFF_METHOD_ERROR_HELPER)",
    )
end

function _UserFunctionEvaluator(
    dimension::Integer,
    f::Function,
    ::Type{T} = Float64;
    name = "",
) where {T}
    g = x -> f(x...)
    cfg = ForwardDiff.GradientConfig(g, zeros(T, dimension))
    ∇f = function (out, y)
        try
            ForwardDiff.gradient!(out, g, y, cfg)
        catch err
            _intercept_ForwardDiff_MethodError(err, name)
        end
    end
    return _UserFunctionEvaluator(g, ∇f, dimension)
end

"""
    _validate_register_assumptions(
        f::Function,
        name::Symbol,
        dimension::Integer,
    )

A function that attempts to check if `f` is suitable for registration via
[`register`](@ref) and throws an informative error if it is not.

Because we don't know the domain of `f`, this function may encounter false
negatives. But it should catch the majority of cases in which users supply
non-differentiable functions that rely on `::Float64` assumptions.
"""
function _validate_register_assumptions(
    f::Function,
    name::Symbol,
    dimension::Integer,
)
    # Assumption 1: check that `f` can be called with `Float64` arguments.
    y = 0.0
    try
        if dimension == 1
            y = f(0.0)
        else
            y = f(zeros(dimension))
        end
    catch
        # We hit some other error, perhaps we called a function like log(0).
        # Ignore for now, and hope that a useful error is shown to the user
        # during the solve.
    end
    if !(y isa Real)
        error(
            "Expected return type of `Float64` from the user-defined " *
            "function :$(name), but got `$(typeof(y))`.",
        )
    end
    # Assumption 2: check that `f` can be differentiated using `ForwardDiff`.
    try
        if dimension == 1
            ForwardDiff.derivative(f, 0.0)
        else
            ForwardDiff.gradient(x -> f(x...), zeros(dimension))
        end
    catch err
        if err isa MethodError
            error(
                "Unable to register the function :$name because it does not " *
                "support differentiation via ForwardDiff.\n\n" *
                _FORWARD_DIFF_METHOD_ERROR_HELPER,
            )
        end
        # We hit some other error, perhaps we called a function like log(0).
        # Ignore for now, and hope that a useful error is shown to the user
        # during the solve.
    end
    return
end

function _checked_derivative(f::F, s) where {F}
    return function (x)
        try
            return ForwardDiff.derivative(f, x)
        catch err
            _intercept_ForwardDiff_MethodError(err, s)
        end
    end
end

"""
    register(
        model::Model,
        s::Symbol,
        dimension::Integer,
        f::Function;
        autodiff:Bool = false,
    )

Register the user-defined function `f` that takes `dimension` arguments in
`model` as the symbol `s`.

The function `f` must support all subtypes of `Real` as arguments. Do not assume
that the inputs are `Float64`.

## Notes

 * For this method, you must explicitly set `autodiff = true`, because no
   user-provided gradient function `∇f` is given.
 * Second-derivative information is only computed if `dimension == 1`.
 * `s` does not have to be the same symbol as `f`, but it is generally more
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
    m::Model,
    s::Symbol,
    dimension::Integer,
    f::Function;
    autodiff::Bool = false,
)
    autodiff == true ||
        error("If only the function is provided, must set autodiff=true")
    _validate_register_assumptions(f, s, dimension)
    _init_NLP(m)
    if dimension == 1
        fprime = _checked_derivative(f, s)
        fprimeprime = x -> ForwardDiff.derivative(fprime, x)
        _Derivatives.register_univariate_operator!(
            m.nlp_data.user_operators,
            s,
            f,
            fprime,
            fprimeprime,
        )
    else
        m.nlp_data.largest_user_input_dimension =
            max(m.nlp_data.largest_user_input_dimension, dimension)
        _Derivatives.register_multivariate_operator!(
            m.nlp_data.user_operators,
            s,
            _UserFunctionEvaluator(dimension, f; name = s),
        )
    end
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
    m::Model,
    s::Symbol,
    dimension::Integer,
    f::Function,
    ∇f::Function;
    autodiff::Bool = false,
)
    _init_NLP(m)
    if dimension == 1
        autodiff == true || error(
            "Currently must provide 2nd order derivatives of univariate functions. Try setting autodiff=true.",
        )
        _validate_register_assumptions(∇f, s, dimension)
        fprimeprime = _checked_derivative(∇f, s)
        _Derivatives.register_univariate_operator!(
            m.nlp_data.user_operators,
            s,
            f,
            ∇f,
            fprimeprime,
        )
    else
        autodiff == false ||
            @warn("autodiff=true ignored since gradient is already provided.")
        m.nlp_data.largest_user_input_dimension =
            max(m.nlp_data.largest_user_input_dimension, dimension)
        d = _UserFunctionEvaluator(
            x -> f(x...),
            (g, x) -> ∇f(g, x...),
            dimension,
        )
        _Derivatives.register_multivariate_operator!(
            m.nlp_data.user_operators,
            s,
            d,
        )
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
    m::Model,
    s::Symbol,
    dimension::Integer,
    f::Function,
    ∇f::Function,
    ∇²f::Function,
)
    dimension == 1 || error(
        "Providing hessians for multivariate functions is not yet supported",
    )
    _init_NLP(m)
    return _Derivatives.register_univariate_operator!(
        m.nlp_data.user_operators,
        s,
        f,
        ∇f,
        ∇²f,
    )
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
    return NonlinearExpression(model, _NonlinearExprData(model, ex))
end

function add_NL_expression(::Model, ::Expr)
    return error(
        "`add_NL_expression` has been renamed to `add_nonlinear_expression`",
    )
end

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
    return set_objective(model, sense, _NonlinearExprData(model, x))
end

function set_NL_objective(::Model, ::MOI.OptimizationSense, ::Expr)
    return error(
        "`set_NL_objective` has been renamed to `set_nonlinear_objective`",
    )
end

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
    nl_constraints = model.nlp_data.nlconstr
    if isexpr(ex, :call) # One-sided constraint.
        # Simple comparison - move everything to the LHS
        op = ex.args[1]
        if op == :(==)
            lb = 0.0
            ub = 0.0
        elseif op == :(<=) || op == :(≤)
            lb = -Inf
            ub = 0.0
        elseif op == :(>=) || op == :(≥)
            lb = 0.0
            ub = Inf
        else
            error(
                "In expression ($ex): expected comparison operator (<=, >=," *
                " or ==).",
            )
        end
        lhs = :($(ex.args[2]) - $(ex.args[3]))
        c = _NonlinearConstraint(_NonlinearExprData(model, lhs), lb, ub)
        push!(nl_constraints, c)
        return ConstraintRef(
            model,
            NonlinearConstraintIndex(length(nl_constraints)),
            ScalarShape(),
        )
    elseif isexpr(ex, :comparison)
        # ranged row
        if (ex.args[2] != :<= && ex.args[2] != :≤) ||
           (ex.args[4] != :<= && ex.args[4] != :≤)
            error(
                "In expression ($ex): only ranged rows of the form lb <= " *
                "expr <= ub are supported.",
            )
        end
        lb = ex.args[1]
        ub = ex.args[5]
        if !isa(lb, Number) || !isa(ub, Number)
            error(
                "Interval constraint contains non-constant left- or " *
                "right-hand sides. Reformulate as two separate " *
                "constraints, or move all variables into the central term.",
            )
        end
        c = _NonlinearConstraint(_NonlinearExprData(model, ex.args[3]), lb, ub)
        push!(nl_constraints, c)
        return ConstraintRef(
            model,
            NonlinearConstraintIndex(length(nl_constraints)),
            ScalarShape(),
        )
    else
        # Unknown
        error(
            "In expression ($ex): constraints must be in one of the " *
            "following forms:\n" *
            "       expr1 <= expr2\n       expr1 >= expr2\n" *
            "       expr1 == expr2",
        )
    end
end

function add_NL_constraint(::Model, ::Expr)
    return error(
        "`add_NL_constraint` has been renamed to `add_nonlinear_constraint`",
    )
end
