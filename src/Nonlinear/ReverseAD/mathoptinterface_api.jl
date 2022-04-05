#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

function MOI.features_available(d::NLPEvaluator)
    # Check if we have any user-defined multivariate operators, in which case we
    # need to disable hessians. The result of features_available depends on this.
    d.disable_2ndorder =
        length(d.data.operators.registered_multivariate_operators) > 0
    if d.disable_2ndorder
        return [:Grad, :Jac, :JacVec]
    end
    return [:Grad, :Jac, :JacVec, :Hess, :HessVec]
end

function MOI.initialize(d::NLPEvaluator, requested_features::Vector{Symbol})
    # Check that we support the features requested by the user.
    available_features = MOI.features_available(d)
    for feature in requested_features
        if !(feature in available_features)
            error("Unsupported feature $feature")
        end
    end
    moi_index_to_consecutive_index =
        Dict(x => i for (i, x) in enumerate(d.ordered_variables))
    N = length(moi_index_to_consecutive_index)
    #
    largest_user_input_dimension = 1
    for op in d.data.operators.registered_multivariate_operators
        largest_user_input_dimension = max(largest_user_input_dimension, op.N)
    end
    d.objective = nothing
    d.user_output_buffer = zeros(largest_user_input_dimension)
    d.jac_storage = zeros(max(N, largest_user_input_dimension))
    d.constraints = _FunctionStorage[]
    d.last_x = fill(NaN, N)
    d.want_hess = :Hess in requested_features
    want_hess_storage = (:HessVec in requested_features) || d.want_hess
    coloring_storage = Coloring.IndexedSet(N)
    max_expr_length = 0
    #
    main_expressions = [c.expression.nodes for (_, c) in d.data.constraints]
    if d.data.objective !== nothing
        pushfirst!(main_expressions, d.data.objective.nodes)
    end
    d.subexpression_order, individual_order = _order_subexpressions(
        main_expressions,
        [expr.nodes for expr in d.data.expressions],
    )
    num_subexpressions = length(d.data.expressions)
    d.subexpression_linearity = Vector{Linearity}(undef, num_subexpressions)
    subexpression_variables = Vector{Vector{Int}}(undef, num_subexpressions)
    subexpression_edgelist =
        Vector{Set{Tuple{Int,Int}}}(undef, num_subexpressions)
    d.subexpressions = Vector{_SubexpressionStorage}(undef, num_subexpressions)
    d.subexpression_forward_values = zeros(num_subexpressions)
    d.subexpression_reverse_values = zeros(num_subexpressions)
    for k in d.subexpression_order
        # Only load expressions which actually are used
        d.subexpression_forward_values[k] = NaN
        subex = _SubexpressionStorage(
            d.data.expressions[k],
            d.subexpression_linearity,
            moi_index_to_consecutive_index,
        )
        d.subexpressions[k] = subex
        d.subexpression_linearity[k] = subex.linearity
        if d.want_hess
            empty!(coloring_storage)
            _compute_gradient_sparsity!(coloring_storage, subex.nodes)
            # union with all dependent expressions
            for idx in _list_subexpressions(subex.nodes)
                union!(coloring_storage, subexpression_variables[idx])
            end
            subexpression_variables[k] = collect(coloring_storage)
            empty!(coloring_storage)
            linearity = _classify_linearity(
                subex.nodes,
                subex.adj,
                d.subexpression_linearity,
            )
            edgelist = _compute_hessian_sparsity(
                subex.nodes,
                subex.adj,
                linearity,
                coloring_storage,
                subexpression_edgelist,
                subexpression_variables,
            )
            subexpression_edgelist[k] = edgelist
        end
    end
    max_chunk = 1
    if d.data.objective !== nothing
        d.objective = _FunctionStorage(
            main_expressions[1],
            d.data.objective.values,
            N,
            coloring_storage,
            d.want_hess,
            d.subexpressions,
            individual_order[1],
            d.subexpression_linearity,
            subexpression_edgelist,
            subexpression_variables,
            moi_index_to_consecutive_index,
        )
        max_expr_length = max(max_expr_length, length(d.objective.nodes))
        max_chunk = max(max_chunk, size(d.objective.seed_matrix, 2))
    end
    for (k, (_, constraint)) in enumerate(d.data.constraints)
        idx = d.data.objective !== nothing ? k + 1 : k
        push!(
            d.constraints,
            _FunctionStorage(
                main_expressions[idx],
                constraint.expression.values,
                N,
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
        max_expr_length = max(max_expr_length, length(d.constraints[end].nodes))
        max_chunk = max(max_chunk, size(d.constraints[end].seed_matrix, 2))
    end
    # 10 is hardcoded upper bound to avoid excess memory allocation
    max_chunk = min(max_chunk, 10)
    if d.want_hess || want_hess_storage
        d.input_ϵ = zeros(max_chunk * N)
        d.output_ϵ = zeros(max_chunk * N)
        #
        len = max_chunk * max_expr_length
        d.forward_storage_ϵ = zeros(len)
        d.partials_storage_ϵ = zeros(len)
        d.reverse_storage_ϵ = zeros(len)
        #
        len = max_chunk * length(d.subexpressions)
        d.subexpression_forward_values_ϵ = zeros(len)
        d.subexpression_reverse_values_ϵ = zeros(len)
        #
        for k in d.subexpression_order
            len = max_chunk * length(d.subexpressions[k].nodes)
            d.subexpressions[k].forward_storage_ϵ = zeros(len)
            d.subexpressions[k].partials_storage_ϵ = zeros(len)
            d.subexpressions[k].reverse_storage_ϵ = zeros(len)
        end
        d.max_chunk = max_chunk
        if d.want_hess
            _compute_hessian_lagrangian_structure(d)
        end
    end
    return
end

function MOI.eval_objective(d::NLPEvaluator, x)
    if d.objective === nothing
        error("No nonlinear objective.")
    end
    _reverse_mode(d, x)
    return d.objective.forward_storage[1]
end

function MOI.eval_objective_gradient(d::NLPEvaluator, g, x)
    if d.objective === nothing
        error("No nonlinear objective.")
    end
    _reverse_mode(d, x)
    fill!(g, 0.0)
    _extract_reverse_pass(g, d, d.objective)
    return
end

function MOI.eval_constraint(d::NLPEvaluator, g, x)
    _reverse_mode(d, x)
    for i in 1:length(d.constraints)
        g[i] = d.constraints[i].forward_storage[1]
    end
    return
end

function MOI.jacobian_structure(d::NLPEvaluator)
    J = Tuple{Int64,Int64}[]
    for (row, constraint) in enumerate(d.constraints)
        for col in constraint.grad_sparsity
            push!(J, (row, col))
        end
    end
    return J
end

function MOI.eval_constraint_jacobian(d::NLPEvaluator, J, x)
    _reverse_mode(d, x)
    fill!(J, 0.0)
    offset = 0
    for ex in d.constraints
        for i in ex.grad_sparsity
            d.jac_storage[i] = 0.0
        end
        _extract_reverse_pass(d.jac_storage, d, ex)
        for (k, idx) in enumerate(ex.grad_sparsity)
            J[offset+k] = d.jac_storage[idx]
        end
        offset += length(ex.grad_sparsity)
    end
    return
end

function MOI.eval_constraint_jacobian_product(d::NLPEvaluator, y, x, w)
    fill!(y, 0.0)
    J_struct = MOI.jacobian_structure(d)
    nnz = length(J_struct)
    J = zeros(nnz)
    MOI.eval_constraint_jacobian(d, J, x)
    for (k, (i, j)) in enumerate(J_struct)
        y[i] += J[k] * w[j]
    end
    return
end

function MOI.eval_constraint_jacobian_transpose_product(
    d::NLPEvaluator,
    y::AbstractVector{Float64},
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
)
    fill!(y, 0.0)
    J_struct = MOI.jacobian_structure(d)
    nnz = length(J_struct)
    J = zeros(nnz)
    MOI.eval_constraint_jacobian(d, J, x)
    for (k, (i, j)) in enumerate(J_struct)
        y[j] += J[k] * w[i]
    end
    return y
end

function MOI.hessian_lagrangian_structure(d::NLPEvaluator)
    @assert d.want_hess
    return d.hessian_sparsity
end

function _compute_hessian_lagrangian_structure(d::NLPEvaluator)
    d.hessian_sparsity = Tuple{Int64,Int64}[]
    if d.objective !== nothing
        append!(d.hessian_sparsity, zip(d.objective.hess_I, d.objective.hess_J))
    end
    for c in d.constraints
        append!(d.hessian_sparsity, zip(c.hess_I, c.hess_J))
    end
    return
end

function MOI.eval_hessian_lagrangian(d::NLPEvaluator, H, x, σ, μ)
    @assert d.want_hess
    _reverse_mode(d, x)
    fill!(d.input_ϵ, 0.0)
    offset = 0
    if d.objective !== nothing
        offset += _eval_hessian(d, d.objective, H, σ, offset)::Int
    end
    for (i, ex) in enumerate(d.constraints)
        offset += _eval_hessian(d, ex, H, μ[i], offset)::Int
    end
    return
end

function MOI.eval_hessian_lagrangian_product(d::NLPEvaluator, h, x, v, σ, μ)
    _reverse_mode(d, x)
    fill!(h, 0.0)
    T = ForwardDiff.Partials{1,Float64}
    input_ϵ = reinterpret(T, d.input_ϵ)
    output_ϵ = reinterpret(T, d.output_ϵ)
    for i in 1:length(x)
        input_ϵ[i] = ForwardDiff.Partials((v[i],))
    end
    # forward evaluate all subexpressions once
    subexpr_forward_values_ϵ = reinterpret(T, d.subexpression_forward_values_ϵ)
    for i in d.subexpression_order
        subexpr = d.subexpressions[i]
        subexpr_forward_values_ϵ[i] = _forward_eval_ϵ(
            subexpr,
            reinterpret(T, subexpr.forward_storage_ϵ),
            reinterpret(T, subexpr.partials_storage_ϵ),
            input_ϵ,
            subexpr_forward_values_ϵ,
            d.data.operators,
        )
    end
    # we only need to do one reverse pass through the subexpressions as well
    subexpr_reverse_values_ϵ = reinterpret(T, d.subexpression_reverse_values_ϵ)
    fill!(subexpr_reverse_values_ϵ, zero(T))
    fill!(d.subexpression_reverse_values, 0.0)
    fill!(d.reverse_storage_ϵ, 0.0)
    fill!(output_ϵ, zero(T))
    if d.objective !== nothing
        _forward_eval_ϵ(
            d.objective,
            reinterpret(T, d.forward_storage_ϵ),
            reinterpret(T, d.partials_storage_ϵ),
            input_ϵ,
            subexpr_forward_values_ϵ,
            d.data.operators,
        )
        _reverse_eval_ϵ(
            output_ϵ,
            d.objective,
            reinterpret(T, d.reverse_storage_ϵ),
            reinterpret(T, d.partials_storage_ϵ),
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            σ,
            zero(T),
        )
    end
    for (i, con) in enumerate(d.constraints)
        _forward_eval_ϵ(
            con,
            reinterpret(T, d.forward_storage_ϵ),
            reinterpret(T, d.partials_storage_ϵ),
            input_ϵ,
            subexpr_forward_values_ϵ,
            d.data.operators,
        )
        _reverse_eval_ϵ(
            output_ϵ,
            con,
            reinterpret(T, d.reverse_storage_ϵ),
            reinterpret(T, d.partials_storage_ϵ),
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            μ[i],
            zero(T),
        )
    end
    for i in length(d.subexpression_order):-1:1
        j = d.subexpression_order[i]
        subexpr = d.subexpressions[j]
        _reverse_eval_ϵ(
            output_ϵ,
            subexpr,
            reinterpret(T, subexpr.reverse_storage_ϵ),
            reinterpret(T, subexpr.partials_storage_ϵ),
            d.subexpression_reverse_values,
            subexpr_reverse_values_ϵ,
            d.subexpression_reverse_values[j],
            subexpr_reverse_values_ϵ[j],
        )
    end
    for i in 1:length(x)
        h[i] += output_ϵ[i].values[1]
    end
    return
end
