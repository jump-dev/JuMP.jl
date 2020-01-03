#  Copyright 2019, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modeling language for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

"""
    lp_rhs_perturbation_range(constraint::ConstraintRef;
                              feasibility_tolerance::Float64)
                              ::Tuple{Float64, Float64}

Gives the range by which the rhs coefficient can change and the current LP basis
remains feasible, i.e., where the shadow prices apply.

## Notes
- The rhs coefficient is the value right of the relation, i.e., b for the constraint when of the form a*x □ b, where □ is ≤, =, or ≥.
- The range denotes valid changes, e.g., for a*x <= b + Δ, the LP basis remains feasible for all Δ ∈ [l, u].
- `feasibility_tolerance` is the primal feasibility tolerance, this should
  preferably match the tolerance used by the solver. The default tolerance should
  however apply in most situations (c.f. "Computational Techniques of the
  Simplex Method" by István Maros, section 9.3.4).
"""
function lp_rhs_perturbation_range(constraint::ConstraintRef{Model, <:_MOICON}; feasibility_tolerance::Float64 = 1e-8)
    error("The perturbation range of rhs is not defined or not implemented for this type " *
          "of constraint.")
end

"""
Builds the standard form constraint matrix, variable bounds, and a vector of constraint reference, one for each row.
If Ã is the current constraint matrix, and the corresponding feasible set is of the form
s.t. Row_L <= Ãx <= Row_U,
     Var_L <=  x <= Var_U.
(Equality constraint is interpreted as, e.g., Row_L_i == Row_U_i,
and single sided constraints with an infinite value in the free directions)
Then the function returns a tuple with
    A = [Ã -I]
    var_lower = [Var_L, Row_L]
    var_upper = [Var_U, Row_U]
    affine_constraints::Vector{ConstraintRef} - references corresponding to Row_L <= Ãx <= Row_U
This represents the LP problem feasible set of the form
s.t.            Ax == 0
  var_lower <=   x <= var_upper
"""
function _std_matrix(model::Model)
    con_types = list_of_constraint_types(model)
    con_func_type = first.(con_types)
    if !all(broadcast(<:, AffExpr, con_func_type) .| broadcast(<:, VariableRef, con_func_type))
        error("The requested operation is not supported because the problem is not a linear optimization problem.")
    end
    con_set_type = last.(con_types)
    if !all(broadcast(<:, con_set_type, Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo, MOI.Interval}))
        error("The requested operation is not supported because the problem is not a linear optimization problem.")
    end
    vars = all_variables(model)
    var_to_idx = Dict{VariableRef, Int}()
    for (j, var) in enumerate(vars)
        var_to_idx[var] = j
    end
    col_j = Int[]
    row_i = Int[]
    coeffs = Float64[]

    var_lower = fill(-Inf, length(vars))
    var_upper = fill(Inf, length(vars))
    affine_constraints = ConstraintRef[]

    cur_row_idx = 0

    for (F, S) in con_types
        constraints = all_constraints(model, F, S)
        if F <: VariableRef
            for constraint in constraints
                con_obj = constraint_object(constraint)
                j = var_to_idx[con_obj.func]
                range = MOI.Interval(con_obj.set)
                var_lower[j] = max(var_lower[j], range.lower)
                var_upper[j] = min(var_upper[j], range.upper)
            end
            #issue 1892 makes it hard to know how to handle variable constraints, atm. assume unique variable constraints.
        else
            for constraint in constraints
                cur_row_idx += 1
                con_obj = constraint_object(constraint)
                con_terms = con_obj.func.terms
                push!(affine_constraints, constraint)
                for (var, coeff) in con_terms
                    push!(col_j, var_to_idx[var])
                    push!(row_i, cur_row_idx)
                    push!(coeffs, coeff)
                end
                push!(col_j, cur_row_idx + length(vars))
                push!(row_i, cur_row_idx)
                push!(coeffs, -1.0)

                range = MOI.Interval(con_obj.set)
                push!(var_lower, range.lower)
                push!(var_upper, range.upper)
            end
        end
    end

    A = sparse(row_i, col_j, coeffs, cur_row_idx, cur_row_idx + length(vars))
    return A, var_lower, var_upper, affine_constraints
end

"""
Returns the optimal primal values for the problem in standard form, c.f. `_std_matrix(model)`.
"""
function _std_primal_solution(model::Model, constraints::Vector{ConstraintRef})
    vars = all_variables(model)
    return [value.(vars); value.(constraints)]
end

"""
Returns the basis status of all variables of the problem in standard form, c.f. `_std_matrix(model)`.
"""
function _std_basis_status(model::Model)::Vector{MOI.BasisStatusCode}
    con_types = list_of_constraint_types(model)
    vars = all_variables(model)
    var_to_idx = Dict{VariableRef, Int}()
    for (j, var) in enumerate(vars)
        var_to_idx[var] = j
    end
    basis_status = fill(MOI.BASIC, length(vars))

    for (F, S) in con_types
        constraints = all_constraints(model, F, S)
        if F <: VariableRef
            for constraint in constraints
                con_obj = constraint_object(constraint)
                j = var_to_idx[con_obj.func]
                try
                    basis_status[j] = MOI.get(model, MOI.ConstraintBasisStatus(), constraint)
                catch e
                    error("The perturbation range of rhs is not available because the solver doesn't "*
                          "support ´ConstraintBasisStatus´.")
                end
                if S <: MOI.LessThan && basis_status[j] == MOI.NONBASIC
                    basis_status[j] = MOI.NONBASIC_AT_UPPER
                elseif S <: MOI.GreaterThan && basis_status[j] == MOI.NONBASIC
                    basis_status[j] = MOI.NONBASIC_AT_LOWER
                end
                # ´_std_matrix´ checks that no invalid sets are used.
            end
            #issue 1892 makes it hard to know how to handle variable constraints, atm. assume unique variable constraints.
        else
            for constraint in constraints
                con_basis_status = MOI.get(model, MOI.ConstraintBasisStatus(), constraint)
                if S <: MOI.LessThan && con_basis_status == MOI.NONBASIC
                    con_basis_status = MOI.NONBASIC_AT_UPPER
                elseif S <: MOI.GreaterThan && con_basis_status == MOI.NONBASIC
                    con_basis_status = MOI.NONBASIC_AT_LOWER
                end
                push!(basis_status, con_basis_status)
            end
        end
    end
    return basis_status
end

function lp_rhs_perturbation_range(constraint::ConstraintRef{Model, _MOICON{F, S}}; feasibility_tolerance::Float64 = 1e-8
                      )::Tuple{Float64, Float64} where {S <: Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo}, T, F <: Union{MOI.ScalarAffineFunction{T}, MOI.SingleVariable}}
    model = owner_model(constraint)
    if termination_status(model) != MOI.OPTIMAL
        error("The perturbation range of rhs is not available because the current solution "*
              "is not optimal.")
    end

    try
        con_status = MOI.get(model, MOI.ConstraintBasisStatus(), constraint)
        # The constraint is inactive.
        if con_status == MOI.BASIC
            if S <: MOI.LessThan
                return value(constraint) - constraint_object(constraint).set.upper, Inf
            end
            if S <: MOI.GreaterThan
                return -Inf, value(constraint) - constraint_object(constraint).set.lower
            end
            return 0.0, 0.0
        end
    catch e
        error("The perturbation range of rhs is not available because the solver doesn't "*
              "support ´ConstraintBasisStatus´.")
    end
    # The constraint is active.
    # TODO: This could be made more efficient because actually we only need the basic columns.
    A, L, U, constraints = _std_matrix(model)
    var_basis_status = _std_basis_status(model)
    x = _std_primal_solution(model, constraints)
    # Ax = 0 => B*x_B = -N x_N =: b_N
    # Perturb row i by Δ => ̃x_B = B^-1 (b_N + Δ e_i) = x_B + Δ ρ
    # Perturb var bound j by Δ => ̃x_B = B^-1 (b_N - Δ N_j) = x_B + Δ ρ
    # Find min and max values of Δ : L_B <= ̃x_B <= U_B
    # L_B <= x_B + Δ ρ <= U_B  -> L_B - x_B <= Δ ρ <= U_B - x_B
    # Δ <= (U_B - x_B)_+ ./ ρ_+,  (L_B - x_B)_- ./ p_-
    # Δ >= (U_B - x_B)_- ./ ρ_-,  (L_B - x_B)_+ ./ p_+
    basic = var_basis_status .== Ref(MOI.BASIC)
    B = A[:, basic]
    x_B = x[basic]
    local rho::Vector{Float64}
    if F <: MOI.SingleVariable
        var = constraint_object(constraint).func
        j = findfirst(isequal(var), all_variables(model))
        rho = -B \ Vector(A[:, j])
    else
        i = findfirst(isequal(constraint), constraints)
        ei = zeros(size(A, 1))
        ei[i] = 1.0
        rho = B \ ei
    end
    UmX_B = U[basic] - x_B
    LmX_B = L[basic] - x_B
    pos = rho .> 0.0
    neg = rho .< 0.0
    # Find the first basic variable bound that is strictly violated.
    lower_bounds_delta = [-Inf; UmX_B[neg]  ./ rho[neg]; LmX_B[pos] ./ rho[pos]]
    lower_bounds_delta_strict = lower_bounds_delta + [0.0; feasibility_tolerance ./ rho[neg]; -feasibility_tolerance ./ rho[pos]]
    upper_bounds_delta = [Inf; UmX_B[pos] ./ rho[pos]; LmX_B[neg] ./ rho[neg]]
    upper_bounds_delta_strict = upper_bounds_delta + [0.0; feasibility_tolerance ./ rho[pos]; -feasibility_tolerance ./ rho[neg]]

    lower_max_idx = last(findmax(lower_bounds_delta_strict))
    upper_min_idx = last(findmin(upper_bounds_delta_strict))
    # Return the unperturbed values
    return lower_bounds_delta[lower_max_idx], upper_bounds_delta[upper_min_idx]
end

"""
Returns the optimal reduced costs for the variables of the problem in standard form, c.f. `_std_matrix(model)`
"""
function _std_reduced_costs(model::Model, constraints::Vector{ConstraintRef})
    con_types = list_of_constraint_types(model)
    vars = all_variables(model)
    var_to_idx = Dict{VariableRef, Int}()
    for (j, var) in enumerate(vars)
        var_to_idx[var] = j
    end
    y = [zeros(length(vars)); dual.(constraints)]

    for (F, S) in con_types
        if F <: VariableRef
            constraints = all_constraints(model, F, S)
            for constraint in constraints
                con_obj = constraint_object(constraint)
                j = var_to_idx[con_obj.func]
                y[j] += dual(constraint)
            end
            #issue 1892 makes is hard to know how to handle variable constraints, atm. assume unique variable constraints
        end
    end

    if objective_sense(model) == MOI.MAX_SENSE
        y .*= -1.0
    end
    return y
end

"""
    lp_objective_perturbation_range(var::VariableRef;
                                    optimality_tolerance::Float64)
                                    ::Tuple{Float64, Float64}

Gives the range by which the cost coefficient can change and the current LP basis
remains optimal, i.e., the reduced costs remain valid.

## Notes
- The range denotes valid changes, Δ ∈ [l, u], for which cost[var] += Δ do not
  violate the current optimality conditions.
- `optimality_tolerance` is the dual feasibility tolerance, this should
  preferably match the tolerance used by the solver. The defualt tolerance should
  however apply in most situations (c.f. "Computational Techniques of the
  Simplex Method" by István Maros, section 9.3.4).
"""
function lp_objective_perturbation_range(var::VariableRef; optimality_tolerance::Float64 = 1e-8)::Tuple{Float64, Float64}
    model = owner_model(var)
    if termination_status(model) != MOI.OPTIMAL
        error("The perturbation range of the objective is not available because the current solution "*
              "is not optimal.")
    end
    if objective_sense(model) == MOI.FEASIBILITY_SENSE
        error("The perturbation range of the objective is not applicable on feasibility problems.")
    end
    if !has_duals(model)
        error("The perturbation range of the objective is not available because no dual result is " *
              "available.")
    end
    # TODO: This could be made more efficient for non-basic variables by checking them
    #       before building calling ´_std_matrix´ and ´_std_reduced_costs´ and only check
    #       the reduced cost of that variable. (issue 1892 makes this somewhat non-trivial)

    # The variable is in the basis, for a minimization problem without variable upper bounds we need:
    # c_N - N^T B^-T (c_B + Δ e_j ) >= 0
    # c_red - Δ N^T (B^-T e_j) = c_red - Δ N^T ρ = c_red - Δ N_red >= 0
    # c_red >= Δ N_red
    # maximum( c_red_- ./ N_red_- ) <= Δ <= minimum( c_red_+ ./ N_red_+ )
    # If upper bounds are present (and active), those inequalities are flipped.
    # If the problem is a maximization problem all inequalities are flipped.
    vars = all_variables(model)
    j = findfirst(isequal(var), vars)

    A, var_lower, var_upper, constraints = _std_matrix(model)
    var_basis_status = _std_basis_status(model)
    c_red = _std_reduced_costs(model, constraints)

    if var_basis_status[j] != MOI.BASIC
        if var_lower[j] == var_upper[j]
            return -Inf, Inf
        end
        if (var_basis_status[j] == MOI.NONBASIC_AT_LOWER && objective_sense(model) == MOI.MIN_SENSE) ||
            (var_basis_status[j] == MOI.NONBASIC_AT_UPPER && objective_sense(model) == MOI.MAX_SENSE)
            return -c_red[j], Inf
        end
        return -Inf, -c_red[j]
    end

    basic = var_basis_status .== Ref(MOI.BASIC)
    upper = var_basis_status[.!basic] .== Ref(MOI.NONBASIC_AT_UPPER)
    ej = zeros(size(A, 1))
    ej[findfirst(isequal(var), vars[basic[1:length(vars)]])] = 1.0
    N_red = A[:, .!basic]' * (A[:, basic]' \ ej)
    c_red = c_red[.!basic]
    pos = N_red .> 0.0
    neg = N_red .< 0.0

    in_lb = (neg .& .!upper) .| (pos .& upper)
    in_ub = (pos .& .!upper) .| (neg .& upper)
    if objective_sense(model) == MOI.MAX_SENSE
        in_lb, in_ub = in_ub, in_lb
    end
    # The reduced cost of variables fixed at equality do not be accounted for (they only change binding direction).
    unfixed_vars = (var_lower .< var_upper)[.!basic]
    in_lb .&= unfixed_vars
    in_ub .&= unfixed_vars
    # Find the first reduced cost that is strictly violated
    lower_bounds_delta = [-Inf; c_red[in_lb] ./ N_red[in_lb]]
    lower_bounds_delta_strict = lower_bounds_delta - [0.0; optimality_tolerance ./ abs.(N_red[in_lb])]
    upper_bounds_delta = [Inf; c_red[in_ub] ./ N_red[in_ub]]
    upper_bounds_delta_strict = upper_bounds_delta + [0.0; optimality_tolerance ./ abs.(N_red[in_ub])]

    lower_max_idx = last(findmax(lower_bounds_delta_strict))
    upper_min_idx = last(findmin(upper_bounds_delta_strict))
    return lower_bounds_delta[lower_max_idx], upper_bounds_delta[upper_min_idx]
end
