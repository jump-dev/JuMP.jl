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
    perturbation_range_of_feasibility(constraint::ConstraintRef)

Gives the range by which the rhs coefficient can change and the current lp-basis
remains feasible, i.e., where the shadow prices apply.

## Notes
- The range denote valid changes, e.g., for a*x <= b + Δ, the lp-basis remains feasible for all Δ in [l, u].

"""
function perturbation_range_of_feasibility(constraint::ConstraintRef{Model, <:_MOICON})
    error("The range of feasibility is not defined or not implemented for this type " *
          "of constraint.")
end

"""
Builds the standard form constraint matrix, varaible bounds (rhs := 0, by choice of slack bounds), and a vector of constraint reference one for each row.
i.e., the LP-problem feasible set on the form
s.t.  Ax == 0
L <= x <= U
"""
function _std_matrix(model::Model)
    con_types = list_of_constraint_types(model)
    con_func_type = first.(con_types)
    if !all(broadcast(<:, AffExpr, con_func_type) .| broadcast(<:, VariableRef, con_func_type))
        error("Range of feasibility and optimality require all constraints functions to be affine (the problem to be linear).")
    end
    vars = all_variables(model)
    var2idx = Dict{VariableRef, Int64}()
    for (j, var) in pairs(vars)
        var2idx[var] = j
    end
    col_j = Int64[]
    row_i = Int64[]
    coeffs = Float64[]

    var_lower = fill(-Inf, length(vars))
    var_upper = fill(Inf, length(vars))
    affine_constraints = ConstraintRef[]

    m = 0
    n = length(vars)
    for (F, S) in con_types
        constraints = all_constraints(model, F, S)
        if F <: VariableRef
            for (i, constr) in pairs(constraints)
                constr_obj = constraint_object(constr)
                j = var2idx[constr_obj.func]
                range = MOI.Interval(constr_obj.set)
                var_lower[j] = range.lower
                var_upper[j] = range.upper
            end
            continue
        end
        for (i, constr) in pairs(constraints)
            constr_obj = constraint_object(constr)
            constr_terms = constr_obj.func.terms
            m += 1
            n += 1
            push!(affine_constraints, constr)
            for (var, coeff) in constr_terms
                push!(col_j, var2idx[var])
                push!(row_i, m)
                push!(coeffs, coeff)
            end
            push!(col_j, n)
            push!(row_i, m)
            push!(coeffs, -1.0)

            range = MOI.Interval(constr_obj.set)
            push!(var_lower, range.lower)
            push!(var_upper, range.upper)
        end
    end

    A = sparse(row_i, col_j, coeffs, m, n)
    return A, var_lower, var_upper, affine_constraints
end
"""
Returns the optimal primal values for the problem in standard form, c.f. `_std_matrix(model)`.
"""
function _std_primal_solution(model::Model, constraints::Vector{ConstraintRef})
    vars = all_variables(model)
    x = [value.(vars); value.(constraints)]
    return x
end


"""
Returns the basis status of all variables of the problem in standard form, c.f. `_std_matrix(model)`.
"""
function _std_basis_status(model::Model)::Vector{MOI.BasisStatusCode}
    con_types = list_of_constraint_types(model)
    vars = all_variables(model)
    var2idx = Dict{VariableRef, Int64}()
    for (j, var) in pairs(vars)
        var2idx[var] = j
    end
    basis_status = fill(MOI.BASIC, length(vars))

    for (F, S) in con_types
        constraints = all_constraints(model, F, S)
        if F <: VariableRef
            for (i, constr) in pairs(constraints)
                constr_obj = constraint_object(constr)
                j = var2idx[constr_obj.func]
                basis_status[j] = MOI.get(model, MOI.ConstraintBasisStatus(), constr)
                if S <: MOI.LessThan && basis_status[j] == MOI.NONBASIC
                    basis_status[j] = MOI.NONBASIC_AT_UPPER
                elseif S <: MOI.GreaterThan && basis_status[j] == MOI.NONBASIC
                    basis_status[j] = MOI.NONBASIC_AT_LOWER
                end
            end
            #issue 1892 makes is hard to know how to handle variable constraints, atm. assume unique variable constraints.
            continue
        end
        for constr in constraints
            con_basis_status = MOI.get(model, MOI.ConstraintBasisStatus(), constr)
            if S <: MOI.LessThan && con_basis_status == MOI.NONBASIC
                con_basis_status = MOI.NONBASIC_AT_UPPER
            elseif S <: MOI.GreaterThan && con_basis_status == MOI.NONBASIC
                con_basis_status = MOI.NONBASIC_AT_LOWER
            end
            push!(basis_status, con_basis_status)
        end
    end
    return basis_status
end

function perturbation_range_of_feasibility(constraint::ConstraintRef{Model, _MOICON{F, S}}
                      ) where {S <: Union{MOI.LessThan, MOI.GreaterThan, MOI.EqualTo}, T, F <: Union{MOI.ScalarAffineFunction{T}, MOI.SingleVariable}}
    model = constraint.model
    if termination_status(model) != MOI.OPTIMAL
        error("The range of feasibility is not available because the current solution "*
              "is not optimal.")
    end

    # This throw the error: ... do not support accessing ... ConstraintBasisStatus(), if its not implemented by the solver.
    con_status = MOI.get(model, MOI.ConstraintBasisStatus(), constraint)
    # The constraint is inactive.
    if con_status == MOI.BASIC
        if S <: MOI.LessThan
            return (value(constraint) - constraint_object(constraint).set.upper, Inf)
        end
        if S <: MOI.GreaterThan
            return (- Inf, value(constraint) - constraint_object(constraint).set.lower)
        end
        return (0.0, 0.0)
    end
    # The constraint is active.

    # Optimization possible, only basic columns needed
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
    basic = [var_basis_status[i] == MOI.BASIC for i in 1:length(var_basis_status)] # .== throws ERROR: MethodError: no method matching length(::MathOptInterface.BasisStatusCode)
    B = A[:, basic]
    x_B = x[basic]
    local rho::Vector{Float64}
    if F <: MOI.SingleVariable
        var = constraint_object(constraint).func
        j = findfirst(all_variables(model) .== var)
        rho = -B \ Vector(A[:, j])
    else
        i = findfirst(constraints .== constraint)
        ei = zeros(size(A)[1])
        ei[i] = 1.0
        rho = B \ ei
    end
    UmX_B = U[basic] - x_B
    LmX_B = L[basic] - x_B
    pos = rho .> 1e-7
    neg = rho .< -1e-7
    return (maximum([-Inf; UmX_B[neg] ./ rho[neg]; LmX_B[pos] ./ rho[pos]]),
            minimum([Inf; UmX_B[pos] ./ rho[pos]; LmX_B[neg] ./ rho[neg]]))
end

"""
Returns the optimal reduced costs for the variables of the problem in standard form, c.f. `_std_matrix(model)`
"""
function _std_reduced_costs(model::Model, constraints::Vector{ConstraintRef})
    con_types = list_of_constraint_types(model)
    vars = all_variables(model)
    var2idx = Dict{VariableRef, Int64}()
    for (j, var) in pairs(vars)
        var2idx[var] = j
    end
    y = [zeros(length(vars)); dual.(constraints)]

    for (F, S) in con_types
        if F <: VariableRef
            constraints = all_constraints(model, F, S)
            for (i, constr) in pairs(constraints)
                constr_obj = constraint_object(constr)
                j = var2idx[constr_obj.func]
                y[j] += dual(constr)
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
    perturbation_range_of_optimality(var::VariableRef)

Gives the range by which the cost coefficient can change and the current lp-basis
remains optimal, i.e., the reduced costs remain valid.

## Notes
- The range denote valid changes, Δ in [l, u], for which c[var] += Δ do not violate the current optimality conditions.

"""
function perturbation_range_of_optimality(var)
    model = var.model
    if termination_status(model) != MOI.OPTIMAL
        error("The optimality range is not available because the current solution "*
              "is not optimal.")
    end
    if objective_sense(model) == MOI.FEASIBILITY_SENSE
        error("The optimality range is not applicable on feasibility problems.")
    end
    if !has_duals(model)
        error("The optimality range is not available because no dual result is " *
              "available.")
    end
    # Optimization possible: since the has_upper_bound(var) and has_lower_bound(var)
    # Do not capture all variable bounds (issue 1892), it's hard to directly
    # access the reduced cost of a variable, now the complete matrix is always built.

    # The variable is in the basis for a minimization problem with no variable upper bounds we need:
    # c_N - N^T B^-T (c_B + Δ e_j ) >= 0
    # c_red - Δ N^T (B^-T e_j) = c_red - Δ N^T ρ = c_red - Δ N_red >= 0
    # c_red >= Δ N_red
    # c_red >= Δ N_red
    # maximum( c_red_- ./ N_red_- ) <= Δ <= minimum( c_red_+ ./ N_red_+ )
    # If upper bounds are present (and active), those inequalities are flipped.
    # If the problem is a maximization problem all inequalities are flipped.
    vars = all_variables(model)
    j = findfirst(vars .== var)

    A, L, U, constraints = _std_matrix(model)
    var_basis_status = _std_basis_status(model)
    c_red = _std_reduced_costs(model, constraints)

    if var_basis_status[j] != MOI.BASIC
        if L[j] == U[j]
            return (-Inf, Inf)
        end
        if (var_basis_status[j] == MOI.NONBASIC_AT_LOWER && objective_sense(model) == MOI.MIN_SENSE) ||
            (var_basis_status[j] == MOI.NONBASIC_AT_UPPER && objective_sense(model) == MOI.MAX_SENSE)
            return (-c_red[j], Inf)
        end
        return (-Inf, -c_red[j])
    end

    basic = [var_basis_status[i] == MOI.BASIC for i in 1:length(var_basis_status)] # .== throws ERROR: MethodError: no method matching length(::MathOptInterface.BasisStatusCode)
    upper = [var_basis_status[.!basic][i] == MOI.NONBASIC_AT_UPPER for i in 1:length(var_basis_status[.!basic])] # .== throws ERROR: MethodError: no method matching length(::MathOptInterface.BasisStatusCode)
    ej = zeros(size(A)[1])
    ej[findfirst(vars[basic[1:length(vars)]] .== var)] = 1.0
    N_red = A[:, .!basic]' * (A[:, basic]' \ ej)
    c_red = c_red[.!basic]
    pos = N_red .> 1e-7
    neg = N_red .< -1e-7
    in_lb = (neg .& .!upper) .| (pos .& upper)
    in_ub = (pos .& .!upper) .| (neg .& upper)
    if objective_sense(model) == MOI.MAX_SENSE
        in_lb, in_ub = in_ub, in_lb
    end
    unfixed_vars = (L .< U)[.!basic]
    in_lb .&= unfixed_vars
    in_ub .&= unfixed_vars
    return (maximum([-Inf; c_red[in_lb] ./ N_red[in_lb]]),
            minimum([Inf; c_red[in_ub] ./ N_red[in_ub]]))
end
