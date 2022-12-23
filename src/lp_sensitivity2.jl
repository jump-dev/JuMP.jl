#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    SensitivityReport

See [`lp_sensitivity_report`](@ref).
"""
struct SensitivityReport{T}
    rhs::Dict{ConstraintRef,Tuple{T,T}}
    objective::Dict{GenericVariableRef{T},Tuple{T,T}}
end

Base.getindex(s::SensitivityReport, c::ConstraintRef) = s.rhs[c]
Base.getindex(s::SensitivityReport, x::GenericVariableRef) = s.objective[x]

"""
    lp_sensitivity_report(model::GenericModel{T}; atol::T = Base.rtoldefault(T))::SensitivityReport{T} where {T}

Given a linear program `model` with a current optimal basis, return a
[`SensitivityReport`](@ref) object, which maps:

 - Every variable reference to a tuple `(d_lo, d_hi)::Tuple{T,T}`,
   explaining how much the objective coefficient of the corresponding variable
   can change by, such that the original basis remains optimal.
 - Every constraint reference to a tuple `(d_lo, d_hi)::Tuple{T,T}`,
   explaining how much the right-hand side of the corresponding constraint can
   change by, such that the basis remains optimal.

Both tuples are relative, rather than absolute. So given a objective coefficient
of `1.0` and a tuple `(-0.5, 0.5)`, the objective coefficient can range
between `1.0 - 0.5` an `1.0 + 0.5`.

`atol` is the primal/dual optimality tolerance, and should match the tolerance
of the solver used to compute the basis.

Note: interval constraints are NOT supported.

## Example

```jldoctest
julia> model = Model(HiGHS.Optimizer);

julia> set_silent(model)

julia> @variable(model, -1 <= x <= 2)
x

julia> @objective(model, Min, x)
x

julia> optimize!(model)

julia> report = lp_sensitivity_report(model; atol = 1e-7);

julia> dx_lo, dx_hi = report[x]
(-1.0, Inf)

julia> println(
           "The objective coefficient of `x` can decrease by \$dx_lo or " *
           "increase by \$dx_hi."
       )
The objective coefficient of `x` can decrease by -1.0 or increase by Inf.

julia> dRHS_lo, dRHS_hi = report[LowerBoundRef(x)]
(-Inf, 3.0)

julia> println(
           "The lower bound of `x` can decrease by \$dRHS_lo or increase " *
           "by \$dRHS_hi."
       )
The lower bound of `x` can decrease by -Inf or increase by 3.0.
```
"""
function lp_sensitivity_report(
    model::GenericModel{T};
    atol::T = Base.rtoldefault(T),
) where {T}
    if !_is_lp(model)
        error(
            "Unable to compute LP sensitivity because model is not a linear " *
            "program (or it contains interval constraints).",
        )
    elseif !has_values(model)
        error("Unable to compute LP sensitivity: no primal solution available.")
    elseif !has_duals(model)
        error("Unable to compute LP sensitivity: no dual solution available.")
    end

    std_form = _standard_form_matrix(model)
    basis = _standard_form_basis(model, std_form)
    B = std_form.A[:, basis.basic_cols]
    if size(B, 1) != size(B, 2)
        error(
            "Unable to compute LP sensitivity: problem is degenerate. Try " *
            "adding variable bounds to free variables",
        )
    end

    n = length(std_form.columns)
    is_min = objective_sense(model) == MIN_SENSE

    x = vcat(value.(all_variables(model)), value.(std_form.constraints))
    x_B = @view x[basis.basic_cols]
    l_B = @view std_form.lower[basis.basic_cols]
    u_B = @view std_form.upper[basis.basic_cols]

    B_fact = if size(B, 1) > 0
        LinearAlgebra.lu(B)
    else
        zeros(Float64, (0, 0))
    end
    d = Dict{Int,Vector{T}}(
        # We call `collect` here because some Julia versions are missing sparse
        # matrix \ sparse vector fallbacks.
        j => B_fact \ collect(std_form.A[:, j]) for
        j in 1:length(basis.basic_cols) if basis.basic_cols[j] == false
    )

    report = SensitivityReport(
        Dict{ConstraintRef,Tuple{T,T}}(),
        Dict{GenericVariableRef{T},Tuple{T,T}}(),
    )

    ###
    ### Compute RHS sensitivity
    ###

    # There is an easy case to consider: a constraint is basic, so we can just
    # take the distance between the value of the constraint and the
    # corresponding bound. Otherwise, we need to compute a search direction as
    # in `_compute_rhs_range`. This is just the negative of the search direction
    # computed above. Moreover, we have to be careful with doubly-bounded
    # variables, because our computed range doesn't take into account the
    # inactive bound.

    for (i, con) in enumerate(std_form.constraints)
        if basis.constraints[i] == MOI.BASIC
            report.rhs[con] = _basic_range(con, constraint_object(con).set)
        else
            report.rhs[con] = _compute_rhs_range(-d[i+n], x_B, l_B, u_B, atol)
        end
    end
    for (i, con) in enumerate(std_form.bounds)
        con_obj = constraint_object(con)
        if basis.bounds[i] == MOI.BASIC
            report.rhs[con] = _basic_range(con, con_obj.set)
        else
            col = std_form.columns[con_obj.func]
            t_lo, t_hi = _compute_rhs_range(-d[col], x_B, l_B, u_B, atol)
            if basis.bounds[i] == MOI.NONBASIC_AT_UPPER
                t_lo = max(t_lo, std_form.lower[col] - x[col])
            elseif basis.bounds[i] == MOI.NONBASIC_AT_LOWER
                t_hi = min(t_hi, std_form.upper[col] - x[col])
            end
            report.rhs[con] = (t_lo, t_hi)
        end
    end

    ###
    ### Compute objective sensitivity
    ###

    π = Dict{Int,T}(
        i => reduced_cost(var) for
        (var, i) in std_form.columns if basis.variables[i] != MOI.BASIC
    )
    for (i, c) in enumerate(std_form.constraints)
        if basis.constraints[i] != MOI.BASIC
            π[n+i] = is_min ? dual(c) : -dual(c)
        end
    end

    for (var, i) in std_form.columns
        if basis.variables[i] == MOI.BASIC
            # The variable `i` is basic. Given an optimal basis B, the reduced
            # costs are:
            #   c_bar = π = c_N - c_Bᵀ(B⁻¹N)
            # To maintain optimality, we want to find a δ such that (if
            # minimizing):
            #     c_N - (c_B + δeᵢ)ᵀ(B⁻¹N) ≥ 0
            #     c_N - c_BᵀB⁻¹N - δ(eᵢ)ᵀ(B⁻¹N) ≥ 0
            #     π_N ≥ δ * (eᵢ)ᵀ(B⁻¹N)
            # To do so, we can loop through every nonbasic variable `j`, and
            # compute
            #     dᵢⱼ = (eᵢ)ᵀB⁻¹aⱼ
            # Then, depending on the sign of dᵢⱼ, we can compute bounds on δ.
            t_lo, t_hi = -Inf, Inf
            e_i = sum(basis.basic_cols[ii] for ii in 1:i)
            for j in 1:length(basis.basic_cols)
                if basis.basic_cols[j]
                    continue  # Ignore basic components.
                elseif isapprox(
                    std_form.lower[j],
                    std_form.upper[j];
                    atol = atol,
                )
                    continue  # Fixed variables can be ignored.
                elseif abs(d[j][e_i]) <= atol
                    continue  # Direction is ≈0 so any value of δ is okay.
                end
                # There are three confusing sign switch opportunities. We
                # usually want to be:
                # - minimizing (switch if maximizing)
                # - with d[j][e_i] ≥ 0 (switch if ≤ 0)
                # - and nonbasic at the lower bound (switch if upper)
                # If an odd number of these things is true, then the ratio
                # forms an upper bound for δ. Otherwise, it forms a lower bound.
                stat = j <= n ? basis.variables[j] : basis.constraints[j-n]
                if xor(is_min, d[j][e_i] > atol, stat == MOI.NONBASIC_AT_LOWER)
                    t_hi = min(t_hi, π[j] / d[j][e_i])
                else
                    t_lo = max(t_lo, π[j] / d[j][e_i])
                end
            end
            report.objective[var] = (t_lo, t_hi)
        elseif std_form.lower[i] == -Inf && std_form.upper[i] == Inf
            # The variable is nonbasic with free bounds.
            report.objective[var] = (0.0, 0.0)
        elseif std_form.lower[i] == std_form.upper[i]
            # The VariableIndex-in-EqualTo case.
            # The variable is nonbasic with fixed bounds. Therefore,
            # (δ⁻, δ⁺) = (-∞, ∞) because the variable can be effectively
            # substituted out.
            report.objective[var] = (-Inf, Inf)
        elseif basis.variables[i] == MOI.NONBASIC_AT_LOWER
            # The VariableIndex-in-GreaterThan case.
            # Variable `i` is nonbasic at lower bound. If minimizing, (δ⁻, δ⁺) =
            # (-πᵢ, ∞) because increasing the objective coefficient will only
            # keep it at the bound. If maximizing, the opposite is true.
            report.objective[var] = is_min ? (-π[i], Inf) : (-Inf, -π[i])
        else
            # The VariableIndex-in-LessThan case. Because we don't support
            # interval constraints, this assertion must hold.
            @assert basis.variables[i] == MOI.NONBASIC_AT_UPPER
            # Variable `i` is nonbasic at upper bound. The opposite case of the
            # one above.
            report.objective[var] = is_min ? (-Inf, -π[i]) : (-π[i], Inf)
        end
    end
    return report
end

_basic_range(con, set::MOI.LessThan) = (value(con) - set.upper, Inf)
_basic_range(con, set::MOI.GreaterThan) = (-Inf, value(con) - set.lower)
_basic_range(::Any, ::Any) = (0.0, 0.0)

"""
    _compute_rhs_range(d_B, x_B, l_B, u_B, atol)

Assume we start with the optimal solution `x_old`, we want to compute a step
size `t` in a direction `d` such that `x_new = x_old + t * d` is still
represented by the same optimal basis. This can be computed a la primal simplex
where we use an artificial entering variable.

    A * x_new = A * (x_old + t * d)
                = A * x_old + t * A * d
                = 0         + t * A * d  # Since A * x_old = 0
    =>  A * d = 0
    => B * d_B + N * d_N = 0
    => d_B = B \\ -(N * d_N)

Note we only have to compute the basic component of the direction vector,
because `d_N` is just zeros with a `1` in the component associated with the
artificial entering variable. Therefore, all that remains is to compute the
associated column of `N`.

If we are increasing the bounds associated with the `i`th decision variable,
then our artificial entering variable is a duplicate of the `i`th variable, and
`N * d_N = A[:, i]`.

If we are increasing the bounds associated with the `i`th affine constraint,
then our artificial entering variable is a duplicate of the slack variable
associated with the `i`th constraint, i.e., a `-1` in the `i`th row and zeros
everywhere else.

In either case:

    d_B = -(B \\ A[:, i])

Now, having computed a direction such that `x_new = x_old + t * d`. By ensuring
that `A * d = 0`, we maintained structural feasibility. Now we need to compute
bounds on `t` such that `x_new` maintains bound feasibility. That is, compute
bounds on t such that:

    l_B[j] <= x_B[j] + t * d_B[j] <= u_B[j].
"""
function _compute_rhs_range(d_B, x_B, l_B, u_B, atol)
    t_lo, t_hi = -Inf, Inf
    for j in 1:length(l_B)
        if d_B[j] > atol
            t_lo = max(t_lo, (l_B[j] - x_B[j]) / d_B[j])
            t_hi = min(t_hi, (u_B[j] - x_B[j]) / d_B[j])
        elseif d_B[j] < -atol
            t_lo = max(t_lo, (u_B[j] - x_B[j]) / d_B[j])
            t_hi = min(t_hi, (l_B[j] - x_B[j]) / d_B[j])
        else
            continue  # d_B[j] ≈ 0.0
        end
    end
    return t_lo, t_hi
end

"""
    _is_lp(model::GenericModel)

Return `true` if `model` is a linear program.
"""
function _is_lp(model::GenericModel)
    for (F, S) in list_of_constraint_types(model)
        # TODO(odow): support Interval constraints.
        if !(S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo})
            return false
        elseif !(F <: Union{GenericVariableRef,GenericAffExpr})
            return false
        end
    end
    return true
end

"""
    _standard_form_matrix(model::GenericModel)

Given a problem:

    r_l <= Ax <= r_u
    c_l <=  x <= c_u

Return the standard form:

           [A -I] [x, y] = 0
    [c_l, r_l] <= [x, y] <= [c_u, r_u]

`columns` maps the variable references to column indices.
"""
function _standard_form_matrix(model::GenericModel{T}) where {T}
    columns = Dict(var => i for (i, var) in enumerate(all_variables(model)))
    n = length(columns)
    c_l, c_u = fill(typemin(T), n), fill(typemax(T), n)
    r_l, r_u = T[], T[]
    I, J, V = Int[], Int[], T[]
    bound_constraints = ConstraintRef[]
    affine_constraints = ConstraintRef[]
    for (F, S) in list_of_constraint_types(model)
        _fill_standard_form(
            model,
            columns,
            bound_constraints,
            affine_constraints,
            F,
            S,
            c_l,
            c_u,
            r_l,
            r_u,
            I,
            J,
            V,
        )
    end
    return (
        columns = columns,
        lower = vcat(c_l, r_l),
        upper = vcat(c_u, r_u),
        A = SparseArrays.sparse(I, J, V, length(r_l), n + length(r_l)),
        bounds = bound_constraints,
        constraints = affine_constraints,
    )
end

function _fill_standard_form(
    model::GenericModel{T},
    x::Dict{GenericVariableRef{T},Int},
    bound_constraints::Vector{ConstraintRef},
    ::Vector{ConstraintRef},
    F::Type{GenericVariableRef{T}},
    S::Type,
    c_l::Vector{T},
    c_u::Vector{T},
    ::Vector{T},
    ::Vector{T},
    ::Vector{Int},
    ::Vector{Int},
    ::Vector{T},
) where {T}
    for c in all_constraints(model, F, S)
        push!(bound_constraints, c)
        c_obj = constraint_object(c)
        i = x[c_obj.func]
        set = MOI.Interval(c_obj.set)
        c_l[i] = max(c_l[i], set.lower)
        c_u[i] = min(c_u[i], set.upper)
    end
    return
end

function _fill_standard_form(
    model::GenericModel{T},
    x::Dict{GenericVariableRef{T},Int},
    ::Vector{ConstraintRef},
    affine_constraints::Vector{ConstraintRef},
    F::Type{<:GenericAffExpr},
    S::Type,
    ::Vector{T},
    ::Vector{T},
    r_l::Vector{T},
    r_u::Vector{T},
    I::Vector{Int},
    J::Vector{Int},
    V::Vector{T},
) where {T}
    for c in all_constraints(model, F, S)
        push!(affine_constraints, c)
        c_obj = constraint_object(c)
        @assert iszero(c_obj.func.constant)
        row = length(r_l) + 1
        set = MOI.Interval(c_obj.set)
        push!(r_l, set.lower)
        push!(r_u, set.upper)
        for (var, coef) in c_obj.func.terms
            push!(I, row)
            push!(J, x[var])
            push!(V, coef)
        end
        push!(I, row)
        push!(J, length(x) + row)
        push!(V, -one(T))
    end
    return
end

_convert_nonbasic_status(::MOI.LessThan) = MOI.NONBASIC_AT_UPPER
_convert_nonbasic_status(::MOI.GreaterThan) = MOI.NONBASIC_AT_LOWER
_convert_nonbasic_status(::Any) = MOI.NONBASIC

function _try_get_constraint_basis_status(model::GenericModel, constraint)
    try
        return MOI.get(model, MOI.ConstraintBasisStatus(), constraint)
    catch
        error(
            "Unable to query LP sensitivity information because this solver " *
            "does not support querying the status of constraints in the " *
            "optimal basis.",
        )
    end
end

function _try_get_variable_basis_status(model::GenericModel, variable)
    try
        return MOI.get(model, MOI.VariableBasisStatus(), variable)
    catch
        error(
            "Unable to query LP sensitivity information because this solver " *
            "does not support querying the status of variables in the " *
            "optimal basis.",
        )
    end
end

_nonbasic_at_lower(::MOI.GreaterThan) = MOI.NONBASIC_AT_LOWER
_nonbasic_at_lower(::Any) = MOI.BASIC
_nonbasic_at_upper(::MOI.LessThan) = MOI.NONBASIC_AT_UPPER
_nonbasic_at_upper(::Any) = MOI.BASIC

function _standard_form_basis(model::GenericModel, std_form)
    variable_status = fill(MOI.BASIC, length(std_form.columns))
    bound_status = fill(MOI.BASIC, length(std_form.bounds))
    constraint_status = fill(MOI.BASIC, length(std_form.constraints))
    for (x, col) in std_form.columns
        variable_status[col] = _try_get_variable_basis_status(model, x)
    end
    for (i, c) in enumerate(std_form.bounds)
        c_obj = constraint_object(c)
        col = std_form.columns[c_obj.func]
        status = variable_status[col]
        if status == MOI.NONBASIC_AT_LOWER
            bound_status[i] = _nonbasic_at_lower(c_obj.set)
        elseif status == MOI.NONBASIC_AT_UPPER
            bound_status[i] = _nonbasic_at_upper(c_obj.set)
        elseif status == MOI.NONBASIC
            bound_status[i] = _convert_nonbasic_status(c_obj.set)
        else
            @assert status == MOI.BASIC
            bound_status[i] = MOI.BASIC
        end
    end
    for (i, c) in enumerate(std_form.constraints)
        status = _try_get_constraint_basis_status(model, c)
        if status == MOI.NONBASIC
            status = _convert_nonbasic_status(constraint_object(c).set)
        end
        constraint_status[i] = status
    end
    return (
        variables = variable_status,
        bounds = bound_status,
        constraints = constraint_status,
        basic_cols = [variable_status; constraint_status] .== Ref(MOI.BASIC),
    )
end
