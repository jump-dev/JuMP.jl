struct ParametricTerm{T}
    var::MOI.VariableIndex
    coeff::T
    param::MOI.VariableIndex
end

struct FixParametricVariablesBridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    aff_constr::MOI.ConstraintIndex
    parametric_terms::Vector{ParametricTerm{T}}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:FixParametricVariablesBridge}, 
    F::Type{<:MOI.ScalarQuadraticFunction},
    S::Type{<:MOI.AbstractScalarSet},
)
    return FixParametricVariablesBridge{Float64}
end

function MOI.supports_constraint(
    ::Type{<:FixParametricVariablesBridge},
	::Type{<:MOI.ScalarQuadraticFunction},
	::Type{<:MOI.AbstractScalarSet},
)
	return true
end

function MOI.Bridges.added_constrained_variable_types(::Type{<:FixParametricVariablesBridge}) 
	return Tuple{DataType}[]    # todo
end

function MOI.Bridges.added_constraint_types(::Type{<:FixParametricVariablesBridge})
	return []                   # todo
end

MOI.Bridges.needs_final_touch(::FixParametricVariablesBridge) = true

function MOI.Bridges.final_touch(
    bridge::FixParametricVariablesBridge{T},
    model::MOI.ModelLike,
) where {T}
    for _t in bridge.touches
        coeff = _t.coeff
        fix_constr = MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}}(_t.param.value)
        curr_param_val = MOI.get(model, MOI.ConstraintSet(), fix_constr).value
        attached_var = _t.var
        
        MOI.modify(
            model,
            bridge.aff_constr,
            MOI.ScalarCoefficientChange(attached_var, convert(T, coeff * curr_param_val)),
        )
    end
    # todo: do we need to set the model to "dirty"?
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{<:FixParametricVariablesBridge}, 
    m::MOI.ModelLike, 
    f::MOI.ScalarQuadraticFunction, 
    s::MOI.AbstractScalarSet,       # todo: Union(LT, GT, EQ)
)
    affine_terms = f.affine_terms
    constant = f.constant

    parametric_terms = []
    for term in f.quadratic_terms
        coeff = term.coefficient

        c1 = MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}}(term.variable_1.value)
        if MOI.is_valid(m, c1)
            push!(parametric_terms, Touch(term.variable_2, coeff, term.variable_1))
            coeff *= MOI.get(m, MOI.ConstraintSet(), c1).value
            push!(affine_terms, MOI.ScalarAffineTerm(coeff, term.variable_2))
            continue
        end

        c2 = MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}}(term.variable_2.value)
        if MOI.is_valid(m, c2)
            push!(parametric_terms, Touch(term.variable_1, coeff, term.variable_2))
            coeff *= MOI.get(m, MOI.ConstraintSet(), c2).value
            push!(affine_terms, MOI.ScalarAffineTerm(coeff, term.variable_1))
            continue
        end
        
        @error "nope"
    end

    f = MOI.ScalarAffineFunction(affine_terms, constant)
    aff_constr = MOI.add_constraint(m, f, s)

    return FixParametricVariablesBridge{Float64}(aff_constr, parametric_terms)
end

function wrap_model!(model)
    add_bridge(model, FixParametricVariablesBridge)
end