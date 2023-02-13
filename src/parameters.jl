import SparseArrays
import LinearAlgebra
import OrderedCollections
using JuMP


"""
    ParameterData

Representation of the data associated with a parameter (its `value`). Keeps
the incremental `index` to efficiently construct the affine update calculation.
"""
mutable struct ParameterData
    index::UInt64
    value::Float64
end

"""
    AffineVariableUpdate

An affine calculation of a variable's coefficient in a constraint.
"""
struct AffineVariableUpdate
    # This keeps the constant term.
    constant::Float64
    # This keeps the multiplicative factors for each param in the linear term.
    params::SparseArrays.SparseVector{Float64,UInt64}
end

"""
    AffineUpdateLink

A single update to a given `constraint`, calculating the coefficient of `variable`
following the affine update calculation represented by `update`
"""
struct AffineUpdateLink
    constraint::ConstraintRef
    variable::VariableRef
    update::AffineVariableUpdate
end

"""
    AffineObjectiveUpdateLink

A single update to the objective, calculating the coefficient of `variable`
following the affine update calculation represented by `update`
"""
struct AffineObjectiveUpdateLink
    variable::VariableRef
    update::AffineVariableUpdate
end

struct ParametricModelData
    parameters::OrderedCollections.OrderedDict{VariableRef,ParameterData}
    constraint_links::Vector{AffineUpdateLink}
    objective_links::Vector{AffineObjectiveUpdateLink}
end

"""
    ParametricConstraint

Parametric constraint type, that allows using parameters (in the form of VariableRef)
inside the constraints body. Keeps the affine expression, the set, as well as the
temporal (not bound to a ConstraintRef) `AffineVariableUpdate`s.
"""
struct ParametricConstraint{S} <: AbstractConstraint
    f::AffExpr
    s::S
    temporal_update_links::Dict{VariableRef,AffineVariableUpdate}
end

"""
    JuMP.build_constraint(..., ::GenericAffExpr, ..., ...)

Wrapper that enables constructing a parametric constraint that does not contain
any parameters (therefore constructing an AffExpr instead of a QuadExpr).
"""
function JuMP.build_constraint(
    _error::Function,
    f::GenericAffExpr,
    set::MOI.AbstractScalarSet,
    ::Type{ParametricConstraint};
)
    return JuMP.build_constraint(_error, f, set)
end

"""
    JuMP.build_constraint(..., ::GenericQuadExpr, ..., ...)

Convert the constructed QuadExpr into an AffExpr and prepare proper
updates to later calculate the coefficient of a given `VariableRef`
inside the constraint based on all parameter values.
"""
function JuMP.build_constraint(
    _error::Function,
    f::JuMP.GenericQuadExpr,
    set::MOI.AbstractScalarSet,
    ::Type{ParametricConstraint};
)
    affine = f.aff
    terms = f.terms

    model = nothing # todo: extract the model somehow so that it doesn't default to `nothing`
    n_param = 0

    # We expect as many updates as there are quadratic terms.
    updates = Dict{VariableRef,AffineVariableUpdate}()
    sizehint!(updates, length(terms))

    for (vars, coeff) in terms
        if isnothing(model)
            model = owner_model(vars.a)
            n_param = length(model.ext[:__parameters].parameters)
        end

        # Look up ParameterData for both involved variables.
        v_a = get(model.ext[:__parameters].parameters, vars.a, nothing)
        v_b = get(model.ext[:__parameters].parameters, vars.b, nothing)

        if v_a !== nothing
            if !haskey(updates, vars.b)
                # Construct a AffineVariableUpdate, based on the constant
                # (the affine term of the variabel) and a zeroed sparse vector.
                updates[vars.b] = AffineVariableUpdate(
                    get(affine.terms, vars.b, 0.0),
                    SparseArrays.spzeros(n_param),
                )
            end

            # Add the coefficient of the quadratic term to the entry of
            # the sparse vector corresponding to the correct parameter.
            # Basically, for `2.0 * p * x`, this remembers the mapping
            #   x -> p -> 2.0
            # later telling us that the coefficient of `x` can be
            # calculated as:
            #   constant + 2.0 * `p`
            updates[vars.b].params[v_a.index] += coeff
        elseif v_b !== nothing
            if !haskey(updates, vars.a)
                updates[vars.a] = AffineVariableUpdate(
                    get(affine.terms, vars.a, 0.0),
                    SparseArrays.spzeros(n_param),
                )
            end

            updates[vars.a].params[v_b.index] += coeff
        else
            # Both involved variables are "free" and not a parameter.
            _error("At least one parameter is not properly registered.")
        end
    end

    # The affine function may contain a constant term, move it to the set.
    new_set = typeof(set)(MOI.constant(set) - affine.constant)
    affine.constant = 0

    return ParametricConstraint(affine, new_set, updates)
end

"""
    JuMP.add_constraint(::Model, ::ParametricConstraint, ::String)

Add the given parametric constraint `con` to the model. After adding it,
this inserts the (previously only temporal) update links into the model.
"""
function JuMP.add_constraint(
    model::Model,
    con::ParametricConstraint,
    name::String,
)
    # todo: handle the name

    constr = add_constraint(model, ScalarConstraint(con.f, con.s))

    # Add all update links, now bound to the constructed constraint.
    for (k, v) in con.temporal_update_links
        push!(model.ext[:__parameters].constraint_links, AffineUpdateLink(constr, k, v))
    end

    return constr
end

"""
    _mmdot(x, y)

Compute the dot product of `x` and `y`, allowing for mismatches in the dimension.
The missing entries of the shorter vector are assumed to be zero.
"""
function _mmdot(x, y)
    lx = length(x)
    ly = length(y)

    if lx == ly
        return LinearAlgebra.dot(x, y)
    elseif lx < ly
        return LinearAlgebra.dot(x, @view(y[1:lx]))
    else
        return LinearAlgebra.dot(@view(x[1:ly]), y)
    end
end

"""
    _eval_var_coeff(update::AffineVariableUpdate, param_values::Vector{Float64})

Evaluate the affine update calculation given by `update` on the current
global state of all parameters (as defined by `param_values`)
"""
function _eval_var_coeff(
    update::AffineVariableUpdate,
    param_values::Vector{Float64},
)
    return update.constant + _mmdot(update.params, param_values)
end

"""
    _finalize_parameters(model::JuMP.Model)

Finalize all parameter values by updating the coefficients of each variable
that is subject to a registered `AffineVariableUpdate` in a `ParametricConstraint`.
"""
function _finalize_parameters(model::JuMP.Model)
    # Collect the values of all parameters into a vector. Since
    # `model.ext[:__parameters]` is an `OrderedDict`, this unambiguous.
    param_values = collect(it[2].value for it in model.ext[:__parameters].parameters)

    # Update all coefficients.
    @inbounds @simd for ul in model.ext[:__parameters].constraint_links
        # todo: somehow this results in less allocations but more time (be = backend(model) previously):
        # MOI.modify(
        #     be, ul.constraint,
        #     MOI.ScalarCoefficientChange(ul.variable, eval_var_coeff(model, ul.update))
        # )
        set_normalized_coefficient(
            ul.constraint,
            ul.variable,
            _eval_var_coeff(ul.update, param_values),
        )
    end
    # model.is_model_dirty = true       # this is used if MOI.modify(...) is used

    # Update the objective function
    @inbounds @simd for ul in model.ext[:__parameters].objective_links
        set_objective_coefficient(
            model,
            ul.variable,
            _eval_var_coeff(ul.update, param_values),
        )
    end

    return optimize!(model; ignore_optimize_hook = true)
end

"""
    enable_parameters!(model::JuMP.Model)

Enable parameters for `model` by registering the necessary optimize hook,
as well as creating the intermediate storages that are used by parameters.
"""
function enable_parameters!(model::JuMP.Model)
    haskey(model, :__parameters) && return nothing

    set_optimize_hook(model, _finalize_parameters)

    # todo: these could be part of `Model` (like `set_string_names_on_creation`)
    model.ext[:__parameters] = ParametricModelData(
        OrderedCollections.OrderedDict{VariableRef,ParameterData}(),
        Vector{AffineUpdateLink}(),
        Vector{AffineUpdateLink}()
    )

    return nothing
end

function _set_value(parameter::VariableRef, value::Float64; fix::Bool = true)
    # todo: this could be done only once if this is called by broadcasting
    model = owner_model(parameter)
    if !haskey(model.ext[:__parameters].parameter, parameter)
        error("Can not set value on a non-parameter.")
    end

    # Change the value of the parameter.
    model.ext[:__parameters].parameters[parameter].value = value

    # If we need fixing, update the constraint.
    fix && JuMP.fix(parameter, value; force = true)

    # todo: modify the objective function here to fix all parameter values
    return nothing
end

"""
    set_value(parameter::VariableRef, value::Float64; fix::Bool=true)

Set the value of `parameter` to `value`. If `fix` is `true`, this also
inserts a MOI.EqualTo, that fixes the value of the corresponding `VariableRef`.
Fixing decreases performance and is only necessary for parameters that are not
only used as multiplicative parameters (e.g. on the RHS of a constraint).
"""
function JuMP.set_value(
    parameter::VariableRef,
    value::Float64;
    fix::Bool = true,
)
    return _set_value(parameter, value; fix = fix)
end

"""
    @pexpression(...)

Construct a parametric expression.

Example usage:

    model = JuMP.Model(HiGHS.Optimizer)
    enable_parameters!(model)

    @variable(model, x)
    @parameter(model, p)

    @pexpression(model, e, x + 10)
    add_to_expression!(e, p * x)
"""
macro pexpression(args...)
    # todo: only to the wrapping if it is an AffExpr
    # todo: handle the case were the only arg is a VariableRef
    return esc(:($JuMP.QuadExpr.($JuMP.@expression($(args...)))))
end

"""
    @pconstraint(...)

Construct a parametric constraint.

Using this macro for a constraint makes sure that all involved
parameters are tracked accordingly and automatically update the constraint
as soon as `optimize!(...)` is called. Using this for a constraint that does
not contain any parameters does not involve a performance loss (since it
then just defaults to the standard `build_constraint(...)`); therefore it
can (mostly) be used liberally.

Before a call to `optimize!(...)`, this constraint does not represent
the correct left-hand-side expression. To force an update, `_finalize_parameters(model)`
can be used. This is not recommended!

Example usage:

    model = JuMP.Model(HiGHS.Optimizer)
    enable_parameters!(model)

    @variable(model, x)
    @parameter(model, p)

    @pconstraint(model, p * x >= 0)
"""
macro pconstraint(args...)
    return esc(:($JuMP.@constraint($(args...), ParametricConstraint)))
end

"""
    @parameter(model, ..., value; fix=true)

Construct a parameter (can be named and vectorized like a variable) and
set its initial `value`. If `fix` is `true`, this will create an additional
MOI.EqualTo constraint, that fixes the value of the corresponding `VariableRef`.
Fixing decreases performance and is only necessary for parameters that are not
only used as multiplicative parameters (e.g. on the RHS of a constraint).
"""
macro parameter(args...)
    # todo: this is performing worse than the manual code, improve this!
    _error(str...) = JuMP._macro_error(:parameter, args, __source__, str...)

    args = JuMP._reorder_parameters(args)
    flat_args, kw_args, _ = JuMP.Containers._extract_kw_args(args)
    kw_args = Dict(kw.args[1] => kw.args[2] for kw in kw_args)

    if !(length(flat_args) in [2, 3])
        _error(
            "Wrong number of arguments. Did you miss the initial parameter value?",
        )
    end
    # todo: properly catch misconfigurations of arguments
    # if !(flat_args[end] isa Number) and !(flat_args[end] isa Vector)
    #     _error("Wrong arguments. Did you miss the initial parameter value?")
    # end
    _vector_scalar_get(_scalar::Number, ::Int64) = _scalar
    function _vector_scalar_get(_vector::Vector{T} where {T<:Number}, i::Int64)
        return _vector[i]
    end

    # todo: fix this mess
    if get(kw_args, :fix, true)
        return esc(
            quote
                _var = $JuMP.@variable($(flat_args[1:(end-1)]...))
                if _var isa Vector
                    for i in eachindex(_var)
                        $(flat_args[1]).ext[:__parameters].parameters[_var[i]] =
                            ParameterData(
                                UInt64(
                                    length($(flat_args[1]).ext[:__parameters].parameters) +
                                    1,
                                ),
                                $_vector_scalar_get($(flat_args[end]), i),
                            )
                    end
                else
                    $(flat_args[1]).ext[:__parameters].parameters[_var] = ParameterData(
                        UInt64(length($(flat_args[1]).ext[:__parameters].parameters) + 1),
                        $(flat_args[end]),
                    )
                    $JuMP.fix(_var, $(flat_args[end]); force = true)
                end
                _var
            end,
        )
    else
        return esc(
            quote
                _var = $JuMP.@variable($(flat_args[1:(end-1)]...))
                if _var isa Vector
                    for i in eachindex(_var)
                        $(flat_args[1]).ext[:__parameters].parameters[_var[i]] =
                            ParameterData(
                                UInt64(
                                    length($(flat_args[1]).ext[:__parameters].parameters) +
                                    1,
                                ),
                                $_vector_scalar_get($(flat_args[end]), i),
                            )
                    end
                else
                    $(flat_args[1]).ext[:__parameters].parameters[_var] = ParameterData(
                        UInt64(length($(flat_args[1]).ext[:__parameters].parameters) + 1),
                        $(flat_args[end]),
                    )
                end
                _var
            end,
        )
    end
end

# ==================
# the following is not properly documented since there is probably a
# better way to hook into objective creation
# 
# additionally, we should use `@pobjective`, and that should handle
# converting (e.g.) from `Min` to `ParametricMin` automatically

struct ParametricMin <: MOI.AbstractModelAttribute end
struct ParametricMax <: MOI.AbstractModelAttribute end

function JuMP._moi_sense(_error::Function, sense)
    if sense == :Min
        expr = MIN_SENSE
    elseif sense == :Max
        expr = MAX_SENSE
    else
        # Refers to a variable that holds the sense.
        # TODO: Better document this behavior
        expr = esc(sense)
    end
    # todo: catch invalid senses again and add the sense properly above
    return :($expr)# :(_throw_error_for_invalid_sense($_error, $expr))
end

function _prepare_parametric_objective(model::AbstractModel, func)
    affine = func.aff
    terms = func.terms

    n_param = length(model.ext[:__parameters].parameters)

    # We expect as many updates as there are quadratic terms.
    updates = Dict{VariableRef,AffineVariableUpdate}()
    sizehint!(updates, length(terms))

    for (vars, coeff) in terms
        # Look up ParameterData for both involved variables.
        v_a = get(model.ext[:__parameters].parameters, vars.a, nothing)
        v_b = get(model.ext[:__parameters].parameters, vars.b, nothing)

        if v_a !== nothing
            if !haskey(updates, vars.b)
                # Construct a AffineVariableUpdate, based on the constant
                # (the affine term of the variabel) and a zeroed sparse vector.
                updates[vars.b] = AffineVariableUpdate(
                    get(affine.terms, vars.b, 0.0),
                    SparseArrays.spzeros(n_param),
                )
            end

            # Add the coefficient of the quadratic term to the entry of
            # the sparse vector corresponding to the correct parameter.
            # Basically, for `2.0 * p * x`, this remembers the mapping
            #   x -> p -> 2.0
            # later telling us that the coefficient of `x` can be
            # calculated as:
            #   constant + 2.0 * `p`
            updates[vars.b].params[v_a.index] += coeff
        elseif v_b !== nothing
            if !haskey(updates, vars.a)
                updates[vars.a] = AffineVariableUpdate(
                    get(affine.terms, vars.a, 0.0),
                    SparseArrays.spzeros(n_param),
                )
            end

            updates[vars.a].params[v_b.index] += coeff
        else
            # Both involved variables are "free" and not a parameter.
            _error("At least one parameter is not properly registered.")
        end
    end

    empty!(model.ext[:__parameters].objective_links)
    for (k, v) in updates
        push!(model.ext[:__parameters].objective_links, AffineObjectiveUpdateLink(k, v))
    end

    return set_objective_function(model, affine)
end

function _parametric_min(model::AbstractModel, ::Type{ParametricMin}, func)
    set_objective_sense(model, MOI.MIN_SENSE)
    _prepare_parametric_objective(model, func)
end

function _parametric_max(model::AbstractModel, ::Type{ParametricMax}, func)
    set_objective_sense(model, MOI.MAX_SENSE)
    _prepare_parametric_objective(model, func)
end

JuMP.set_objective(model::AbstractModel, ::Type{ParametricMin}, func) = _parametric_min(model, ParametricMin, func)
JuMP.set_objective(model::AbstractModel, ::Type{ParametricMax}, func) = _parametric_max(model, ParametricMax, func)