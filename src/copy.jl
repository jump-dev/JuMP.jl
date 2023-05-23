#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    copy_extension_data(data, new_model::AbstractModel, model::AbstractModel)

Return a copy of the extension data `data` of the model `model` to the extension
data of the new model `new_model`.

A method should be added for any JuMP extension storing data in the `ext` field.

!!! warning
    Do not engage in type piracy by implementing this method for types of `data`
    that you did not define! JuMP extensions should store types that they
    define in `model.ext`, rather than regular Julia types.
"""
function copy_extension_data(data, ::AbstractModel, ::AbstractModel)
    @warn(
        "Model contains extension data of type $(typeof(data)) that we do " *
        "not know how to copy.\n\nIf you are using a JuMP extension and you " *
        "did not add data to the `model.ext` dictionary. Please open an " *
        "issue on the GitHub repository of the JuMP extension and tell them " *
        "to implement `JuMP.copy_extension_data`.\n\nIf you added things to " *
        "`model.ext`, they have not been copied.",
    )
    return missing
end

"""
    GenericReferenceMap{T}

Mapping between variable and constraint reference of a model and its copy. The
reference of the copied model can be obtained by indexing the map with the
reference of the corresponding reference of the original model.
"""
struct GenericReferenceMap{T}
    model::GenericModel{T}
    index_map::MOIU.IndexMap
end

const ReferenceMap = GenericReferenceMap{Float64}

function Base.getindex(map::GenericReferenceMap, vref::GenericVariableRef)
    return GenericVariableRef(map.model, map.index_map[index(vref)])
end

function Base.getindex(map::GenericReferenceMap, cref::ConstraintRef)
    return ConstraintRef(map.model, map.index_map[index(cref)], cref.shape)
end

function Base.getindex(map::GenericReferenceMap, expr::GenericAffExpr)
    result = zero(expr)
    for (coef, var) in linear_terms(expr)
        add_to_expression!(result, coef, map[var])
    end
    result.constant = expr.constant
    return result
end

function Base.getindex(map::GenericReferenceMap, expr::GenericQuadExpr)
    aff = map[expr.aff]
    terms = [
        UnorderedPair(map[key.a], map[key.b]) => val for
        (key, val) in expr.terms
    ]
    return GenericQuadExpr(aff, terms)
end

function Base.getindex(map::GenericReferenceMap, val::AbstractArray)
    return getindex.(map, val)
end

Base.broadcastable(reference_map::GenericReferenceMap) = Ref(reference_map)

# Return a Boolean if the filtering function (1st argument) indicates that the whole value should
# be copied over.
_should_copy_complete_object(_, _) = true
function _should_copy_complete_object(
    filter_constraints::Function,
    value::ConstraintRef,
)
    return filter_constraints(value)
end
function _should_copy_complete_object(
    filter_constraints::Function,
    value::AbstractArray{T},
) where {T<:ConstraintRef}
    return all([filter_constraints(value[i]) for i in eachindex(value)])
end # all(filter_constraints.(value))

"""
    copy_model(model::GenericModel; filter_constraints::Union{Nothing, Function}=nothing)

Return a copy of the model `model` and a [`GenericReferenceMap`](@ref) that can
be used to obtain the variable and constraint reference of the new model
corresponding to a given `model`'s reference. A
[`Base.copy(::AbstractModel)`](@ref) method has also been implemented, it is
similar to `copy_model` but does not return the reference map.

If the `filter_constraints` argument is given, only the constraints for which
this function returns `true` will be copied. This function is given a
constraint reference as argument.

## Note

Model copy is not supported in `DIRECT` mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`GenericModel`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Example

In the following example, a model `model` is constructed with a variable `x` and
a constraint `cref`. It is then copied into a model `new_model` with the new
references assigned to `x_new` and `cref_new`.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, cref, x == 2)
cref : x = 2.0

julia> new_model, reference_map = copy_model(model);

julia> x_new = reference_map[x]
x

julia> cref_new = reference_map[cref]
cref : x = 2.0
```
"""
function copy_model(
    model::GenericModel{T};
    filter_constraints::Union{Nothing,Function} = nothing,
) where {T}
    if mode(model) == DIRECT
        error(
            "Cannot copy a model in `DIRECT` mode. Use the `Model` ",
            "constructor instead of the `direct_model` constructor to be ",
            "able to copy the constructed model.",
        )
    end
    new_model = GenericModel{T}()

    # At JuMP's level, filter_constraints should work with JuMP.ConstraintRef,
    # whereas MOI.copy_to's filter_constraints works with MOI.ConstraintIndex.
    function moi_filter(cref::MOI.ConstraintIndex)
        return filter_constraints(constraint_ref_with_index(model, cref))
    end
    moi_filter(::Any) = true

    index_map = if filter_constraints === nothing
        MOI.copy_to(backend(new_model), backend(model))
    else
        filtered_src = MOI.Utilities.ModelFilter(moi_filter, backend(model))
        MOI.copy_to(backend(new_model), filtered_src)
    end

    new_model.optimize_hook = model.optimize_hook

    # TODO copy NLP data
    if nonlinear_model(model) !== nothing
        error(
            "copy is not supported yet for models with nonlinear constraints",
            " and/or nonlinear objective function",
        )
    end

    reference_map = GenericReferenceMap(new_model, index_map)

    for (name, value) in object_dictionary(model)
        if _should_copy_complete_object(filter_constraints, value)
            try
                new_model[name] = getindex(reference_map, value)
            catch err
                if err isa MethodError
                    @warn(
                        "Skipping the copy of object `:$(name)` due to " *
                        "unsupported type $(typeof(value)). Please open a " *
                        "GitHub issue at https://github.com/jump-dev/JuMP.jl " *
                        "with this message.",
                    )
                else
                    rethrow(err)
                end
            end
        end
    end

    for (key, data) in model.ext
        new_model.ext[key] = copy_extension_data(data, new_model, model)
    end

    return new_model, reference_map
end

"""
    copy(model::AbstractModel)

Return a copy of the model `model`. It is similar to [`copy_model`](@ref)
except that it does not return the mapping between the references of `model`
and its copy.

## Note

Model copy is not supported in `DIRECT` mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`GenericModel`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Example

In the following example, a model `model` is constructed with a variable `x` and
a constraint `cref`. It is then copied into a model `new_model` with the new
references assigned to `x_new` and `cref_new`.

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> @constraint(model, cref, x == 2)
cref : x = 2.0

julia> new_model = copy(model);

julia> x_new = model[:x]
x

julia> cref_new = model[:cref]
cref : x = 2.0
```
"""
function Base.copy(model::AbstractModel)
    new_model, _ = copy_model(model)
    return new_model
end

"""
    copy_conflict(model::GenericModel)

Return a copy of the current conflict for the model `model` and a
[`GenericReferenceMap`](@ref) that can be used to obtain the variable and
constraint reference of the new model corresponding to a given `model`'s
reference.

This is a convenience function that provides a filtering function for
[`copy_model`](@ref).

## Note

Model copy is not supported in `DIRECT` mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`GenericModel`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Example

In the following example, a model `model` is constructed with a variable `x` and
two constraints `c1` and `c2`. This model has no solution, as the two
constraints are mutually exclusive. The solver is asked to compute a conflict
with [`compute_conflict!`](@ref). The parts of `model` participating in the
conflict are then copied into a model `iis_model`.

```julia
julia> using JuMP

julia> import Gurobi

julia> model = Model(Gurobi.Optimizer);

julia> set_silent(model)

julia> @variable(model, x >= 0)
x

julia> @constraint(model, c1, x >= 2)
c1 : x ≥ 2.0

julia> @constraint(model, c2, x <= 1)
c2 : x ≤ 1.0

julia> optimize!(model)

julia> compute_conflict!(model)

julia> if get_attribute(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
           iis_model, reference_map = copy_conflict(model)
           print(iis_model)
       end
Feasibility
Subject to
 c1 : x ≥ 2.0
 c2 : x ≤ 1.0
```
"""
function copy_conflict(model::GenericModel)
    filter_constraints =
        (cref) ->
            MOI.get(model, MOI.ConstraintConflictStatus(), cref) !=
            MOI.NOT_IN_CONFLICT
    new_model, reference_map =
        copy_model(model; filter_constraints = filter_constraints)
    return new_model, reference_map
end

# Calling `deepcopy` over a JuMP model is not supported, nor planned to be
# supported, because it would involve making a deep copy of the underlying
# solver (behind a C pointer).
function Base.deepcopy(::GenericModel)
    return error(
        "`JuMP.Model` does not support `deepcopy` as the reference to the underlying solver cannot be deep copied, use `copy` instead.",
    )
end

function MOI.copy_to(dest::MOI.ModelLike, src::GenericModel)
    if nonlinear_model(src) !== nothing
        # Re-set the NLP block in-case things have changed since last
        # solve.
        evaluator = MOI.Nonlinear.Evaluator(
            nonlinear_model(src),
            MOI.Nonlinear.SparseReverseMode(),
            index.(all_variables(src)),
        )
        MOI.set(src, MOI.NLPBlock(), MOI.NLPBlockData(evaluator))
    end
    return MOI.copy_to(dest, backend(src))
end

function MOI.copy_to(dest::GenericModel, src::MOI.ModelLike)
    index_map = MOI.copy_to(backend(dest), src)
    if MOI.NLPBlock() in MOI.get(src, MOI.ListOfModelAttributesSet())
        block = MOI.get(src, MOI.NLPBlock())
        dest.nlp_model = _nlp_model_from_nlpblock(block, block.evaluator)
    end
    return index_map
end

_lift_variable_from_expression(expr) = expr

function _lift_variable_from_expression(expr::Expr)
    if isexpr(expr, :ref, 2) && expr.args[1] == :x
        return expr.args[2]
    end
    for i in 1:length(expr.args)
        expr.args[i] = _lift_variable_from_expression(expr.args[i])
    end
    return expr
end

function _set_from_nlpboundspair(bound::MOI.NLPBoundsPair)
    if bound.lower == bound.upper
        return MOI.EqualTo(bound.lower)
    elseif isfinite(bound.lower) && !isfinite(bound.upper)
        return MOI.GreaterThan(bound.lower)
    elseif isfinite(bound.upper) && !isfinite(bound.lower)
        return MOI.LessThan(bound.upper)
    else
        return MOI.Interval(bound.lower, bound.upper)
    end
end

function _nlp_model_from_nlpblock(block::MOI.NLPBlockData, evaluator)
    if !(:ExprGraph in MOI.features_available(evaluator))
        error(
            "Unable to copy model because the nonlinear evaluator doesn't " *
            "support :ExprGraph.",
        )
    end
    MOI.initialize(evaluator, [:ExprGraph])
    model = MOI.Nonlinear.Model()
    for (i, bound) in enumerate(block.constraint_bounds)
        expr = MOI.constraint_expr(evaluator, i)
        f = isexpr(expr, :comparison) ? expr.args[3] : expr.args[2]
        MOI.Nonlinear.add_constraint(
            model,
            _lift_variable_from_expression(f),
            _set_from_nlpboundspair(bound),
        )
    end
    if block.has_objective
        expr = _lift_variable_from_expression(MOI.objective_expr(evaluator))
        MOI.Nonlinear.set_objective(model, expr)
    end
    return model
end
