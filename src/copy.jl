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
    ReferenceMap

Mapping between variable and constraint reference of a model and its copy. The
reference of the copied model can be obtained by indexing the map with the
reference of the corresponding reference of the original model.
"""
struct ReferenceMap
    model::Model
    index_map::MOIU.IndexMap
end

function Base.getindex(map::ReferenceMap, vref::VariableRef)
    return VariableRef(map.model, map.index_map[index(vref)])
end

function Base.getindex(map::ReferenceMap, cref::ConstraintRef)
    return ConstraintRef(map.model, map.index_map[index(cref)], cref.shape)
end

function Base.getindex(map::ReferenceMap, expr::GenericAffExpr)
    result = zero(expr)
    for (coef, var) in linear_terms(expr)
        add_to_expression!(result, coef, map[var])
    end
    result.constant = expr.constant
    return result
end

function Base.getindex(map::ReferenceMap, expr::GenericQuadExpr)
    aff = map[expr.aff]
    terms = [
        UnorderedPair(map[key.a], map[key.b]) => val for
        (key, val) in expr.terms
    ]
    return GenericQuadExpr(aff, terms)
end

Base.getindex(map::ReferenceMap, val::AbstractArray) = getindex.(map, val)

Base.broadcastable(reference_map::ReferenceMap) = Ref(reference_map)

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
    copy_model(model::Model; filter_constraints::Union{Nothing, Function}=nothing)

Return a copy of the model `model` and a [`ReferenceMap`](@ref) that can be used
to obtain the variable and constraint reference of the new model corresponding
to a given `model`'s reference. A [`Base.copy(::AbstractModel)`](@ref) method
has also been implemented, it is similar to `copy_model` but does not return
the reference map.

If the `filter_constraints` argument is given, only the constraints for which
this function returns `true` will be copied. This function is given a
constraint reference as argument.

## Note

Model copy is not supported in `DIRECT` mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`Model`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Examples

In the following example, a model `model` is constructed with a variable `x` and
a constraint `cref`. It is then copied into a model `new_model` with the new
references assigned to `x_new` and `cref_new`.
```julia
model = Model()
@variable(model, x)
@constraint(model, cref, x == 2)

new_model, reference_map = copy_model(model)
x_new = reference_map[x]
cref_new = reference_map[cref]
```
"""
function copy_model(
    model::Model;
    filter_constraints::Union{Nothing,Function} = nothing,
)
    if mode(model) == DIRECT
        error(
            "Cannot copy a model in `DIRECT` mode. Use the `Model` ",
            "constructor instead of the `direct_model` constructor to be ",
            "able to copy the constructed model.",
        )
    end
    new_model = Model()

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
    if model.nlp_data !== nothing
        error(
            "copy is not supported yet for models with nonlinear constraints",
            " and/or nonlinear objective function",
        )
    end

    reference_map = ReferenceMap(new_model, index_map)

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
using the [`direct_model`](@ref) constructor instead of the [`Model`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Examples

In the following example, a model `model` is constructed with a variable `x` and
a constraint `cref`. It is then copied into a model `new_model` with the new
references assigned to `x_new` and `cref_new`.
```julia
model = Model()
@variable(model, x)
@constraint(model, cref, x == 2)

new_model = copy(model)
x_new = model[:x]
cref_new = model[:cref]
```
"""
function Base.copy(model::AbstractModel)
    new_model, _ = copy_model(model)
    return new_model
end

"""
    copy_conflict(model::Model)

Return a copy of the current conflict for the model `model` and a
[`ReferenceMap`](@ref) that can be used to obtain the variable and constraint
reference of the new model corresponding to a given `model`'s reference.

This is a convenience function that provides a filtering function for
[`copy_model`](@ref).

## Note

Model copy is not supported in `DIRECT` mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`Model`](@ref)
constructor. Moreover, independently on whether an optimizer was provided at
model construction, the new model will have no optimizer, i.e., an optimizer
will have to be provided to the new model in the [`optimize!`](@ref) call.

## Examples

In the following example, a model `model` is constructed with a variable `x` and
two constraints `cref` and `cref2`. This model has no solution, as the two
constraints are mutually exclusive. The solver is asked to compute a conflict
with [`compute_conflict!`](@ref). The parts of `model` participating in the
conflict are then copied into a model `new_model`.
```julia
model = Model() # You must use a solver that supports conflict refining/IIS
# computation, like CPLEX or Gurobi
@variable(model, x)
@constraint(model, cref, x >= 2)
@constraint(model, cref2, x <= 1)

compute_conflict!(model)
if MOI.get(model, MOI.ConflictStatus()) != MOI.CONFLICT_FOUND
    error("No conflict could be found for an infeasible model.")
end

new_model, reference_map = copy_conflict(model)
```
"""
function copy_conflict(model::Model)
    filter_constraints =
        (cref) ->
            MOI.get(model, MOI.ConstraintConflictStatus(), cref) !=
            MOI.NOT_IN_CONFLICT
    new_model, reference_map =
        copy_model(model, filter_constraints = filter_constraints)
    return new_model, reference_map
end

# Calling `deepcopy` over a JuMP model is not supported, nor planned to be
# supported, because it would involve making a deep copy of the underlying
# solver (behind a C pointer).
function Base.deepcopy(::Model)
    return error(
        "`JuMP.Model` does not support `deepcopy` as the reference to the underlying solver cannot be deep copied, use `copy` instead.",
    )
end

function MOI.copy_to(dest::MOI.ModelLike, src::Model)
    if src.nlp_data !== nothing
        # Re-set the NLP block in-case things have changed since last
        # solve.
        block = MOI.NLPBlockData(src.nlp_data, index.(all_variables(src)))
        MOI.set(src, MOI.NLPBlock(), block)
    end
    return MOI.copy_to(dest, backend(src))
end

function MOI.copy_to(dest::Model, src::MOI.ModelLike)
    return MOI.copy_to(backend(dest), src)
end
