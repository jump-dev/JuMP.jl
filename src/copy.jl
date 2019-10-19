#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    copy_extension_data(data, new_model::AbstractModel, model::AbstractModel)

Return a copy of the extension data `data` of the model `model` to the extension
data of the new model `new_model`. A method should be added for any JuMP
extension storing data in the `ext` field.
"""
function copy_extension_data end

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
function Base.getindex(reference_map::ReferenceMap, vref::VariableRef)
    return VariableRef(reference_map.model,
                       reference_map.index_map[index(vref)])
end
function Base.getindex(reference_map::ReferenceMap, cref::ConstraintRef)
    return ConstraintRef(reference_map.model,
                         reference_map.index_map[index(cref)],
                         cref.shape)
end
Base.broadcastable(reference_map::ReferenceMap) = Ref(reference_map)


"""
    copy_model(model::Model)

Return a copy of the model `model` and a [`ReferenceMap`](@ref) that can be used
to obtain the variable and constraint reference of the new model corresponding
to a given `model`'s reference. A [`Base.copy(::AbstractModel)`](@ref) method
has also been implemented, it is similar to `copy_model` but does not return
the reference map.

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
function copy_model(model::Model)
    if mode(model) == DIRECT
        error("Cannot copy a model in `DIRECT` mode. Use the `Model` ",
              "constructor instead of the `direct_model` constructor to be ",
              "able to copy the constructed model.")
    end
    caching_mode = backend(model).mode
    new_model = Model(caching_mode = caching_mode)

    # Copy the MOI backend, note that variable and constraint indices may have
    # changed, the `index_map` gives the map between the indices of
    # `backend(model` and the indices of `backend(new_model)`.
    index_map = MOI.copy_to(backend(new_model), backend(model),
                            copy_names = true)

    new_model.optimize_hook = model.optimize_hook

    # TODO copy NLP data
    if model.nlp_data !== nothing
        error("copy is not supported yet for models with nonlinear constraints",
              " and/or nonlinear objective function")
    end

    reference_map = ReferenceMap(new_model, index_map)

    for (name, value) in object_dictionary(model)
        new_model.obj_dict[name] = getindex.(reference_map, value)
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

# Calling `deepcopy` over a JuMP model is not supported, nor planned to be
# supported, because it would involve making a deep copy of the underlying
# solver (behind a C pointer).
function Base.deepcopy(::Model)
    error("`JuMP.Model` does not support `deepcopy` as the reference to the underlying solver cannot be deep copied, use `copy` instead.")
end
