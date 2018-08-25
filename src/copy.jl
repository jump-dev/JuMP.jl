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
    copy_variablewise_constraints(dest::Dict{MOIVAR,
                                             MOICON{MOI.SingleVariable, S}},
                                  src::Dict{MOIVAR,
                                            MOICON{MOI.SingleVariable, S}},
                                  index_map) where S

Copy the variablewise constraint indices of `src` into `dest` mapping variable
and constraint indices using `index_map`.
"""
function copy_variablewise_constraints(dest::Dict{MOIVAR,
                                                  MOICON{MOI.SingleVariable, S}},
                                       src::Dict{MOIVAR,
                                                 MOICON{MOI.SingleVariable, S}},
                                       index_map) where S
    for (variable_index, constraint_index) in src
        dest[index_map[variable_index]] = index_map[constraint_index]
    end
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
function Base.getindex(reference_map::ReferenceMap, vref::VariableRef)
    return VariableRef(reference_map.model,
                       reference_map.index_map[index(vref)])
end
function Base.getindex(reference_map::ReferenceMap, cref::ConstraintRef)
    return ConstraintRef(reference_map.model,
                         reference_map.index_map[index(cref)],
                         cref.shape)
end

"""
    copy(model::Model)

Return a copy of the model `model` and a [`ReferenceMap`](@ref) that can be used
to obtain the variable and constraint reference of the new model corresponding
to a given `model`'s reference.

## Note

Model copy is not supported in Direct mode, i.e. when a model is constructed
using the [`direct_model`](@ref) constructor instead of the [`Model`](@ref)
constructor.

## Examples

In the following example, a model `model` is constructed with a variable `x` and
a constraint `cref`. It is then copied into a model `new_model` with the new
references assigned to `x_new` and `cref_new`.
```julia
model = Model()
@variable(model, x)
@constraint(model, cref, x == 2)

new_model, reference_map = JuMP.copy(model)
x_new = reference_map[x]
cref_new = reference_map[cref]
```
"""
function copy(model::Model)
    if mode(model) == Direct
        error("Cannot copy a model in Direct mode. Use the `Model` constructor",
              " instead of the `direct_model` constructor to be able to copy",
              " the constructed model.")
    end
    caching_mode = caching_optimizer(model).mode
    # TODO add bridges added to the bridge optimizer that are not part of the
    #      fullbridgeoptimizer
    bridge_constraints = model.moibackend isa MOI.Bridges.LazyBridgeOptimizer{<:MOIU.CachingOptimizer}
    new_model = Model(caching_mode = caching_mode,
                      bridge_constraints = bridge_constraints)

    # Copy the MOI backend, note that variable and constraint indices may have
    # changed, the `index_map` gives the map between the indices of
    # `model.moibackend` and the indices of `new_model.moibackend`.
    index_map = MOI.copy!(new_model.moibackend, model.moibackend,
                          copynames = true)
    # TODO copynames is needed because of https://github.com/JuliaOpt/MathOptInterface.jl/issues/494
    #      we can remove it when this is fixed and released

    copy_variablewise_constraints(new_model.variabletolowerbound,
                                  model.variabletolowerbound, index_map)
    copy_variablewise_constraints(new_model.variabletoupperbound,
                                  model.variabletoupperbound, index_map)
    copy_variablewise_constraints(new_model.variabletofix,
                                  model.variabletofix, index_map)
    copy_variablewise_constraints(new_model.variabletointegrality,
                                  model.variabletointegrality, index_map)
    copy_variablewise_constraints(new_model.variabletozeroone,
                                  model.variabletozeroone, index_map)

    new_model.optimizehook = model.optimizehook

    # TODO copy NLP data
    @assert model.nlpdata === nothing

    # TODO copy objdict

    for (key, data) in model.ext
        new_model.ext[key] = copy_extension_data(data, new_model, model)
    end

    return new_model, ReferenceMap(new_model, index_map)
end
