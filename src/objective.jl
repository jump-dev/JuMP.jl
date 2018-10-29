"""
    objective_bound(model::Model)

Return the best known bound on the optimal objective value after a call to
`optimize!model)`.
"""
objective_bound(model::Model) = MOI.get(model, MOI.ObjectiveBound())

"""
    objective_value(model::Model)

Return the objective value after a call to `optimize!model)`.
"""
objective_value(model::Model) = MOI.get(model, MOI.ObjectiveValue())

"""
    objective_sense(model::Model)::MathOptInterface.OptimizationSense

Return the objective sense.
"""
function objective_sense(model::Model)
    return MOI.get(model, MOI.ObjectiveSense())
end

"""
    set_objective_sense(model::Model, sense::MathOptInterface.OptimizationSense)

Sets the objective sense of the model to the given sense. See
[`set_objective_function`](@ref) to set the objective function. These are
low-level functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_sense(model::Model, sense::MOI.OptimizationSense)
    MOI.set(model, MOI.ObjectiveSense(), sense)
end

"""
    set_objective_function(model::Model,
                           func::Union{AbstractJuMPScalar,
                                       MathOptInterface.AbstractScalarFunction})

Sets the objective function of the model to the given function. See
[`set_objective_sense`](@ref) to set the objective sense. These are low-level
functions; the recommended way to set the objective is with the
[`@objective`](@ref) macro.
"""
function set_objective_function end

function set_objective_function(model::Model, func::MOI.AbstractScalarFunction)
    attr = MOI.ObjectiveFunction{typeof(func)}()
    if !MOI.supports(model.moi_backend, attr)
        error("The solver does not support an objective function of type ",
              typeof(func), ".")
    end
    MOI.set(model, attr, func)
    # Keeping the explicit `return` is helpful for type inference because we
    # don't know what `MOI.set` will return.
    return
end

function set_objective_function(model::Model, func::AbstractJuMPScalar)
    set_objective_function(model, moi_function(func))
end

function set_objective(model::Model, sense::MOI.OptimizationSense,
                       func::AbstractJuMPScalar)
    set_objective_sense(model, sense)
    set_objective_function(model, func)
end
