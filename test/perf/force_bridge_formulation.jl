using JuMP
import Clp
import GLPK

### No bridging

function foo(optimizer, force_bridge_formulation)
    model = Model(optimizer; bridge_constraints = force_bridge_formulation)
    set_silent(model)
    @variable(model, x >= 0)
    @constraint(model, 2x + 1 <= 1)
    @objective(model, Max, 1.0 * x)
    return optimize!(model)
end

### ObjectiveFunction bridging

function foo_obj(optimizer)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x >= 0)
    @constraint(model, 2x + 1 <= 1)
    @objective(model, Max, x)
    return optimize!(model)
end

### variable bridging

function foo_var(optimizer)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:1] in MOI.Nonnegatives(1))
    @constraint(model, 2x[1] + 1 <= 1)
    @objective(model, Max, 1.0 * x[1])
    return optimize!(model)
end

### constraint bridging

function foo_con(optimizer)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @constraint(model, x in MOI.Nonpositives(2))
    @objective(model, Max, 1.0 * x[1])
    return optimize!(model)
end

const optimizer = ARGS[1] == "clp" ? Clp.Optimizer : GLPK.Optimizer

if ARGS[2] == "--bridge"
    @eval @time foo(optimizer, true)
    @eval @time foo(optimizer, true)
elseif ARGS[2] == "--no-bridge"
    @eval @time foo(optimizer, false)
    @eval @time foo(optimizer, false)
elseif ARGS[2] == "--var"
    @eval @time foo_var(optimizer)
    @eval @time foo_var(optimizer)
elseif ARGS[2] == "--con"
    @eval @time foo_con(optimizer)
    @eval @time foo_con(optimizer)
else
    @assert ARGS[2] == "--obj"
    @eval @time foo_obj(optimizer)
    @eval @time foo_obj(optimizer)
end
