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
    set_callbacks(model::Model; lazy = nothing, heuristic = nothing)

Set an (optional) lazy and heuristic callback. This can be used only when
`JuMP.mode(mode) == DIRECT`.

Each callback must be a function that takes one argument, a solver-specific type
`cb_data`.

### Example

    set_callbacks(model;
        lazy = (cb_data) -> println("Called from lazy callback.")
    )
"""
function set_callbacks(model::Model; lazy = nothing, heuristic = nothing)
    if JuMP.mode(model) != JuMP.DIRECT
        error("You must use a solver in DIRECT mode to use callbacks in JuMP.")
    end
    MOI.set(JuMP.backend(model), MOI.Callbacks(lazy=lazy, heuristic=heuristic))
end

"""
    @lazy_constraint(model, cb_data, expr)

Add a lazy constraint (where the constraint is given by `expr`) to `model`.

This can be called *only* from a lazy callback set by `JuMP.set_callbacks`.

### Examples

    @lazy_constraint(model, cb_data, 2x + y <+ 1)
"""
macro lazy_constraint(model, cb_data, expr)
    code = quote
        lazy_con = @build_constraint $expr
        MOI.add_lazy_constraint(
            JuMP.backend($model),
            $cb_data,
            JuMP.moi_function(lazy_con.func),
            lazy_con.set
        )

    end
    quote
        let
            $(esc(code))
        end
    end
end

"""
    add_heuristic_solution(model, cb_data, sol::Dict{JuMP.VariableRef, Float64})

Provide the heuristic solution given by the variable-value mapping of `sol` to
`model`.

This can be called only from a heuristic callback set by `set_callbacks`.
"""
function add_heuristic_solution(
        model::JuMP.Model, cb_data, sol::Dict{JuMP.VariableRef, Float64})
    MOI.add_heuristic_solution(JuMP.backend(model), cb_data, Dict(
        index(variable) => value for (variable, value) in sol))
end
