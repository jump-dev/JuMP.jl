#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

module JuMPUnitfulExt

import JuMP
import MathOptInterface as MOI
import Unitful

# Similarly to complex numbers, the units are stored in the coefficient type of
# the expressions. A variable with units is therefore represented by a
# `GenericAffExpr` with a single term whose coefficient is `1` in the unit of the
# variable.

function JuMP._complex_convert_type(
    ::Type{T},
    ::Type{<:Unitful.Quantity{S,D,U}},
) where {T,S,D,U}
    return Unitful.Quantity{JuMP._complex_convert_type(T, S),D,U}
end

function JuMP._complex_convert(::Type{T}, x::Unitful.Quantity) where {T}
    return convert(JuMP._complex_convert_type(T, typeof(x)), x)
end

const _UnitfulExpr = Union{
    JuMP.GenericAffExpr{<:Unitful.Quantity},
    JuMP.GenericQuadExpr{<:Unitful.Quantity},
}

_unit(::JuMP.GenericAffExpr{C}) where {C} = Unitful.unit(C)

_unit(::JuMP.GenericQuadExpr{C}) where {C} = Unitful.unit(C)

function _ustrip(unit::Unitful.Units, f::_UnitfulExpr)
    return JuMP.map_coefficients(Base.Fix1(Unitful.ustrip, unit), f)
end

_ustrip(f::_UnitfulExpr) = _ustrip(_unit(f), f)

# Bounds without units are interpreted in the unit of the variable so that,
# for example, `@variable(model, x >= 0, u"m")` is supported.
function _ustrip_bound(unit::Unitful.Units, x::Unitful.Quantity)
    return Unitful.ustrip(unit, x)
end

_ustrip_bound(::Unitful.Units, x) = x

struct _UnitfulVariable{U<:Unitful.Units} <: JuMP.AbstractVariable
    variable::JuMP.ScalarVariable
    unit::U
end

function JuMP.build_variable(
    ::Function,
    info::JuMP.VariableInfo,
    unit::Unitful.Units,
)
    variable = JuMP.ScalarVariable(info)
    return _UnitfulVariable(
        JuMP._mapinfo(Base.Fix1(_ustrip_bound, unit), variable),
        unit,
    )
end

function JuMP.add_variable(
    model::JuMP.GenericModel{T},
    v::_UnitfulVariable,
    name::String = "",
) where {T}
    x = JuMP.add_variable(model, v.variable, name)
    return (one(T) * v.unit) * x
end

function JuMP.build_constraint(
    error_fn::Function,
    f::_UnitfulExpr,
    set::MOI.AbstractScalarSet,
)
    return JuMP.build_constraint(error_fn, f, set, _unit(f))
end

function JuMP.build_constraint(
    error_fn::Function,
    f::_UnitfulExpr,
    set::MOI.AbstractScalarSet,
    unit::Unitful.Units,
)
    return JuMP.build_constraint(error_fn, _ustrip(unit, f), set)
end

function JuMP.build_constraint(
    error_fn::Function,
    f::_UnitfulExpr,
    lb::Unitful.Quantity,
    ub::Unitful.Quantity,
)
    return JuMP.build_constraint(error_fn, f, lb, ub, _unit(f))
end

function JuMP.build_constraint(
    error_fn::Function,
    f::_UnitfulExpr,
    lb::Unitful.Quantity,
    ub::Unitful.Quantity,
    unit::Unitful.Units,
)
    return JuMP.build_constraint(
        error_fn,
        _ustrip(unit, f),
        Unitful.ustrip(unit, lb),
        Unitful.ustrip(unit, ub),
    )
end

function JuMP.set_objective_function(model::JuMP.GenericModel, f::_UnitfulExpr)
    return JuMP.set_objective_function(model, _ustrip(f))
end

function JuMP.function_string(mode::MIME, f::_UnitfulExpr)
    return string("(", JuMP.function_string(mode, _ustrip(f)), ") ", _unit(f))
end

end
