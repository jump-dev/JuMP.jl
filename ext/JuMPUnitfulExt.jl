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
# `GenericAffExpr` with a single term whose coefficient is the scale of the
# variable, for example, `1 m/s` for `u"m/s"` or `10 m/s` for `10u"m/s"`. The
# scale can be used to improve the scaling of the problem.

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

struct _UnitfulVariable{S<:Union{Unitful.Units,Unitful.Quantity}} <:
       JuMP.AbstractVariable
    variable::JuMP.ScalarVariable
    scale::S
end

function _strip_bound(::Function, scale, x::Unitful.Quantity)
    return Unitful.uconvert(Unitful.NoUnits, x / scale)
end

function _strip_bound(error_fn::Function, scale, x)
    return error_fn(
        """
        Unable to use `$x` as a bound or start value because it has no units.

        Use a quantity with the same dimension as `$scale`.
        """,
    )
end

function _strip_bound(error_fn::Function, scale, has_value::Bool, x)
    if has_value
        return _strip_bound(error_fn, scale, x)
    end
    return x
end

function JuMP.build_variable(
    error_fn::Function,
    info::JuMP.VariableInfo,
    scale::Union{Unitful.Units,Unitful.Quantity},
)
    new_info = JuMP.VariableInfo(
        info.has_lb,
        _strip_bound(error_fn, scale, info.has_lb, info.lower_bound),
        info.has_ub,
        _strip_bound(error_fn, scale, info.has_ub, info.upper_bound),
        info.has_fix,
        _strip_bound(error_fn, scale, info.has_fix, info.fixed_value),
        info.has_start,
        _strip_bound(error_fn, scale, info.has_start, info.start),
        info.binary,
        info.integer,
    )
    return _UnitfulVariable(JuMP.ScalarVariable(new_info), scale)
end

function JuMP.add_variable(
    model::JuMP.GenericModel{T},
    v::_UnitfulVariable,
    name::String = "",
) where {T}
    x = JuMP.add_variable(model, v.variable, name)
    return (one(T) * v.scale) * x
end

# The units are stripped from the constraint before it is passed to MOI. The
# unit is stored in the shape of the constraint so that `value` and `dual` can
# add the units back.
struct _UnitfulShape{U<:Unitful.Units} <: JuMP.AbstractShape
    unit::U
end

function JuMP.reshape_vector(
    x::Union{Number,JuMP.AbstractJuMPScalar},
    shape::_UnitfulShape,
)
    return x * (1 * shape.unit)
end

JuMP.vectorize(x, shape::_UnitfulShape) = Unitful.ustrip(shape.unit, x)

JuMP.vectorize(::Nothing, ::_UnitfulShape) = nothing

# The objective has no units so the dual has the inverse unit of the
# constraint.
JuMP.dual_shape(shape::_UnitfulShape) = _UnitfulShape(inv(shape.unit))

function JuMP.in_set_string(
    mode,
    con::JuMP.ScalarConstraint{<:Any,<:Any,<:_UnitfulShape},
)
    set_str = JuMP.in_set_string(mode, JuMP.moi_set(con))
    return string(set_str, " ", con.shape.unit)
end

function _with_unit(con::JuMP.ScalarConstraint, unit::Unitful.Units)
    return JuMP.ScalarConstraint(con.func, con.set, _UnitfulShape(unit))
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
    con = JuMP.build_constraint(error_fn, _ustrip(unit, f), set)
    return _with_unit(con, unit)
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
    con = JuMP.build_constraint(
        error_fn,
        _ustrip(unit, f),
        Unitful.ustrip(unit, lb),
        Unitful.ustrip(unit, ub),
    )
    return _with_unit(con, unit)
end

function JuMP.set_objective_function(model::JuMP.GenericModel, f::_UnitfulExpr)
    return JuMP.set_objective_function(model, _ustrip(f))
end

function JuMP.function_string(mode::MIME, f::_UnitfulExpr)
    return string("(", JuMP.function_string(mode, _ustrip(f)), ") ", _unit(f))
end

end
