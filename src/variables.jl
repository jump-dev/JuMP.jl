#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    AbstractVariable

Variable returned by [`build_variable`](@ref). It represents a variable that has
not been added yet to any model. It can be added to a given `model` with
[`add_variable`](@ref).
"""
abstract type AbstractVariable end

# Any fields can usually be either a number or an expression
mutable struct _VariableInfoExpr
    has_lb::Bool
    lower_bound::Any
    has_ub::Bool
    upper_bound::Any
    has_fix::Bool
    fixed_value::Any
    has_start::Bool
    start::Any
    binary::Any
    integer::Any
end

function _set_lower_bound_or_error(
    _error::Function,
    info::_VariableInfoExpr,
    lower,
)
    info.has_lb && _error("Cannot specify variable lower_bound twice")
    info.has_lb = true
    return info.lower_bound = lower
end
function _set_upper_bound_or_error(
    _error::Function,
    info::_VariableInfoExpr,
    upper,
)
    info.has_ub && _error("Cannot specify variable upper_bound twice")
    info.has_ub = true
    return info.upper_bound = upper
end
function _fix_or_error(_error::Function, info::_VariableInfoExpr, value)
    info.has_fix && _error("Cannot specify variable fixed value twice")
    info.has_fix = true
    return info.fixed_value = value
end
function _set_binary_or_error(_error::Function, info::_VariableInfoExpr)
    info.binary === false ||
        _error("'Bin' and 'binary' keyword argument cannot both be specified.")
    return info.binary = true
end
function _set_integer_or_error(_error::Function, info::_VariableInfoExpr)
    info.integer === false ||
        _error("'Int' and 'integer' keyword argument cannot both be specified.")
    return info.integer = true
end

function _is_info_keyword(kw::Expr)
    return kw.args[1] in [:lower_bound, :upper_bound, :start, :binary, :integer]
end
# :(start = 0)     -> (:start, 0)
# :(start = i + 1) -> (:start, :($(Expr(:escape, :(i + 1)))))
function _keywordify(kw::Expr)
    return (kw.args[1], _esc_non_constant(kw.args[2]))
end
function _VariableInfoExpr(;
    lower_bound = NaN,
    upper_bound = NaN,
    start = NaN,
    binary = false,
    integer = false,
)
    # isnan(::Expr) is not defined so we need to do !== NaN
    return _VariableInfoExpr(
        lower_bound !== NaN,
        lower_bound,
        upper_bound !== NaN,
        upper_bound,
        false,
        NaN,
        start !== NaN,
        start,
        binary,
        integer,
    )
end

# It isn't sufficient to use `isfinite` below, because some bounds are given as
# matrices. As a fallback, we define `_isfinite`, because overloading `isfinite`
# would be type piracy.
_isfinite(x::Number) = isfinite(x)
_isfinite(x) = true

"""
    VariableInfo{S,T,U,V}

A struct by JuMP internally when creating variables. This may also be used by
JuMP extensions to create new types of variables.

See also: [`ScalarVariable`](@ref).
"""
struct VariableInfo{S,T,U,V}
    has_lb::Bool
    lower_bound::S
    has_ub::Bool
    upper_bound::T
    has_fix::Bool
    fixed_value::U
    has_start::Bool
    start::V
    binary::Bool
    integer::Bool
    function VariableInfo(
        has_lb::Bool,
        lower_bound::S,
        has_ub::Bool,
        upper_bound::T,
        has_fix::Bool,
        fixed_value::U,
        has_start::Bool,
        start::V,
        binary::Bool,
        integer::Bool,
    ) where {S,T,U,V}
        if has_lb && !_isfinite(lower_bound)
            has_lb = false
            lower_bound = NaN
        end
        if has_ub && !_isfinite(upper_bound)
            has_ub = false
            upper_bound = NaN
        end
        if has_fix && !_isfinite(fixed_value)
            error("Unable to fix variable to $(fixed_value)")
        end
        return new{typeof(lower_bound),typeof(upper_bound),U,V}(
            has_lb,
            lower_bound,
            has_ub,
            upper_bound,
            has_fix,
            fixed_value,
            has_start,
            start,
            binary,
            integer,
        )
    end
end

function _constructor_expr(info::_VariableInfoExpr)
    return :(VariableInfo(
        $(info.has_lb),
        $(info.lower_bound),
        $(info.has_ub),
        $(info.upper_bound),
        $(info.has_fix),
        $(info.fixed_value),
        $(info.has_start),
        $(info.start),
        $(info.binary),
        $(info.integer),
    ))
end

"""
    ScalarVariable{S,T,U,V} <: AbstractVariable

A struct used when adding variables.

See also: [`add_variable`](@ref).
"""
struct ScalarVariable{S,T,U,V} <: AbstractVariable
    info::VariableInfo{S,T,U,V}
end

"""
    AbstractVariableRef

Variable returned by [`add_variable`](@ref). Affine (resp. quadratic) operations
with variables of type `V<:AbstractVariableRef` and coefficients of type `T`
    create a `GenericAffExpr{T,V}` (resp. `GenericQuadExpr{T,V}`).
"""
abstract type AbstractVariableRef <: AbstractJuMPScalar end

variable_ref_type(v::AbstractVariableRef) = typeof(v)
Base.conj(v::AbstractVariableRef) = v
Base.real(v::AbstractVariableRef) = v
Base.imag(v::AbstractVariableRef) = v
Base.abs2(v::AbstractVariableRef) = v^2

"""
    VariableRef <: AbstractVariableRef

Holds a reference to the model and the corresponding MOI.VariableIndex.
"""
struct VariableRef <: AbstractVariableRef
    model::Model
    index::MOI.VariableIndex
end

# `AbstractVariableRef` types must override the default `owner_model` if the field
#  name is not `model`.
"""
    owner_model(v::AbstractVariableRef)

Returns the model to which `v` belongs.

# Example
```jldoctest; setup=:(using JuMP)
julia> model = Model();

julia> x = @variable(model)
_[1]

julia> owner_model(x) === model
true
```
"""
owner_model(v::AbstractVariableRef) = v.model

"""
    struct VariableNotOwned{V <: AbstractVariableRef} <: Exception
        variable::V
    end

The variable `variable` was used in a model different to
`owner_model(variable)`.
"""
struct VariableNotOwned{V<:AbstractVariableRef} <: Exception
    variable::V
end

"""
    check_belongs_to_model(func::AbstractJuMPScalar, model::AbstractModel)

Throw [`VariableNotOwned`](@ref) if the [`owner_model`](@ref) of one of the
variables of the function `func` is not `model`.

    check_belongs_to_model(constraint::AbstractConstraint, model::AbstractModel)

Throw [`VariableNotOwned`](@ref) if the [`owner_model`](@ref) of one of the
variables of the constraint `constraint` is not `model`.
"""
function check_belongs_to_model end

function check_belongs_to_model(v::AbstractVariableRef, model::AbstractModel)
    if owner_model(v) !== model
        throw(VariableNotOwned(v))
    end
end

Base.iszero(::VariableRef) = false
Base.copy(v::VariableRef) = VariableRef(v.model, v.index)
Base.broadcastable(v::VariableRef) = Ref(v)

"""
    coefficient(v1::VariableRef, v2::VariableRef)

Return `1.0` if `v1 == v2`, and `0.0` otherwise.

This is a fallback for other [`coefficient`](@ref) methods to simplify code in
which the expression may be a single variable.
"""
coefficient(v1::VariableRef, v2::VariableRef) = (v1 == v2 ? 1.0 : 0.0)
coefficient(v1::VariableRef, v2::VariableRef, v3::VariableRef) = 0.0

isequal_canonical(v::VariableRef, other::VariableRef) = isequal(v, other)

"""
    delete(model::Model, variable_ref::VariableRef)

Delete the variable associated with `variable_ref` from the model `model`.

Note that `delete` does not unregister the name from the model, so adding a new
variable of the same name will throw an error. Use [`unregister`](@ref) to
unregister the name after deletion as follows:
```julia
@variable(model, x)
delete(model, x)
unregister(model, :x)
```

See also: [`unregister`](@ref)
"""
function delete(model::Model, variable_ref::VariableRef)
    if model !== owner_model(variable_ref)
        error(
            "The variable reference you are trying to delete does not " *
            "belong to the model.",
        )
    end
    model.is_model_dirty = true
    return MOI.delete(backend(model), variable_ref.index)
end

"""
    delete(model::Model, variable_refs::Vector{VariableRef})

Delete the variables associated with `variable_refs` from the model `model`.
Solvers may implement methods for deleting multiple variables that are
more efficient than repeatedly calling the single variable delete method.

See also: [`unregister`](@ref)
"""
function delete(model::Model, variable_refs::Vector{VariableRef})
    if any(model !== owner_model(v) for v in variable_refs)
        error(
            "A variable reference you are trying to delete does not " *
            "belong to the model.",
        )
    end
    model.is_model_dirty = true
    MOI.delete(backend(model), index.(variable_refs))
    return
end

"""
    is_valid(model::Model, variable_ref::VariableRef)

Return `true` if `variable` refers to a valid variable in `model`.
"""
function is_valid(model::Model, variable_ref::VariableRef)
    return (
        model === owner_model(variable_ref) &&
        MOI.is_valid(backend(model), variable_ref.index)
    )
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/jump-dev/MathOptInterface.jl/issues/234#issuecomment-366868878
function Base.hash(v::VariableRef, h::UInt)
    return hash(objectid(owner_model(v)), hash(v.index.value, h))
end
function Base.isequal(v1::VariableRef, v2::VariableRef)
    return owner_model(v1) === owner_model(v2) && v1.index == v2.index
end

"""
    index(v::VariableRef)::MOI.VariableIndex

Return the index of the variable that corresponds to `v` in the MOI backend.
"""
index(v::VariableRef) = v.index

function VariableRef(m::Model)
    index = MOI.add_variable(backend(m))
    return VariableRef(m, index)
end

"""
    VariableRef(c::ConstraintRef)

Get the variable associated with a `ConstraintRef`, if `c` is a constraint on a
single variable.

## Examples

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0)
x

julia> c = LowerBoundRef(x)
x ≥ 0.0

julia> VariableRef(c) == x
true
```
"""
function VariableRef(
    c::ConstraintRef{<:AbstractModel,<:MOI.ConstraintIndex{MOI.VariableIndex}},
)
    return VariableRef(owner_model(c), MOI.VariableIndex(index(c).value))
end

# Name setter/getters
# These functions need to be implemented for all `AbstractVariableRef`s
"""
    name(v::VariableRef)::String

Get a variable's name attribute.
"""
function name(v::VariableRef)
    model = owner_model(v)
    if !MOI.supports(backend(model), MOI.VariableName(), MOI.VariableIndex)
        return ""
    end
    return MOI.get(model, MOI.VariableName(), v)::String
end

"""
    set_name(v::VariableRef, s::AbstractString)

Set a variable's name attribute.
"""
function set_name(v::VariableRef, s::String)
    return MOI.set(owner_model(v), MOI.VariableName(), v, s)
end

"""
    variable_by_name(model::AbstractModel,
                     name::String)::Union{AbstractVariableRef, Nothing}

Returns the reference of the variable with name attribute `name` or `Nothing` if
no variable has this name attribute. Throws an error if several variables have
`name` as their name attribute.

```jldoctest objective_function; setup = :(using JuMP), filter = r"Stacktrace:.*"s
julia> model = Model();

julia> @variable(model, x)
x

julia> variable_by_name(model, "x")
x

julia> @variable(model, base_name="x")
x

julia> variable_by_name(model, "x")
ERROR: Multiple variables have the name x.
Stacktrace:
 [1] error(::String) at ./error.jl:33
 [2] get(::MOIU.Model{Float64}, ::Type{MathOptInterface.VariableIndex}, ::String) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/model.jl:222
 [3] get at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/universalfallback.jl:201 [inlined]
 [4] get(::MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.AbstractOptimizer,MathOptInterface.Utilities.UniversalFallback{MOIU.Model{Float64}}}, ::Type{MathOptInterface.VariableIndex}, ::String) at /home/blegat/.julia/dev/MathOptInterface/src/Utilities/cachingoptimizer.jl:490
 [5] variable_by_name(::Model, ::String) at /home/blegat/.julia/dev/JuMP/src/variables.jl:268
 [6] top-level scope at none:0

julia> var = @variable(model, base_name="y")
y

julia> variable_by_name(model, "y")
y

julia> set_name(var, "z")

julia> variable_by_name(model, "y")

julia> variable_by_name(model, "z")
z

julia> @variable(model, u[1:2])
2-element Array{VariableRef,1}:
 u[1]
 u[2]

julia> variable_by_name(model, "u[2]")
u[2]
```
"""
function variable_by_name(model::Model, name::String)
    index = MOI.get(backend(model), MOI.VariableIndex, name)
    if index isa Nothing
        return nothing
    else
        return VariableRef(model, index)
    end
end

MOI.VariableIndex(v::VariableRef) = index(v)

moi_function(variable::AbstractVariableRef) = index(variable)

moi_function_type(::Type{<:AbstractVariableRef}) = MOI.VariableIndex

# Note: No validation is performed that the variables belong to the same model.
function MOI.VectorOfVariables(vars::Vector{VariableRef})
    return MOI.VectorOfVariables(index.(vars))
end
function moi_function(variables::Vector{<:AbstractVariableRef})
    return MOI.VectorOfVariables(variables)
end
function moi_function_type(::Type{<:Vector{<:AbstractVariableRef}})
    return MOI.VectorOfVariables
end
function jump_function_type(::Model, ::Type{MOI.VectorOfVariables})
    return Vector{VariableRef}
end
function jump_function(model::Model, variables::MOI.VectorOfVariables)
    return VariableRef[VariableRef(model, v) for v in variables.variables]
end

jump_function_type(::Model, ::Type{MOI.VariableIndex}) = VariableRef
function jump_function(model::Model, variable::MOI.VariableIndex)
    return VariableRef(model, variable)
end

## Bound setter/getters

# lower bounds

"""
    has_lower_bound(v::VariableRef)

Return `true` if `v` has a lower bound. If `true`, the lower bound can be
queried with [`lower_bound`](@ref).

See also [`LowerBoundRef`](@ref), [`lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function has_lower_bound(v::VariableRef)
    return _moi_has_lower_bound(backend(owner_model(v)), v)
end

# _moi_* methods allow us to work around the type instability of the backend of
# a model.
function _moi_has_lower_bound(moi_backend, v::VariableRef)
    return MOI.is_valid(moi_backend, _lower_bound_index(v))
end

function _lower_bound_index(v::VariableRef)
    return _MOICON{MOI.VariableIndex,MOI.GreaterThan{Float64}}(index(v).value)
end

"""
    set_lower_bound(v::VariableRef, lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower
bound constraint.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function set_lower_bound(v::VariableRef, lower::Number)
    if !isfinite(lower)
        error(
            "Unable to set lower bound to $(lower). To remove the bound, use " *
            "`delete_lower_bound`.",
        )
    end
    model = owner_model(v)
    model.is_model_dirty = true
    return _moi_set_lower_bound(backend(model), v, lower)
end

function _moi_set_lower_bound(moi_backend, v::VariableRef, lower::Number)
    new_set = MOI.GreaterThan(convert(Float64, lower))
    if _moi_has_lower_bound(moi_backend, v)
        cindex = _lower_bound_index(v)
        MOI.set(moi_backend, MOI.ConstraintSet(), cindex, new_set)
    else
        @assert !_moi_is_fixed(moi_backend, v)
        _moi_add_constraint(moi_backend, index(v), new_set)
    end
    return
end

"""
    LowerBoundRef(v::VariableRef)

Return a constraint reference to the lower bound constraint of `v`. Errors if
one does not exist.

See also [`has_lower_bound`](@ref), [`lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function LowerBoundRef(v::VariableRef)
    moi_lb = _MOICON{MOI.VariableIndex,MOI.GreaterThan{Float64}}
    return ConstraintRef{Model,moi_lb,ScalarShape}(
        owner_model(v),
        _lower_bound_index(v),
        ScalarShape(),
    )
end

"""
    delete_lower_bound(v::VariableRef)

Delete the lower bound constraint of a variable.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`lower_bound`](@ref), [`set_lower_bound`](@ref).
"""
function delete_lower_bound(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), LowerBoundRef(variable_ref))
    return
end

"""
    lower_bound(v::VariableRef)

Return the lower bound of a variable. Error if one does not exist.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function lower_bound(v::VariableRef)
    if !has_lower_bound(v)
        error("Variable $(v) does not have a lower bound.")
    end
    cset = MOI.get(
        owner_model(v),
        MOI.ConstraintSet(),
        LowerBoundRef(v),
    )::MOI.GreaterThan{Float64}
    return cset.lower
end

# upper bounds

"""
    has_upper_bound(v::VariableRef)

Return `true` if `v` has a upper bound. If `true`, the upper bound can be
queried with [`upper_bound`](@ref).

See also [`UpperBoundRef`](@ref), [`upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function has_upper_bound(v::VariableRef)
    return _moi_has_upper_bound(backend(owner_model(v)), v)
end

function _moi_has_upper_bound(moi_backend, v::VariableRef)
    return MOI.is_valid(moi_backend, _upper_bound_index(v))
end

function _upper_bound_index(v::VariableRef)
    return _MOICON{MOI.VariableIndex,MOI.LessThan{Float64}}(index(v).value)
end

"""
    set_upper_bound(v::VariableRef, upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound
constraint.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function set_upper_bound(v::VariableRef, upper::Number)
    if !isfinite(upper)
        error(
            "Unable to set upper bound to $(upper). To remove the bound, use " *
            "`delete_upper_bound`.",
        )
    end
    model = owner_model(v)
    model.is_model_dirty = true
    return _moi_set_upper_bound(backend(model), v, upper)
end

function _moi_set_upper_bound(moi_backend, v::VariableRef, upper::Number)
    new_set = MOI.LessThan(convert(Float64, upper))
    if _moi_has_upper_bound(moi_backend, v)
        cindex = _upper_bound_index(v)
        MOI.set(moi_backend, MOI.ConstraintSet(), cindex, new_set)
    else
        @assert !_moi_is_fixed(moi_backend, v)
        _moi_add_constraint(moi_backend, index(v), new_set)
    end
    return
end

"""
    UpperBoundRef(v::VariableRef)

Return a constraint reference to the upper bound constraint of `v`. Errors if
one does not exist.

See also [`has_upper_bound`](@ref), [`upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function UpperBoundRef(v::VariableRef)
    moi_ub = _MOICON{MOI.VariableIndex,MOI.LessThan{Float64}}
    return ConstraintRef{Model,moi_ub,ScalarShape}(
        owner_model(v),
        _upper_bound_index(v),
        ScalarShape(),
    )
end

"""
    delete_upper_bound(v::VariableRef)

Delete the upper bound constraint of a variable.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`upper_bound`](@ref), [`set_upper_bound`](@ref).
"""
function delete_upper_bound(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), UpperBoundRef(variable_ref))
    return
end

"""
    upper_bound(v::VariableRef)

Return the upper bound of a variable. Error if one does not exist.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function upper_bound(v::VariableRef)
    if !has_upper_bound(v)
        error("Variable $(v) does not have an upper bound.")
    end
    cset = MOI.get(
        owner_model(v),
        MOI.ConstraintSet(),
        UpperBoundRef(v),
    )::MOI.LessThan{Float64}
    return cset.upper
end

# fixed value

"""
    is_fixed(v::VariableRef)

Return `true` if `v` is a fixed variable. If `true`, the fixed value can be
queried with [`fix_value`](@ref).

See also [`FixRef`](@ref), [`fix_value`](@ref), [`fix`](@ref), [`unfix`](@ref).
"""
function is_fixed(v::VariableRef)
    return _moi_is_fixed(backend(owner_model(v)), v)
end

function _moi_is_fixed(moi_backend, v::VariableRef)
    return MOI.is_valid(moi_backend, _fix_index(v))
end

function _fix_index(v::VariableRef)
    return _MOICON{MOI.VariableIndex,MOI.EqualTo{Float64}}(index(v).value)
end

"""
    fix(v::VariableRef, value::Number; force::Bool = false)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise
create a new one.

If the variable already has variable bounds and `force=false`, calling `fix`
will throw an error. If `force=true`, existing variable bounds will be deleted,
and the fixing constraint will be added. Note a variable will have no bounds
after a call to [`unfix`](@ref).

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix_value`](@ref),
[`unfix`](@ref).
"""
function fix(variable::VariableRef, value::Number; force::Bool = false)
    if !isfinite(value)
        error("Unable to fix variable to $(value)")
    end
    model = owner_model(variable)
    model.is_model_dirty = true
    return _moi_fix(backend(model), variable, value, force)
end

function _moi_fix(
    moi_backend,
    variable::VariableRef,
    value::Number,
    force::Bool,
)
    new_set = MOI.EqualTo(convert(Float64, value))
    if _moi_is_fixed(moi_backend, variable)  # Update existing fixing constraint.
        c_index = _fix_index(variable)
        MOI.set(moi_backend, MOI.ConstraintSet(), c_index, new_set)
    else  # Add a new fixing constraint.
        if _moi_has_upper_bound(moi_backend, variable) ||
           _moi_has_lower_bound(moi_backend, variable)
            if !force
                error(
                    "Unable to fix $(variable) to $(value) because it has " *
                    "existing variable bounds. Consider calling " *
                    "`JuMP.fix(variable, value; force=true)` which will " *
                    "delete existing bounds before fixing the variable.",
                )
            end
            if _moi_has_upper_bound(moi_backend, variable)
                MOI.delete(moi_backend, _upper_bound_index(variable))
            end
            if _moi_has_lower_bound(moi_backend, variable)
                MOI.delete(moi_backend, _lower_bound_index(variable))
            end
        end
        _moi_add_constraint(moi_backend, index(variable), new_set)
    end
    return
end

"""
    unfix(v::VariableRef)

Delete the fixing constraint of a variable.

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix_value`](@ref),
[`fix`](@ref).
"""
function unfix(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), FixRef(variable_ref))
    return
end

"""
    fix_value(v::VariableRef)

Return the value to which a variable is fixed. Error if one does not exist.

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix`](@ref), [`unfix`](@ref).
"""
function fix_value(v::VariableRef)
    cset = MOI.get(
        owner_model(v),
        MOI.ConstraintSet(),
        FixRef(v),
    )::MOI.EqualTo{Float64}
    return cset.value
end

"""
    FixRef(v::VariableRef)

Return a constraint reference to the constraint fixing the value of `v`. Errors
if one does not exist.

See also [`is_fixed`](@ref), [`fix_value`](@ref), [`fix`](@ref),
[`unfix`](@ref).
"""
function FixRef(v::VariableRef)
    moi_fix = _MOICON{MOI.VariableIndex,MOI.EqualTo{Float64}}
    return ConstraintRef{Model,moi_fix,ScalarShape}(
        owner_model(v),
        _fix_index(v),
        ScalarShape(),
    )
end

"""
    is_integer(v::VariableRef)

Return `true` if `v` is constrained to be integer.

See also [`IntegerRef`](@ref), [`set_integer`](@ref), [`unset_integer`](@ref).
"""
function is_integer(v::VariableRef)
    return _moi_is_integer(backend(owner_model(v)), v)
end

function _moi_is_integer(moi_backend, v::VariableRef)
    return MOI.is_valid(moi_backend, _integer_index(v))
end

function _integer_index(v::VariableRef)
    return _MOICON{MOI.VariableIndex,MOI.Integer}(index(v).value)
end

"""
    set_integer(variable_ref::VariableRef)

Add an integrality constraint on the variable `variable_ref`.

See also [`IntegerRef`](@ref), [`is_integer`](@ref), [`unset_integer`](@ref).
"""
function set_integer(v::VariableRef)
    model = owner_model(v)
    model.is_model_dirty = true
    return _moi_set_integer(backend(model), v)
end

function _moi_set_integer(moi_backend, variable_ref::VariableRef)
    if _moi_is_integer(moi_backend, variable_ref)
        return
    elseif _moi_is_binary(moi_backend, variable_ref)
        error(
            "Cannot set the variable_ref $(variable_ref) to integer as it " *
            "is already binary.",
        )
    end
    _moi_add_constraint(moi_backend, index(variable_ref), MOI.Integer())
    return
end

"""
    unset_integer(variable_ref::VariableRef)

Remove the integrality constraint on the variable `variable_ref`.

See also [`IntegerRef`](@ref), [`is_integer`](@ref), [`set_integer`](@ref).
"""
function unset_integer(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), IntegerRef(variable_ref))
    return
end

"""
    IntegerRef(v::VariableRef)

Return a constraint reference to the constraint constraining `v` to be integer.
Errors if one does not exist.

See also [`is_integer`](@ref), [`set_integer`](@ref), [`unset_integer`](@ref).
"""
function IntegerRef(v::VariableRef)
    moi_int = _MOICON{MOI.VariableIndex,MOI.Integer}
    return ConstraintRef{Model,moi_int,ScalarShape}(
        owner_model(v),
        _integer_index(v),
        ScalarShape(),
    )
end

"""
    is_binary(v::VariableRef)

Return `true` if `v` is constrained to be binary.

See also [`BinaryRef`](@ref), [`set_binary`](@ref), [`unset_binary`](@ref).
"""
function is_binary(v::VariableRef)
    return _moi_is_binary(backend(owner_model(v)), v)
end

function _moi_is_binary(moi_backend, v::VariableRef)
    return MOI.is_valid(moi_backend, _binary_index(v))
end

function _binary_index(v::VariableRef)
    return _MOICON{MOI.VariableIndex,MOI.ZeroOne}(index(v).value)
end

"""
    set_binary(v::VariableRef)

Add a constraint on the variable `v` that it must take values in the set
``\\{0,1\\}``.

See also [`BinaryRef`](@ref), [`is_binary`](@ref), [`unset_binary`](@ref).
"""
function set_binary(v::VariableRef)
    model = owner_model(v)
    model.is_model_dirty = true
    return _moi_set_binary(backend(model), v)
end

function _moi_set_binary(moi_backend, variable_ref)
    if _moi_is_binary(moi_backend, variable_ref)
        return
    elseif _moi_is_integer(moi_backend, variable_ref)
        error(
            "Cannot set the variable_ref $(variable_ref) to binary as it " *
            "is already integer.",
        )
    end
    _moi_add_constraint(moi_backend, index(variable_ref), MOI.ZeroOne())
    return
end

"""
    unset_binary(variable_ref::VariableRef)

Remove the binary constraint on the variable `variable_ref`.

See also [`BinaryRef`](@ref), [`is_binary`](@ref), [`set_binary`](@ref).
"""
function unset_binary(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), BinaryRef(variable_ref))
    return
end

"""
    BinaryRef(v::VariableRef)

Return a constraint reference to the constraint constraining `v` to be binary.
Errors if one does not exist.

See also [`is_binary`](@ref), [`set_binary`](@ref), [`unset_binary`](@ref).
"""
function BinaryRef(v::VariableRef)
    moi_bin = _MOICON{MOI.VariableIndex,MOI.ZeroOne}
    return ConstraintRef{Model,moi_bin,ScalarShape}(
        owner_model(v),
        _binary_index(v),
        ScalarShape(),
    )
end

"""
    start_value(v::VariableRef)

Return the start value (MOI attribute `VariablePrimalStart`) of the variable
`v`.

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".

See also [`set_start_value`](@ref).
"""
function start_value(v::VariableRef)::Union{Nothing,Float64}
    return MOI.get(owner_model(v), MOI.VariablePrimalStart(), v)
end

"""
    set_start_value(variable::VariableRef, value::Union{Real,Nothing})

Set the start value (MOI attribute `VariablePrimalStart`) of the `variable` to
`value`.

Pass `nothing` to unset the start value.

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".

See also [`start_value`](@ref).
"""
function set_start_value(variable::VariableRef, value::Union{Nothing,Float64})
    MOI.set(owner_model(variable), MOI.VariablePrimalStart(), variable, value)
    return
end
set_start_value(x::VariableRef, v::Number) = set_start_value(x, Float64(v))

"""
    value(v::VariableRef; result = 1)

Return the value of variable `v` associated with result index `result` of the
most-recent returned by the solver.

Use [`has_values`](@ref) to check if a result exists before asking for values.

See also: [`result_count`](@ref).
"""
function value(v::VariableRef; result::Int = 1)::Float64
    return MOI.get(owner_model(v), MOI.VariablePrimal(result), v)
end

"""
    value(var_value::Function, v::VariableRef)

Evaluate the value of the variable `v` as `var_value(v)`.
"""
function value(var_value::Function, v::VariableRef)
    return var_value(v)
end

"""
    has_values(model::Model; result::Int = 1)

Return `true` if the solver has a primal solution in result index `result`
available to query, otherwise return `false`.

See also [`value`](@ref) and [`result_count`](@ref).
"""
function has_values(model::Model; result::Int = 1)
    return primal_status(model; result = result) != MOI.NO_SOLUTION
end

"""
    add_variable(m::Model, v::AbstractVariable, name::String="")

Add a variable `v` to `Model m` and sets its name.
"""
function add_variable end

function add_variable(model::Model, v::ScalarVariable, name::String = "")
    model.is_model_dirty = true
    return _moi_add_variable(backend(model), model, v, name)
end

function _moi_add_variable(moi_backend, model, v::ScalarVariable, name::String)
    index = MOI.add_variable(moi_backend)
    var_ref = VariableRef(model, index)
    _moi_constrain_variable(moi_backend, index, v.info)
    if !isempty(name) &&
       MOI.supports(moi_backend, MOI.VariableName(), MOI.VariableIndex)
        set_name(var_ref, name)
    end
    return var_ref
end

function _moi_constrain_variable(moi_backend::MOI.ModelLike, index, info)
    # We don't call the _moi* versions (e.g., _moi_set_lower_bound) because they
    # have extra checks that are not necessary for newly created variables.
    if info.has_lb
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.GreaterThan{Float64}(info.lower_bound),
        )
    end
    if info.has_ub
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.LessThan{Float64}(info.upper_bound),
        )
    end
    if info.has_fix
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.EqualTo{Float64}(info.fixed_value),
        )
    end
    if info.binary
        _moi_add_constraint(moi_backend, index, MOI.ZeroOne())
    end
    if info.integer
        _moi_add_constraint(moi_backend, index, MOI.Integer())
    end
    if info.has_start && info.start !== nothing
        MOI.set(
            moi_backend,
            MOI.VariablePrimalStart(),
            index,
            convert(Float64, info.start),
        )
    end
end

"""
    VariablesConstrainedOnCreation <: AbstractVariable

Variable `scalar_variables` constrained to belong to `set`.
Adding this variable can be understood as doing:
```julia
function JuMP.add_variable(model::Model, variable::JuMP.VariableConstrainedOnCreation, names)
    var_ref = JuMP.add_variable(model, variable.scalar_variable, name)
    JuMP.add_constraint(model, JuMP.VectorConstraint(var_ref, variable.set))
    return var_ref
end
```
but adds the variables with `MOI.add_constrained_variable(model, variable.set)`
instead. See [the MOI documentation](https://jump.dev/MathOptInterface.jl/v0.9.3/apireference/#Variables-1)
for the difference between adding the variables with `MOI.add_constrained_variable`
and adding them with `MOI.add_variable` and adding the constraint separately.
"""
struct VariableConstrainedOnCreation{
    S<:MOI.AbstractScalarSet,
    ScalarVarType<:AbstractVariable,
} <: AbstractVariable
    scalar_variable::ScalarVarType
    set::S
end

function add_variable(
    model::Model,
    variable::VariableConstrainedOnCreation,
    name::String,
)
    var_index = _moi_add_constrained_variable(
        backend(model),
        variable.scalar_variable,
        variable.set,
        name,
    )
    return VariableRef(model, var_index)
end

function add_variable(
    model::Model,
    variables::AbstractArray{<:VariableConstrainedOnCreation},
    names::AbstractArray{<:String},
)
    return add_variable.(model, variables, names)
end

function _moi_add_constrained_variable(
    moi_backend::MOI.ModelLike,
    scalar_variable::ScalarVariable,
    set::MOI.AbstractScalarSet,
    name::String,
)
    var_index, con_index = MOI.add_constrained_variable(moi_backend, set)
    _moi_constrain_variable(moi_backend, var_index, scalar_variable.info)
    if !isempty(name)
        MOI.set(moi_backend, MOI.VariableName(), var_index, name)
    end
    return var_index
end

"""
    VariablesConstrainedOnCreation <: AbstractVariable

Vector of variables `scalar_variables` constrained to belong to `set`.
Adding this variable can be thought as doing:
```julia
function JuMP.add_variable(model::Model, variable::JuMP.VariablesConstrainedOnCreation, names)
    var_refs = JuMP.add_variable.(model, variable.scalar_variables,
                                  JuMP.vectorize(names, variable.shape))
    JuMP.add_constraint(model, JuMP.VectorConstraint(var_refs, variable.set))
    return JuMP.reshape_vector(var_refs, variable.shape)
end
```
but adds the variables with `MOI.add_constrained_variables(model, variable.set)`
instead. See [the MOI documentation](https://jump.dev/MathOptInterface.jl/v0.9.3/apireference/#Variables-1)
for the difference between adding the variables with `MOI.add_constrained_variables`
and adding them with `MOI.add_variables` and adding the constraint separately.
"""
struct VariablesConstrainedOnCreation{
    S<:MOI.AbstractVectorSet,
    Shape<:AbstractShape,
    ScalarVarType<:AbstractVariable,
} <: AbstractVariable
    scalar_variables::Vector{ScalarVarType}
    set::S
    shape::Shape
end

function VariablesConstrainedOnCreation(
    variables::Vector{<:AbstractVariable},
    set::MOI.AbstractVectorSet,
)
    return VariablesConstrainedOnCreation(variables, set, VectorShape())
end

_vectorize_names(names, shape) = vectorize(names, shape)

# In some cases `names` may be a "" passed in from the macros. For now, this is
# the only case we support, so throw an assertion error if not empty.
function _vectorize_names(name::String, ::Any)
    @assert isempty(name)
    return nothing
end

function add_variable(
    model::Model,
    variable::VariablesConstrainedOnCreation,
    names,
)
    var_indices = _moi_add_constrained_variables(
        backend(model),
        variable.scalar_variables,
        variable.set,
        _vectorize_names(names, variable.shape),
    )
    var_refs = [VariableRef(model, var_index) for var_index in var_indices]
    return reshape_vector(var_refs, variable.shape)
end

function _moi_add_constrained_variables(
    moi_backend::MOI.ModelLike,
    scalar_variables::Vector{<:ScalarVariable},
    set::MOI.AbstractVectorSet,
    names::Union{Vector{String},Nothing},
)
    if set isa MOI.Reals
        var_indices = MOI.add_variables(moi_backend, MOI.dimension(set))
    else
        var_indices, con_index = MOI.add_constrained_variables(moi_backend, set)
    end
    for (index, variable) in zip(var_indices, scalar_variables)
        _moi_constrain_variable(moi_backend, index, variable.info)
    end
    if names !== nothing
        for (var_index, name) in zip(var_indices, names)
            if !isempty(name)
                MOI.set(moi_backend, MOI.VariableName(), var_index, name)
            end
        end
    end
    return var_indices
end

"""
    ComplexPlane

Complex plane object that can be used to create a complex variable in the
[`@variable`](@ref) macro.

## Examples

Consider the following example:

```jldoctest; setup = :(using JuMP)
julia> model = Model();

julia> @variable(model, x in ComplexPlane())
real(x) + (0.0 + 1.0im) imag(x)

julia> all_variables(model)
2-element Vector{VariableRef}:
 real(x)
 imag(x)
```

We see in the output of the last command that two real variables were created.
The Julia variable `x` binds to an affine expression in terms of these two
variables that parametrize the complex plane.
"""
struct ComplexPlane end

"""
    ComplexVariable{S,T,U,V} <: AbstractVariable

A struct used when adding complex variables.

See also: [`ComplexPlane`](@ref).
"""
struct ComplexVariable{S,T,U,V} <: AbstractVariable
    info::VariableInfo{S,T,U,V}
end

function build_variable(_error::Function, v::ScalarVariable, ::ComplexPlane)
    if _is_binary(v) || _is_integer(v)
        # We would then need to fix the imaginary value to zero. Let's wait to
        # see if there is need for such complication first.
        _error(
            "Creation of binary or integer complex variable is not supported.",
        )
    end
    return ComplexVariable(v.info)
end

function _mapinfo(f::Function, v::JuMP.ScalarVariable)
    info = v.info
    return ScalarVariable(
        VariableInfo(
            info.has_lb,
            f(info.lower_bound),
            info.has_ub,
            f(info.upper_bound),
            info.has_fix,
            f(info.fixed_value),
            info.has_start,
            f(info.start),
            info.binary,
            info.integer,
        ),
    )
end

function _real(s::String)
    if isempty(s)
        return s
    end
    return string("real(", s, ")")
end

function _imag(s::String)
    if isempty(s)
        return s
    end
    return string("imag(", s, ")")
end

_real(v::ScalarVariable) = _mapinfo(real, v)

_imag(v::ScalarVariable) = _mapinfo(imag, v)

_conj(v::ScalarVariable) = _mapinfo(conj, v)

function _isreal(v::ScalarVariable)
    info = v.info
    return isreal(info.lower_bound) &&
           isreal(info.upper_bound) &&
           isreal(info.fixed_value) &&
           isreal(info.start)
end

_is_binary(v::ScalarVariable) = v.info.binary

_is_integer(v::ScalarVariable) = v.info.integer

function JuMP.add_variable(
    model::JuMP.Model,
    v::ComplexVariable,
    name::String = "",
)
    model.is_model_dirty = true
    var = JuMP.ScalarVariable(v.info)
    real_part = add_variable(model, _real(var), _real(name))
    imag_part = add_variable(model, _imag(var), _imag(name))
    # Efficiently build `real_part + imag_part * im`
    return GenericAffExpr{ComplexF64,JuMP.VariableRef}(
        zero(ComplexF64),
        real_part => one(ComplexF64),
        imag_part => convert(ComplexF64, im),
    )
end

"""
    reduced_cost(x::VariableRef)::Float64

Return the reduced cost associated with variable `x`.

Equivalent to querying the shadow price of the active variable bound
(if one exists and is active).

See also: [`shadow_price`](@ref).
"""
function reduced_cost(x::VariableRef)::Float64
    model = owner_model(x)
    if !has_duals(model)
        error(
            "Unable to query reduced cost of variable because model does" *
            " not have duals available.",
        )
    end
    sign = objective_sense(model) == MIN_SENSE ? 1.0 : -1.0
    if is_fixed(x)
        return sign * dual(FixRef(x))
    end
    rc = 0.0
    if has_upper_bound(x)
        rc += dual(UpperBoundRef(x))
    end
    if has_lower_bound(x)
        rc += dual(LowerBoundRef(x))
    end
    return sign * rc
end

"""
    all_variables(model::Model)::Vector{VariableRef}

Returns a list of all variables currently in the model. The variables are
ordered by creation time.

# Example
```jldoctest; setup=:(using JuMP)
model = Model()
@variable(model, x)
@variable(model, y)
all_variables(model)

# output

2-element Array{VariableRef,1}:
 x
 y
```
"""
function all_variables(model::Model)
    all_indices =
        MOI.get(model, MOI.ListOfVariableIndices())::Vector{MOI.VariableIndex}
    return VariableRef[VariableRef(model, idx) for idx in all_indices]
end

function dual(vref::VariableRef)
    return error(
        "To query the dual variables associated with a variable bound, first " *
        "obtain a constraint reference using one of `UpperBoundRef`, `LowerBoundRef`, " *
        "or `FixRef`, and then call `dual` on the returned constraint reference.\nFor " *
        "example, if `x <= 1`, instead of `dual(x)`, call `dual(UpperBoundRef(x))`.",
    )
end

function value(::AbstractArray{<:AbstractJuMPScalar})
    return error(
        "`JuMP.value` is not defined for collections of JuMP types. Use" *
        " Julia's broadcast syntax instead: `JuMP.value.(x)`.",
    )
end

value(::_MA.Zero) = 0.0
value(x::Number) = x

function _info_from_variable(v::VariableRef)
    has_lb = has_lower_bound(v)
    lb = has_lb ? lower_bound(v) : -Inf
    has_ub = has_upper_bound(v)
    ub = has_ub ? upper_bound(v) : Inf
    has_fix = is_fixed(v)
    fixed_value = has_fix ? fix_value(v) : NaN
    has_start, start = false, NaN
    if MOI.supports(
        backend(owner_model(v)),
        MOI.VariablePrimalStart(),
        MOI.VariableIndex,
    )
        start = start_value(v)
        has_start = start !== nothing
    end
    binary = is_binary(v)
    integer = is_integer(v)
    return VariableInfo(
        has_lb,
        lb,
        has_ub,
        ub,
        has_fix,
        fixed_value,
        has_start,
        start,
        binary,
        integer,
    )
end

"""
    relax_integrality(model::Model)

Modifies `model` to "relax" all binary and integrality constraints on
variables. Specifically,

- Binary constraints are deleted, and variable bounds are tightened if
  necessary to ensure the variable is constrained to the interval ``[0, 1]``.
- Integrality constraints are deleted without modifying variable bounds.
- An error is thrown if semi-continuous or semi-integer constraints are
  present (support may be added for these in the future).
- All other constraints are ignored (left in place). This includes discrete
  constraints like SOS and indicator constraints.

Returns a function that can be called without any arguments to restore the
original model. The behavior of this function is undefined if additional
changes are made to the affected variables in the meantime.

# Example
```jldoctest; setup=:(using JuMP)
julia> model = Model();

julia> @variable(model, x, Bin);

julia> @variable(model, 1 <= y <= 10, Int);

julia> @objective(model, Min, x + y);

julia> undo_relax = relax_integrality(model);

julia> print(model)
Min x + y
Subject to
 x ≥ 0.0
 y ≥ 1.0
 x ≤ 1.0
 y ≤ 10.0

julia> undo_relax()

julia> print(model)
Min x + y
Subject to
 y ≥ 1.0
 y ≤ 10.0
 y integer
 x binary
```
"""
function relax_integrality(model::Model)
    if num_constraints(model, VariableRef, MOI.Semicontinuous{Float64}) > 0
        error(
            "Support for relaxing semicontinuous constraints is not " *
            "yet implemented.",
        )
    end
    if num_constraints(model, VariableRef, MOI.Semiinteger{Float64}) > 0
        error(
            "Support for relaxing semi-integer constraints is not " *
            "yet implemented.",
        )
    end

    bin_int_constraints = vcat(
        all_constraints(model, VariableRef, MOI.ZeroOne),
        all_constraints(model, VariableRef, MOI.Integer),
    )
    bin_int_variables = VariableRef.(bin_int_constraints)
    info_pre_relaxation =
        map(v -> (v, _info_from_variable(v)), bin_int_variables)
    # We gather the info first because some solvers perform poorly when you
    # interleave queries and changes. See, e.g.,
    # https://github.com/jump-dev/Gurobi.jl/pull/301.
    for (v, info) in info_pre_relaxation
        if info.integer
            unset_integer(v)
        elseif info.binary
            unset_binary(v)
            if !info.has_fix
                set_lower_bound(v, max(0.0, info.lower_bound))
                set_upper_bound(v, min(1.0, info.upper_bound))
            elseif info.fixed_value < 0 || info.fixed_value > 1
                error(
                    "The model has no valid relaxation: binary variable " *
                    "fixed out of bounds.",
                )
            end
        end
    end
    function unrelax()
        for (v, info) in info_pre_relaxation
            if info.integer
                set_integer(v)
            elseif info.binary
                set_binary(v)
                if !info.has_fix
                    if info.has_lb
                        set_lower_bound(v, info.lower_bound)
                    else
                        delete_lower_bound(v)
                    end
                    if info.has_ub
                        set_upper_bound(v, info.upper_bound)
                    else
                        delete_upper_bound(v)
                    end
                end
            end
        end
        return
    end
    return unrelax
end

###
### Error messages for common incorrect usages
###

for sym in (:(<=), :(>=), :(<), :(>))
    msg = """Cannot evaluate `$(sym)` between a variable and a number.

    There are two common mistakes that lead to this.

     * You tried to write a constraint that depends on the value of a variable

       For example:
       ```julia
       model = Model()
       @variable(model, x[1:2])
       if x[1] $(sym) 1
           @constraint(model, x[2] == 0)
       end
       ```

       You cannot write a model like this. You must formulate your problem as a
       single optimization problem. Unfortunately, the way to do this is
       problem-specific and depends on your choice of solver. You may be able to
       use indicator constraints, or some other mixed-integer linear
       reformulation. If stuck, post your problem on the community forum:
       https://discourse.julialang.org/c/domain/opt/13

     * You wrote a function that expected the value of a variable, but it was
       passed the variable instead

       For example:
       ```julia
       foo(x) = x $(sym) 1 ? 0 : 1 - x
       model = Model()
       @variable(model, x)
       @objective(model, foo(x))
       ```

       To fix this, create a nonlinear model with a user-defined function:
       ```julia
       foo(x) = x $(sym) 1 ? 0 : 1 - x
       model = Model()
       register(model, :foo, 1, foo; autodiff = true)
       @variable(model, x)
       @NLobjective(model, foo(x))
       ```
    """
    @eval begin
        Base.$(sym)(::VariableRef, ::Number) = error($(msg))
        Base.$(sym)(::Number, ::VariableRef) = error($(msg))
    end
end
