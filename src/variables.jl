#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""
    num_variables(model::GenericModel)::Int64

Returns number of variables in `model`.
"""
num_variables(model::GenericModel)::Int64 =
    MOI.get(model, MOI.NumberOfVariables())

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
    if info.has_lb
        _error("Cannot specify variable lower_bound twice")
    end
    info.has_lb = true
    info.lower_bound = lower
    return
end

function _set_upper_bound_or_error(
    _error::Function,
    info::_VariableInfoExpr,
    upper,
)
    if info.has_ub
        _error("Cannot specify variable upper_bound twice")
    end
    info.has_ub = true
    info.upper_bound = upper
    return
end

function _fix_or_error(_error::Function, info::_VariableInfoExpr, value)
    if info.has_fix
        _error("Cannot specify variable fixed value twice")
    end
    info.has_fix = true
    info.fixed_value = value
    return
end

function _set_binary_or_error(_error::Function, info::_VariableInfoExpr)
    if info.binary
        _error("'Bin' and 'binary' keyword argument cannot both be specified.")
    end
    info.binary = true
    return
end

function _set_integer_or_error(_error::Function, info::_VariableInfoExpr)
    if info.integer
        _error("'Int' and 'integer' keyword argument cannot both be specified.")
    end
    info.integer = true
    return
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

"""
    variable_ref_type(::Union{F,Type{F}}) where {F}

A helper function used internally by JuMP and some JuMP extensions. Returns the
variable type associated with the model or expression type `F`.
"""
variable_ref_type(::F) where {F} = variable_ref_type(F)

variable_ref_type(::Type{V}) where {V<:AbstractVariableRef} = V

value_type(::Type{<:AbstractVariableRef}) = Float64

Base.conj(v::AbstractVariableRef) = v
Base.real(v::AbstractVariableRef) = v
Base.imag(v::AbstractVariableRef) = zero(v)
Base.abs2(v::AbstractVariableRef) = v^2
Base.isreal(::AbstractVariableRef) = true

"""
    GenericVariableRef <: AbstractVariableRef

Holds a reference to the model and the corresponding MOI.VariableIndex.
"""
struct GenericVariableRef{T} <: AbstractVariableRef
    model::GenericModel{T}
    index::MOI.VariableIndex
end

const VariableRef = GenericVariableRef{Float64}

value_type(::Type{GenericVariableRef{T}}) where {T} = T

variable_ref_type(::Type{GenericModel{T}}) where {T} = GenericVariableRef{T}

# `AbstractVariableRef` types must override the default `owner_model` if the field
#  name is not `model`.
"""
    owner_model(v::AbstractVariableRef)

Returns the model to which `v` belongs.

## Example

```jldoctest
julia> model = Model();

julia> x = @variable(model)
_[1]

julia> owner_model(x) === model
true
```
"""
owner_model(v::AbstractVariableRef) = v.model

"""
    struct VariableNotOwned{V<:AbstractVariableRef} <: Exception
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

Base.iszero(::GenericVariableRef) = false
function Base.copy(v::GenericVariableRef{T}) where {T}
    return GenericVariableRef{T}(v.model, v.index)
end
Base.broadcastable(v::GenericVariableRef) = Ref(v)

Base.zero(v::AbstractVariableRef) = zero(typeof(v))

function Base.zero(::Type{V}) where {V<:AbstractVariableRef}
    return zero(GenericAffExpr{value_type(V),V})
end

Base.one(v::AbstractVariableRef) = one(typeof(v))

function Base.one(::Type{V}) where {V<:AbstractVariableRef}
    return one(GenericAffExpr{value_type(V),V})
end

"""
    coefficient(v1::GenericVariableRef{T}, v2::GenericVariableRef{T}) where {T}

Return `one(T)` if `v1 == v2`, and `zero(T)` otherwise.

This is a fallback for other [`coefficient`](@ref) methods to simplify code in
which the expression may be a single variable.
"""
function coefficient(
    v1::GenericVariableRef{T},
    v2::GenericVariableRef{T},
) where {T}
    if v1 == v2
        return one(T)
    else
        return zero(T)
    end
end

function coefficient(
    ::GenericVariableRef{T},
    ::GenericVariableRef{T},
    ::GenericVariableRef{T},
) where {T}
    return zero(T)
end

function isequal_canonical(v::GenericVariableRef, other::GenericVariableRef)
    return isequal(v, other)
end

"""
    delete(model::GenericModel, variable_ref::GenericVariableRef)

Delete the variable associated with `variable_ref` from the model `model`.

Note that `delete` does not unregister the name from the model, so adding a new
variable of the same name will throw an error. Use [`unregister`](@ref) to
unregister the name after deletion.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x)
x

julia> delete(model, x)

julia> unregister(model, :x)

julia> print(model)
Feasibility
Subject to

julia> model[:x]
ERROR: KeyError: key :x not found
Stacktrace:
[...]
```
"""
function delete(model::GenericModel, variable_ref::GenericVariableRef)
    if model !== owner_model(variable_ref)
        error(
            "The variable reference you are trying to delete does not " *
            "belong to the model.",
        )
    end
    model.is_model_dirty = true
    MOI.delete(backend(model), variable_ref.index)
    return
end

"""
    delete(model::GenericModel, variable_refs::Vector{<:GenericVariableRef})

Delete the variables associated with `variable_refs` from the model `model`.
Solvers may implement methods for deleting multiple variables that are
more efficient than repeatedly calling the single variable delete method.

See also: [`unregister`](@ref)
"""
function delete(
    model::GenericModel,
    variable_refs::Vector{<:GenericVariableRef},
)
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
    is_valid(model::GenericModel, variable_ref::GenericVariableRef)

Return `true` if `variable` refers to a valid variable in `model`.
"""
function is_valid(model::GenericModel, variable_ref::GenericVariableRef)
    return model === owner_model(variable_ref) &&
           MOI.is_valid(backend(model), variable_ref.index)
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/jump-dev/MathOptInterface.jl/issues/234#issuecomment-366868878
function Base.hash(v::GenericVariableRef, h::UInt)
    return hash(objectid(owner_model(v)), hash(v.index.value, h))
end

function Base.isequal(v1::GenericVariableRef, v2::GenericVariableRef)
    return owner_model(v1) === owner_model(v2) && v1.index == v2.index
end

"""
    index(v::GenericVariableRef)::MOI.VariableIndex

Return the index of the variable that corresponds to `v` in the MOI backend.
"""
index(v::GenericVariableRef) = v.index

function GenericVariableRef{T}(model::GenericModel{T}) where {T}
    index = MOI.add_variable(backend(model))
    return GenericVariableRef{T}(model, index)
end

"""
    GenericVariableRef(c::ConstraintRef)

Get the variable associated with a `ConstraintRef`, if `c` is a constraint on a
single variable.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x >= 0)
x

julia> c = LowerBoundRef(x)
x ≥ 0.0

julia> GenericVariableRef(c) == x
true
```
"""
function GenericVariableRef{T}(
    c::ConstraintRef{GenericModel{T},<:MOI.ConstraintIndex{MOI.VariableIndex}},
) where {T}
    vi = MOI.VariableIndex(index(c).value)
    return GenericVariableRef{T}(owner_model(c), vi)
end

# Name setter/getters
# These functions need to be implemented for all `AbstractVariableRef`s
"""
    name(v::GenericVariableRef)::String

Get a variable's name attribute.
"""
function name(v::GenericVariableRef)
    model = owner_model(v)
    if !MOI.supports(backend(model), MOI.VariableName(), MOI.VariableIndex)
        return ""
    end
    return MOI.get(model, MOI.VariableName(), v)::String
end

"""
    set_name(v::GenericVariableRef, s::AbstractString)

Set a variable's name attribute.
"""
function set_name(v::GenericVariableRef, s::String)
    MOI.set(owner_model(v), MOI.VariableName(), v, s)
    return
end

"""
    variable_by_name(
        model::AbstractModel,
        name::String,
    )::Union{AbstractVariableRef,Nothing}

Returns the reference of the variable with name attribute `name` or `Nothing` if
no variable has this name attribute. Throws an error if several variables have
`name` as their name attribute.

## Examples

```jldoctest objective_function; filter = r"Stacktrace:.*"s
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
 [5] variable_by_name(::GenericModel, ::String) at /home/blegat/.julia/dev/JuMP/src/variables.jl:268
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
function variable_by_name(model::GenericModel, name::String)
    index = MOI.get(backend(model), MOI.VariableIndex, name)
    if index === nothing
        return nothing
    end
    return GenericVariableRef(model, index)
end

MOI.VariableIndex(v::GenericVariableRef) = index(v)

moi_function(variable::AbstractVariableRef) = index(variable)

moi_function_type(::Type{<:AbstractVariableRef}) = MOI.VariableIndex

# Note: No validation is performed that the variables belong to the same model.
function MOI.VectorOfVariables(vars::Vector{<:GenericVariableRef})
    return MOI.VectorOfVariables(index.(vars))
end

function moi_function(variables::Vector{<:AbstractVariableRef})
    return MOI.VectorOfVariables(variables)
end

function moi_function_type(::Type{<:Vector{<:AbstractVariableRef}})
    return MOI.VectorOfVariables
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VectorOfVariables},
) where {T}
    return Vector{GenericVariableRef{T}}
end

function jump_function(
    model::GenericModel{T},
    variables::MOI.VectorOfVariables,
) where {T}
    return GenericVariableRef{T}[
        GenericVariableRef{T}(model, v) for v in variables.variables
    ]
end

function jump_function_type(
    ::GenericModel{T},
    ::Type{MOI.VariableIndex},
) where {T}
    return GenericVariableRef{T}
end

function jump_function(
    model::GenericModel{T},
    variable::MOI.VariableIndex,
) where {T}
    return GenericVariableRef{T}(model, variable)
end

## Bound setter/getters

# lower bounds

"""
    has_lower_bound(v::GenericVariableRef)

Return `true` if `v` has a lower bound. If `true`, the lower bound can be
queried with [`lower_bound`](@ref).

See also [`LowerBoundRef`](@ref), [`lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function has_lower_bound(v::GenericVariableRef)
    return _moi_has_lower_bound(backend(owner_model(v)), v)
end

# _moi_* methods allow us to work around the type instability of the backend of
# a model.
function _moi_has_lower_bound(moi_backend, v::GenericVariableRef)
    return MOI.is_valid(moi_backend, _lower_bound_index(v))
end

function _lower_bound_index(v::GenericVariableRef{T}) where {T}
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{T}}(
        index(v).value,
    )
end

"""
    set_lower_bound(v::GenericVariableRef, lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower
bound constraint.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function set_lower_bound(v::GenericVariableRef, lower::Number)
    if !isfinite(lower)
        error(
            "Unable to set lower bound to $(lower). To remove the bound, use " *
            "`delete_lower_bound`.",
        )
    end
    model = owner_model(v)
    model.is_model_dirty = true
    _moi_set_lower_bound(backend(model), v, lower)
    return
end

function _moi_set_lower_bound(
    moi_backend,
    v::GenericVariableRef{T},
    lower::Number,
) where {T}
    new_set = MOI.GreaterThan(convert(T, lower))
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
    LowerBoundRef(v::GenericVariableRef)

Return a constraint reference to the lower bound constraint of `v`. Errors if
one does not exist.

See also [`has_lower_bound`](@ref), [`lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function LowerBoundRef(v::GenericVariableRef)
    return ConstraintRef(owner_model(v), _lower_bound_index(v), ScalarShape())
end

"""
    delete_lower_bound(v::GenericVariableRef)

Delete the lower bound constraint of a variable.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`lower_bound`](@ref), [`set_lower_bound`](@ref).
"""
function delete_lower_bound(variable_ref::GenericVariableRef)
    delete(owner_model(variable_ref), LowerBoundRef(variable_ref))
    return
end

"""
    lower_bound(v::GenericVariableRef)

Return the lower bound of a variable. Error if one does not exist.

See also [`LowerBoundRef`](@ref), [`has_lower_bound`](@ref),
[`set_lower_bound`](@ref), [`delete_lower_bound`](@ref).
"""
function lower_bound(v::GenericVariableRef{T}) where {T}
    if !has_lower_bound(v)
        error("Variable $(v) does not have a lower bound.")
    end
    set = MOI.get(owner_model(v), MOI.ConstraintSet(), LowerBoundRef(v))
    return set.lower::T
end

# upper bounds

"""
    has_upper_bound(v::GenericVariableRef)

Return `true` if `v` has a upper bound. If `true`, the upper bound can be
queried with [`upper_bound`](@ref).

See also [`UpperBoundRef`](@ref), [`upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function has_upper_bound(v::GenericVariableRef)
    return _moi_has_upper_bound(backend(owner_model(v)), v)
end

function _moi_has_upper_bound(moi_backend, v::GenericVariableRef)
    return MOI.is_valid(moi_backend, _upper_bound_index(v))
end

function _upper_bound_index(v::GenericVariableRef{T}) where {T}
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{T}}(
        index(v).value,
    )
end

"""
    set_upper_bound(v::GenericVariableRef, upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound
constraint.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function set_upper_bound(v::GenericVariableRef, upper::Number)
    if !isfinite(upper)
        error(
            "Unable to set upper bound to $(upper). To remove the bound, use " *
            "`delete_upper_bound`.",
        )
    end
    model = owner_model(v)
    model.is_model_dirty = true
    _moi_set_upper_bound(backend(model), v, upper)
    return
end

function _moi_set_upper_bound(
    moi_backend,
    v::GenericVariableRef{T},
    upper::Number,
) where {T}
    new_set = MOI.LessThan(convert(T, upper))
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
    UpperBoundRef(v::GenericVariableRef)

Return a constraint reference to the upper bound constraint of `v`. Errors if
one does not exist.

See also [`has_upper_bound`](@ref), [`upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function UpperBoundRef(v::GenericVariableRef)
    return ConstraintRef(owner_model(v), _upper_bound_index(v), ScalarShape())
end

"""
    delete_upper_bound(v::GenericVariableRef)

Delete the upper bound constraint of a variable.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`upper_bound`](@ref), [`set_upper_bound`](@ref).
"""
function delete_upper_bound(variable_ref::GenericVariableRef)
    delete(owner_model(variable_ref), UpperBoundRef(variable_ref))
    return
end

"""
    upper_bound(v::GenericVariableRef)

Return the upper bound of a variable. Error if one does not exist.

See also [`UpperBoundRef`](@ref), [`has_upper_bound`](@ref),
[`set_upper_bound`](@ref), [`delete_upper_bound`](@ref).
"""
function upper_bound(v::GenericVariableRef{T}) where {T}
    if !has_upper_bound(v)
        error("Variable $(v) does not have an upper bound.")
    end
    set = MOI.get(owner_model(v), MOI.ConstraintSet(), UpperBoundRef(v))
    return set.upper::T
end

# fixed value

"""
    is_fixed(v::GenericVariableRef)

Return `true` if `v` is a fixed variable. If `true`, the fixed value can be
queried with [`fix_value`](@ref).

See also [`FixRef`](@ref), [`fix_value`](@ref), [`fix`](@ref), [`unfix`](@ref).
"""
function is_fixed(v::GenericVariableRef)
    return _moi_is_fixed(backend(owner_model(v)), v)
end

function _moi_is_fixed(moi_backend, v::GenericVariableRef)
    return MOI.is_valid(moi_backend, _fix_index(v))
end

function _fix_index(v::GenericVariableRef{T}) where {T}
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{T}}(index(v).value)
end

"""
    fix(v::GenericVariableRef, value::Number; force::Bool = false)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise
create a new one.

If the variable already has variable bounds and `force=false`, calling `fix`
will throw an error. If `force=true`, existing variable bounds will be deleted,
and the fixing constraint will be added. Note a variable will have no bounds
after a call to [`unfix`](@ref).

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix_value`](@ref),
[`unfix`](@ref).
"""
function fix(variable::GenericVariableRef, value::Number; force::Bool = false)
    if !isfinite(value)
        error("Unable to fix variable to $(value)")
    end
    model = owner_model(variable)
    model.is_model_dirty = true
    _moi_fix(backend(model), variable, value, force)
    return
end

function _moi_fix(
    moi_backend,
    variable::GenericVariableRef{T},
    value::Number,
    force::Bool,
) where {T}
    new_set = MOI.EqualTo(convert(T, value))
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
    unfix(v::GenericVariableRef)

Delete the fixing constraint of a variable.

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix_value`](@ref),
[`fix`](@ref).
"""
function unfix(variable_ref::GenericVariableRef)
    delete(owner_model(variable_ref), FixRef(variable_ref))
    return
end

"""
    fix_value(v::GenericVariableRef)

Return the value to which a variable is fixed. Error if one does not exist.

See also [`FixRef`](@ref), [`is_fixed`](@ref), [`fix`](@ref), [`unfix`](@ref).
"""
function fix_value(v::GenericVariableRef{T}) where {T}
    set = MOI.get(owner_model(v), MOI.ConstraintSet(), FixRef(v))
    return set.value::T
end

"""
    FixRef(v::GenericVariableRef)

Return a constraint reference to the constraint fixing the value of `v`. Errors
if one does not exist.

See also [`is_fixed`](@ref), [`fix_value`](@ref), [`fix`](@ref),
[`unfix`](@ref).
"""
function FixRef(v::GenericVariableRef)
    return ConstraintRef(owner_model(v), _fix_index(v), ScalarShape())
end

"""
    is_integer(v::GenericVariableRef)

Return `true` if `v` is constrained to be integer.

See also [`IntegerRef`](@ref), [`set_integer`](@ref), [`unset_integer`](@ref).
"""
function is_integer(v::GenericVariableRef)
    return _moi_is_integer(backend(owner_model(v)), v)
end

function _moi_is_integer(moi_backend, v::GenericVariableRef)
    return MOI.is_valid(moi_backend, _integer_index(v))
end

function _integer_index(v::GenericVariableRef)
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(index(v).value)
end

"""
    set_integer(variable_ref::GenericVariableRef)

Add an integrality constraint on the variable `variable_ref`.

See also [`IntegerRef`](@ref), [`is_integer`](@ref), [`unset_integer`](@ref).
"""
function set_integer(v::GenericVariableRef)
    model = owner_model(v)
    model.is_model_dirty = true
    _moi_set_integer(backend(model), v)
    return
end

function _moi_set_integer(moi_backend, variable_ref::GenericVariableRef)
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
    unset_integer(variable_ref::GenericVariableRef)

Remove the integrality constraint on the variable `variable_ref`.

See also [`IntegerRef`](@ref), [`is_integer`](@ref), [`set_integer`](@ref).
"""
function unset_integer(variable_ref::GenericVariableRef)
    delete(owner_model(variable_ref), IntegerRef(variable_ref))
    return
end

"""
    IntegerRef(v::GenericVariableRef)

Return a constraint reference to the constraint constraining `v` to be integer.
Errors if one does not exist.

See also [`is_integer`](@ref), [`set_integer`](@ref), [`unset_integer`](@ref).
"""
function IntegerRef(v::GenericVariableRef)
    return ConstraintRef(owner_model(v), _integer_index(v), ScalarShape())
end

"""
    is_binary(v::GenericVariableRef)

Return `true` if `v` is constrained to be binary.

See also [`BinaryRef`](@ref), [`set_binary`](@ref), [`unset_binary`](@ref).
"""
function is_binary(v::GenericVariableRef)
    return _moi_is_binary(backend(owner_model(v)), v)
end

function _moi_is_binary(moi_backend, v::GenericVariableRef)
    return MOI.is_valid(moi_backend, _binary_index(v))
end

function _binary_index(v::GenericVariableRef)
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(index(v).value)
end

"""
    set_binary(v::GenericVariableRef)

Add a constraint on the variable `v` that it must take values in the set
``\\{0,1\\}``.

See also [`BinaryRef`](@ref), [`is_binary`](@ref), [`unset_binary`](@ref).
"""
function set_binary(v::GenericVariableRef)
    model = owner_model(v)
    model.is_model_dirty = true
    _moi_set_binary(backend(model), v)
    return
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
    unset_binary(variable_ref::GenericVariableRef)

Remove the binary constraint on the variable `variable_ref`.

See also [`BinaryRef`](@ref), [`is_binary`](@ref), [`set_binary`](@ref).
"""
function unset_binary(variable_ref::GenericVariableRef)
    delete(owner_model(variable_ref), BinaryRef(variable_ref))
    return
end

"""
    BinaryRef(v::GenericVariableRef)

Return a constraint reference to the constraint constraining `v` to be binary.
Errors if one does not exist.

See also [`is_binary`](@ref), [`set_binary`](@ref), [`unset_binary`](@ref).
"""
function BinaryRef(v::GenericVariableRef)
    return ConstraintRef(owner_model(v), _binary_index(v), ScalarShape())
end

"""
    start_value(v::GenericVariableRef)

Return the start value (MOI attribute `VariablePrimalStart`) of the variable
`v`.

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".

See also [`set_start_value`](@ref).
"""
function start_value(v::GenericVariableRef{T})::Union{Nothing,T} where {T}
    return MOI.get(owner_model(v), MOI.VariablePrimalStart(), v)
end

"""
    has_start_value(variable::AbstractVariableRef)

Return `true` if the variable has a start value set otherwise return `false`.

See also [`set_start_value`](@ref).
"""
has_start_value(v::AbstractVariableRef)::Bool = start_value(v) !== nothing

_convert_if_something(::Type{T}, x) where {T} = convert(T, x)
_convert_if_something(::Type, ::Nothing) = nothing

"""
    set_start_value(variable::GenericVariableRef, value::Union{Real,Nothing})

Set the start value (MOI attribute `VariablePrimalStart`) of the `variable` to
`value`.

Pass `nothing` to unset the start value.

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".

See also [`start_value`](@ref).
"""
function set_start_value(
    variable::GenericVariableRef{T},
    value::Union{Nothing,Real},
) where {T}
    MOI.set(
        owner_model(variable),
        MOI.VariablePrimalStart(),
        variable,
        _convert_if_something(T, value),
    )
    return
end

function set_start_value(x::GenericVariableRef{T}, v::Number) where {T}
    set_start_value(x, convert(T, v))
    return
end

"""
    value(v::GenericVariableRef; result = 1)

Return the value of variable `v` associated with result index `result` of the
most-recent returned by the solver.

Use [`has_values`](@ref) to check if a result exists before asking for values.

See also: [`result_count`](@ref).
"""
function value(v::GenericVariableRef{T}; result::Int = 1)::T where {T}
    return MOI.get(owner_model(v), MOI.VariablePrimal(result), v)
end

"""
    value(var_value::Function, v::GenericVariableRef)

Evaluate the value of the variable `v` as `var_value(v)`.
"""
function value(var_value::Function, v::GenericVariableRef)
    return var_value(v)
end

"""
    has_values(model::GenericModel; result::Int = 1)

Return `true` if the solver has a primal solution in result index `result`
available to query, otherwise return `false`.

See also [`value`](@ref) and [`result_count`](@ref).
"""
function has_values(model::GenericModel; result::Int = 1)
    return primal_status(model; result = result) != MOI.NO_SOLUTION
end

"""
    add_variable(m::GenericModel, v::AbstractVariable, name::String="")

Add a variable `v` to `Model m` and sets its name.
"""
function add_variable end

function add_variable(model::GenericModel, v::ScalarVariable, name::String = "")
    model.is_model_dirty = true
    return _moi_add_variable(backend(model), model, v, name)
end

function _moi_add_variable(
    moi_backend,
    model::GenericModel{T},
    v::ScalarVariable,
    name::String,
) where {T}
    index = MOI.add_variable(moi_backend)
    var_ref = GenericVariableRef(model, index)
    _moi_constrain_variable(moi_backend, index, v.info, T)
    if !isempty(name) &&
       MOI.supports(moi_backend, MOI.VariableName(), MOI.VariableIndex)
        set_name(var_ref, name)
    end
    return var_ref
end

function _moi_constrain_variable(
    moi_backend::MOI.ModelLike,
    index,
    info,
    ::Type{T},
) where {T}
    # We don't call the _moi* versions (e.g., _moi_set_lower_bound) because they
    # have extra checks that are not necessary for newly created variables.
    if info.has_lb
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.GreaterThan{T}(info.lower_bound),
        )
    end
    if info.has_ub
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.LessThan{T}(info.upper_bound),
        )
    end
    if info.has_fix
        _moi_add_constraint(
            moi_backend,
            index,
            MOI.EqualTo{T}(info.fixed_value),
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
            convert(T, info.start),
        )
    end
end

"""
    VariablesConstrainedOnCreation <: AbstractVariable

Variable `scalar_variables` constrained to belong to `set`.

Adding this variable can be understood as doing:
```julia
function JuMP.add_variable(
    model::GenericModel,
    variable::VariableConstrainedOnCreation,
    names,
)
    var_ref = add_variable(model, variable.scalar_variable, name)
    add_constraint(model, VectorConstraint(var_ref, variable.set))
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
    model::GenericModel{T},
    variable::VariableConstrainedOnCreation,
    name::String,
) where {T}
    var_index = _moi_add_constrained_variable(
        backend(model),
        variable.scalar_variable,
        variable.set,
        name,
        T,
    )
    return GenericVariableRef(model, var_index)
end

function add_variable(
    model::GenericModel,
    variables::AbstractArray{<:VariableConstrainedOnCreation},
    names::AbstractArray{<:String},
)
    return add_variable.(model, variables, names)
end

function add_variable(
    model::GenericModel,
    variables::AbstractArray{<:VariableConstrainedOnCreation},
    name::String,
)
    return add_variable.(model, variables, Ref(name))
end

function _moi_add_constrained_variable(
    moi_backend::MOI.ModelLike,
    scalar_variable::ScalarVariable,
    set::MOI.AbstractScalarSet,
    name::String,
    ::Type{T},
) where {T}
    var_index, con_index = MOI.add_constrained_variable(moi_backend, set)
    _moi_constrain_variable(moi_backend, var_index, scalar_variable.info, T)
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
function JuMP.add_variable(
    model::GenericModel,
    variable::VariablesConstrainedOnCreation,
    names,
)
    v_names = vectorize(names, variable.shape)
    var_refs = add_variable.(model, variable.scalar_variables, v_names)
    add_constraint(model, VectorConstraint(var_refs, variable.set))
    return reshape_vector(var_refs, variable.shape)
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
    return
end

function add_variable(
    model::GenericModel{T},
    variable::VariablesConstrainedOnCreation,
    names,
) where {T}
    var_indices = _moi_add_constrained_variables(
        backend(model),
        variable.scalar_variables,
        variable.set,
        _vectorize_names(names, variable.shape),
        T,
    )
    var_refs =
        [GenericVariableRef{T}(model, var_index) for var_index in var_indices]
    return reshape_vector(var_refs, variable.shape)
end

function _moi_add_constrained_variables(
    moi_backend::MOI.ModelLike,
    scalar_variables::Vector{<:ScalarVariable},
    set::MOI.AbstractVectorSet,
    names::Union{Vector{String},Nothing},
    ::Type{T},
) where {T}
    if set isa MOI.Reals
        var_indices = MOI.add_variables(moi_backend, MOI.dimension(set))
    else
        var_indices, con_index = MOI.add_constrained_variables(moi_backend, set)
    end
    for (index, variable) in zip(var_indices, scalar_variables)
        _moi_constrain_variable(moi_backend, index, variable.info, T)
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

## Example

Consider the following example:

```jldoctest
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

function _mapinfo(f::Function, v::ScalarVariable)
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
_real(scalar::AbstractJuMPScalar) = real(scalar)

_imag(v::ScalarVariable) = _mapinfo(imag, v)
_imag(scalar::AbstractJuMPScalar) = imag(scalar)

_conj(v::ScalarVariable) = _mapinfo(conj, v)

function _isreal(v::ScalarVariable)
    return isreal(v.info.lower_bound) &&
           isreal(v.info.upper_bound) &&
           isreal(v.info.fixed_value) &&
           isreal(v.info.start)
end

_is_binary(v::ScalarVariable) = v.info.binary

_is_integer(v::ScalarVariable) = v.info.integer

function add_variable(
    model::GenericModel{T},
    v::ComplexVariable,
    name::String = "",
) where {T}
    model.is_model_dirty = true
    var = ScalarVariable(v.info)
    real_part = add_variable(model, _real(var), _real(name))
    imag_part = add_variable(model, _imag(var), _imag(name))
    # Efficiently build `real_part + imag_part * im`
    return GenericAffExpr{ComplexF64,GenericVariableRef{T}}(
        zero(ComplexF64),
        real_part => one(ComplexF64),
        imag_part => convert(ComplexF64, im),
    )
end

function build_variable(
    _error::Function,
    variables::AbstractArray{<:ScalarVariable},
    set::ComplexPlane,
)
    return build_variable.(_error, variables, Ref(set))
end

function add_variable(
    model::GenericModel,
    variables::AbstractArray{<:ComplexVariable},
    name::Union{<:AbstractArray{String},String} = "",
)
    return add_variable.(model, variables, name)
end

"""
    reduced_cost(x::GenericVariableRef{T})::T where {T}

Return the reduced cost associated with variable `x`.

Equivalent to querying the shadow price of the active variable bound
(if one exists and is active).

See also: [`shadow_price`](@ref).
"""
function reduced_cost(x::GenericVariableRef{T})::T where {T}
    model = owner_model(x)
    if !has_duals(model)
        error(
            "Unable to query reduced cost of variable because model does" *
            " not have duals available.",
        )
    end
    sign = objective_sense(model) == MIN_SENSE ? one(T) : -one(T)
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
    all_variables(model::GenericModel{T})::Vector{GenericVariableRef{T}} where {T}

Returns a list of all variables currently in the model. The variables are
ordered by creation time.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x);

julia> @variable(model, y);

julia> all_variables(model)
2-element Vector{VariableRef}:
 x
 y
```
"""
function all_variables(model::GenericModel{T}) where {T}
    all_indices =
        MOI.get(model, MOI.ListOfVariableIndices())::Vector{MOI.VariableIndex}
    return GenericVariableRef{T}[
        GenericVariableRef(model, idx) for idx in all_indices
    ]
end

function dual(::GenericVariableRef)
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

function _info_from_variable(v::GenericVariableRef)
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
    relax_integrality(model::GenericModel)

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

## Example

```jldoctest
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
function relax_integrality(model::GenericModel)
    return _relax_or_fix_integrality(nothing, model)
end

"""
    fix_discrete_variables([var_value::Function = value,] model::GenericModel)

Modifies `model` to convert all binary and integer variables to continuous
variables with fixed bounds of `var_value(x)`.

## Return

Returns a function that can be called without any arguments to restore the
original model. The behavior of this function is undefined if additional
changes are made to the affected variables in the meantime.

## Notes

- An error is thrown if semi-continuous or semi-integer constraints are
  present (support may be added for these in the future).
- All other constraints are ignored (left in place). This includes discrete
  constraints like SOS and indicator constraints.

## Example

```jldoctest
julia> model = Model();

julia> @variable(model, x, Bin, start = 1);

julia> @variable(model, 1 <= y <= 10, Int, start = 2);

julia> @objective(model, Min, x + y);

julia> undo_relax = fix_discrete_variables(start_value, model);

julia> print(model)
Min x + y
Subject to
 x = 1.0
 y = 2.0

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
function fix_discrete_variables(var_value::Function, model::GenericModel)
    return _relax_or_fix_integrality(var_value, model)
end

function fix_discrete_variables(model::GenericModel)
    return fix_discrete_variables(value, model)
end

function _relax_or_fix_integrality(
    var_value::Union{Nothing,Function},
    model::GenericModel{T},
) where {T}
    if num_constraints(model, GenericVariableRef{T}, MOI.Semicontinuous{T}) > 0
        error(
            "Support for relaxing semicontinuous constraints is not " *
            "yet implemented.",
        )
    end
    if num_constraints(model, GenericVariableRef{T}, MOI.Semiinteger{T}) > 0
        error(
            "Support for relaxing semi-integer constraints is not " *
            "yet implemented.",
        )
    end
    discrete_variable_constraints = vcat(
        all_constraints(model, GenericVariableRef{T}, MOI.ZeroOne),
        all_constraints(model, GenericVariableRef{T}, MOI.Integer),
    )
    # We gather the info first because we cannot modify-then-query.
    info_pre_relaxation = map(discrete_variable_constraints) do c
        v = GenericVariableRef{T}(c)
        solution = var_value === nothing ? nothing : var_value(v)
        return (v, solution, _info_from_variable(v))
    end
    # Now we can modify.
    for (v, solution, info) in info_pre_relaxation
        if info.integer
            unset_integer(v)
        elseif info.binary
            unset_binary(v)
            if !info.has_fix
                set_lower_bound(v, max(zero(T), info.lower_bound))
                set_upper_bound(v, min(one(T), info.upper_bound))
            elseif info.fixed_value < 0 || info.fixed_value > 1
                error(
                    "The model has no valid relaxation: binary variable " *
                    "fixed out of bounds.",
                )
            end
        end
        if solution !== nothing
            fix(v, solution; force = true)
        end
    end
    function unrelax()
        for (v, solution, info) in info_pre_relaxation
            if solution !== nothing
                unfix(v)
            end
            if info.has_lb
                set_lower_bound(v, info.lower_bound)
            end
            if info.has_ub
                set_upper_bound(v, info.upper_bound)
            end
            if info.integer
                set_integer(v)
            end
            if info.binary
                set_binary(v)
            end
            # Now a special case: when binary variables are relaxed, we add
            # [0, 1] bounds, but only if the variable was not previously fixed
            # and we did not provide a fixed value, and a bound did not already
            # exist. In this case, delete the new bounds that we added.
            if solution === nothing && info.binary && !info.has_fix
                if !info.has_lb
                    delete_lower_bound(v)
                end
                if !info.has_ub
                    delete_upper_bound(v)
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
       https://jump.dev/forum

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
        Base.$(sym)(::GenericVariableRef, ::Number) = error($(msg))
        Base.$(sym)(::Number, ::GenericVariableRef) = error($(msg))
    end
end
