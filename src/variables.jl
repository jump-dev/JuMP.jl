#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

function _set_lower_bound_or_error(_error::Function, info::_VariableInfoExpr, lower)
    info.has_lb && _error("Cannot specify variable lower_bound twice")
    info.has_lb = true
    info.lower_bound = lower
end
function _set_upper_bound_or_error(_error::Function, info::_VariableInfoExpr, upper)
    info.has_ub && _error("Cannot specify variable upper_bound twice")
    info.has_ub = true
    info.upper_bound = upper
end
function _fix_or_error(_error::Function, info::_VariableInfoExpr, value)
    info.has_fix && _error("Cannot specify variable fixed value twice")
    info.has_fix = true
    info.fixed_value = value
end
function _set_binary_or_error(_error::Function, info::_VariableInfoExpr)
    info.binary === false || _error("'Bin' and 'binary' keyword argument cannot both be specified.")
    info.binary = true
end
function _set_integer_or_error(_error::Function, info::_VariableInfoExpr)
    info.integer === false || _error("'Int' and 'integer' keyword argument cannot both be specified.")
    info.integer = true
end

function _is_info_keyword(kw::Expr)
    kw.args[1] in [:lower_bound, :upper_bound, :start, :binary, :integer]
end
# :(start = 0)     -> (:start, 0)
# :(start = i + 1) -> (:start, :($(Expr(:escape, :(i + 1)))))
function _keywordify(kw::Expr)
    (kw.args[1], _esc_non_constant(kw.args[2]))
end
function _VariableInfoExpr(; lower_bound=NaN, upper_bound=NaN, start=NaN, binary=false, integer=false)
    # isnan(::Expr) is not defined so we need to do !== NaN
    _VariableInfoExpr(lower_bound !== NaN, lower_bound, upper_bound !== NaN, upper_bound, false, NaN, start !== NaN, start, binary, integer)
end

struct VariableInfo{S, T, U, V}
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
end

function _constructor_expr(info::_VariableInfoExpr)
    return :(VariableInfo($(info.has_lb), $(info.lower_bound), $(info.has_ub),
             $(info.upper_bound), $(info.has_fix), $(info.fixed_value),
             $(info.has_start), $(info.start), $(info.binary), $(info.integer)))
end

struct ScalarVariable{S, T, U, V} <: AbstractVariable
    info::VariableInfo{S, T, U, V}
end

"""
    AbstractVariableRef

Variable returned by [`add_variable`](@ref). Affine (resp. quadratic) operations
with variables of type `V<:AbstractVariableRef` and coefficients of type `T`
    create a `GenericAffExpr{T,V}` (resp. `GenericQuadExpr{T,V}`).
"""
abstract type AbstractVariableRef <: AbstractJuMPScalar end

variable_ref_type(v::AbstractVariableRef) = typeof(v)

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
```jldoctest
julia> model = Model()

julia> x = @variable(model)

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
struct VariableNotOwned{V <: AbstractVariableRef} <: Exception
    variable::V
end

"""
    check_belongs_to_model(func::AbstractJuMPScalar, model::AbstractModel)

Throw `VariableNotOwned` if the `owner_model` of one of the variables of the
function `func` is not `model`.

    check_belongs_to_model(constraint::AbstractConstraint, model::AbstractModel)

Throw `VariableNotOwned` if the `owner_model` of one of the variables of the
constraint `constraint` is not `model`.
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

isequal_canonical(v::VariableRef, other::VariableRef) = isequal(v, other)

"""
    delete(model::Model, variable_ref::VariableRef)

Delete the variable associated with `variable_ref` from the model `model`.
"""
function delete(model::Model, variable_ref::VariableRef)
    if model !== owner_model(variable_ref)
        error("The variable reference you are trying to delete does not " *
              "belong to the model.")
    end
    MOI.delete(backend(model), variable_ref.index)
end

"""
    is_valid(model::Model, variable_ref::VariableRef)

Return `true` if `variable` refers to a valid variable in `model`.
"""
function is_valid(model::Model, variable_ref::VariableRef)
    return (model === owner_model(variable_ref) &&
            MOI.is_valid(backend(model), variable_ref.index))
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/JuliaOpt/MathOptInterface.jl/issues/234#issuecomment-366868878
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

# Name setter/getters
# These functions need to be implemented for all `AbstractVariableRef`s
"""
    name(v::VariableRef)::String

Get a variable's name attribute.
"""
name(v::VariableRef) = MOI.get(owner_model(v), MOI.VariableName(), v)::String

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
julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

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

MOI.SingleVariable(v::VariableRef) = MOI.SingleVariable(index(v))
function moi_function(variable::AbstractVariableRef)
    return MOI.SingleVariable(variable)
end
function moi_function_type(::Type{<:AbstractVariableRef})
    return MOI.SingleVariable
end

# Note: No validation is performed that the variables belong to the same model.
MOI.VectorOfVariables(vars::Vector{VariableRef}) = MOI.VectorOfVariables(index.(vars))
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

function VariableRef(model::Model, f::MOI.SingleVariable)
    return VariableRef(model, f.variable)
end
jump_function_type(::Model, ::Type{MOI.SingleVariable}) = VariableRef
function jump_function(model::Model, variable::MOI.SingleVariable)
    return VariableRef(model, variable)
end

## Bound setter/getters

# lower bounds

"""
    has_lower_bound(v::VariableRef)

Return `true` if `v` has a lower bound. If `true`, the lower bound can be
queried with [`lower_bound`](@ref). See also [`LowerBoundRef`](@ref).
"""
function has_lower_bound(v::VariableRef)
    return _moi_has_lower_bound(backend(owner_model(v)), v)
end

# _moi_* methods allow us to work around the type instability of the backend of
# a model.
function _moi_has_lower_bound(backend, v::VariableRef)
    return MOI.is_valid(backend, _lower_bound_index(v))
end

function _lower_bound_index(v::VariableRef)
    return _MOICON{MOI.SingleVariable, MOI.GreaterThan{Float64}}(index(v).value)
end



"""
    set_lower_bound(v::VariableRef, lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower
bound constraint. See also [`delete_lower_bound`](@ref).
"""
function set_lower_bound(v::VariableRef, lower::Number)
    return _moi_set_lower_bound(backend(owner_model(v)), v, lower)
end

function _moi_set_lower_bound(backend, v::VariableRef, lower::Number)
    new_set = MOI.GreaterThan(convert(Float64, lower))
    if _moi_has_lower_bound(backend, v)
        cindex = _lower_bound_index(v)
        MOI.set(backend, MOI.ConstraintSet(), cindex, new_set)
    else
        @assert !_moi_is_fixed(backend, v)
        MOI.add_constraint(backend, MOI.SingleVariable(index(v)), new_set)
    end
    return
end

"""
    LowerBoundRef(v::VariableRef)

Return a constraint reference to the lower bound constraint of `v`. Errors if
one does not exist.
"""
function LowerBoundRef(v::VariableRef)
    moi_lb = _MOICON{MOI.SingleVariable, MOI.GreaterThan{Float64}}
    return ConstraintRef{Model, moi_lb, ScalarShape}(owner_model(v),
                                                     _lower_bound_index(v),
                                                     ScalarShape())
end

"""
    delete_lower_bound(v::VariableRef)

Delete the lower bound constraint of a variable.
"""
function delete_lower_bound(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), LowerBoundRef(variable_ref))
    return
end

"""
    lower_bound(v::VariableRef)

Return the lower bound of a variable. Error if one does not exist. See also
[`has_lower_bound`](@ref).
"""
function lower_bound(v::VariableRef)
    if !has_lower_bound(v)
        error("Variable $(v) does not have a lower bound.")
    end
    cset = MOI.get(owner_model(v), MOI.ConstraintSet(),
                   LowerBoundRef(v))::MOI.GreaterThan{Float64}
    return cset.lower
end

# upper bounds

"""
    has_upper_bound(v::VariableRef)

Return `true` if `v` has a upper bound. If `true`, the upper bound can be
queried with [`upper_bound`](@ref). See also [`UpperBoundRef`](@ref).
"""
function has_upper_bound(v::VariableRef)
    return _moi_has_upper_bound(backend(owner_model(v)), v)
end

function _moi_has_upper_bound(backend, v::VariableRef)
    return MOI.is_valid(backend, _upper_bound_index(v))
end

function _upper_bound_index(v::VariableRef)
    return _MOICON{MOI.SingleVariable, MOI.LessThan{Float64}}(index(v).value)
end

"""
    set_upper_bound(v::VariableRef,upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound
constraint. See also [`delete_upper_bound`](@ref).
"""
function set_upper_bound(v::VariableRef, upper::Number)
    return _moi_set_upper_bound(backend(owner_model(v)), v, upper)
end

function _moi_set_upper_bound(backend, v::VariableRef, upper::Number)
    new_set = MOI.LessThan(convert(Float64,upper))
    if _moi_has_upper_bound(backend, v)
        cindex = _upper_bound_index(v)
        MOI.set(backend, MOI.ConstraintSet(), cindex, new_set)
    else
        @assert !_moi_is_fixed(backend, v)
        MOI.add_constraint(backend, MOI.SingleVariable(index(v)), new_set)
    end
    return
end

"""
    UpperBoundRef(v::VariableRef)

Return a constraint reference to the upper bound constraint of `v`. Errors if
one does not exist.
"""
function UpperBoundRef(v::VariableRef)
    moi_ub = _MOICON{MOI.SingleVariable, MOI.LessThan{Float64}}
    return ConstraintRef{Model, moi_ub, ScalarShape}(owner_model(v),
                                                     _upper_bound_index(v),
                                                     ScalarShape())
end

"""
    delete_upper_bound(v::VariableRef)

Delete the upper bound constraint of a variable.
"""
function delete_upper_bound(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), UpperBoundRef(variable_ref))
    return
end

"""
    upper_bound(v::VariableRef)

Return the upper bound of a variable. Error if one does not exist. See also
[`has_upper_bound`](@ref).
"""
function upper_bound(v::VariableRef)
    if !has_upper_bound(v)
        error("Variable $(v) does not have an upper bound.")
    end
    cset = MOI.get(owner_model(v), MOI.ConstraintSet(),
                   UpperBoundRef(v))::MOI.LessThan{Float64}
    return cset.upper
end

# fixed value

"""
    is_fixed(v::VariableRef)

Return `true` if `v` is a fixed variable. If `true`, the fixed value can be
queried with [`fix_value`](@ref). See also [`FixRef`](@ref).
"""
function is_fixed(v::VariableRef)
    return _moi_is_fixed(backend(owner_model(v)), v)
end

function _moi_is_fixed(backend, v::VariableRef)
    return MOI.is_valid(backend, _fix_index(v))
end

function _fix_index(v::VariableRef)
    return _MOICON{MOI.SingleVariable, MOI.EqualTo{Float64}}(index(v).value)
end

"""
    fix(v::VariableRef, value::Number; force::Bool = false)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise
create a new one. See also [`unfix`](@ref).

If the variable already has variable bounds and `force=false`, calling `fix`
will throw an error. If `force=true`, existing variable bounds will be deleted,
and the fixing constraint will be added. Note a variable will have no bounds
after a call to [`unfix`](@ref).
"""
function fix(variable::VariableRef, value::Number; force::Bool = false)
    return _moi_fix(backend(owner_model(variable)), variable, value, force)
end

function _moi_fix(backend, variable::VariableRef, value::Number, force::Bool)
    new_set = MOI.EqualTo(convert(Float64, value))
    if _moi_is_fixed(backend, variable)  # Update existing fixing constraint.
        c_index = _fix_index(variable)
        MOI.set(backend, MOI.ConstraintSet(), c_index, new_set)
    else  # Add a new fixing constraint.
        if _moi_has_upper_bound(backend, variable) ||
            _moi_has_lower_bound(backend, variable)
            if !force
                error("Unable to fix $(variable) to $(value) because it has " *
                      "existing variable bounds. Consider calling " *
                      "`JuMP.fix(variable, value; force=true)` which will " *
                      "delete existing bounds before fixing the variable.")
            end
            if _moi_has_upper_bound(backend, variable)
                delete_upper_bound(variable)
            end
            if _moi_has_lower_bound(backend, variable)
                delete_lower_bound(variable)
            end
        end
        MOI.add_constraint(backend, MOI.SingleVariable(index(variable)),
                           new_set)
    end
    return
end

"""
    unfix(v::VariableRef)

Delete the fixing constraint of a variable.
"""
function unfix(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), FixRef(variable_ref))
    return
end

"""
    fix_value(v::VariableRef)

Return the value to which a variable is fixed. Error if one does not exist. See
also [`is_fixed`](@ref).
"""
function fix_value(v::VariableRef)
    cset = MOI.get(owner_model(v), MOI.ConstraintSet(),
                   FixRef(v))::MOI.EqualTo{Float64}
    return cset.value
end

"""
    FixRef(v::VariableRef)

Return a constraint reference to the constraint fixing the value of `v`. Errors
if one does not exist.
"""
function FixRef(v::VariableRef)
    moi_fix = _MOICON{MOI.SingleVariable, MOI.EqualTo{Float64}}
    return ConstraintRef{Model, moi_fix, ScalarShape}(owner_model(v),
                                                      _fix_index(v),
                                                      ScalarShape())
end

"""
    is_integer(v::VariableRef)

Return `true` if `v` is constrained to be integer. See also
[`IntegerRef`](@ref).
"""
function is_integer(v::VariableRef)
    return _moi_is_integer(backend(owner_model(v)), v)
end

function _moi_is_integer(backend, v::VariableRef)
    return MOI.is_valid(backend, _integer_index(v))
end

function _integer_index(v::VariableRef)
    return _MOICON{MOI.SingleVariable, MOI.Integer}(index(v).value)
end

"""
    set_integer(variable_ref::VariableRef)

Add an integrality constraint on the variable `variable_ref`. See also
[`unset_integer`](@ref).
"""
function set_integer(variable_ref::VariableRef)
    return _moi_set_integer(backend(owner_model(variable_ref)), variable_ref)
end

function _moi_set_integer(backend, variable_ref::VariableRef)
    if _moi_is_integer(backend, variable_ref)
        return
    elseif _moi_is_binary(backend, variable_ref)
        error("Cannot set the variable_ref $(variable_ref) to integer as it " *
              "is already binary.")
    end
    MOI.add_constraint(backend, MOI.SingleVariable(index(variable_ref)),
                       MOI.Integer())
    return
end

"""
    unset_integer(variable_ref::VariableRef)

Remove the integrality constraint on the variable `variable_ref`.
"""
function unset_integer(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), IntegerRef(variable_ref))
    return
end

"""
    IntegerRef(v::VariableRef)

Return a constraint reference to the constraint constrainting `v` to be integer.
Errors if one does not exist.
"""
function IntegerRef(v::VariableRef)
    moi_int = _MOICON{MOI.SingleVariable, MOI.Integer}
    return ConstraintRef{Model, moi_int, ScalarShape}(
        owner_model(v), _integer_index(v), ScalarShape())
end

"""
    is_binary(v::VariableRef)

Return `true` if `v` is constrained to be binary. See also [`BinaryRef`](@ref).
"""
function is_binary(v::VariableRef)
    return _moi_is_binary(backend(owner_model(v)), v)
end

function _moi_is_binary(backend, v::VariableRef)
    return MOI.is_valid(backend, _binary_index(v))
end

function _binary_index(v::VariableRef)
    return _MOICON{MOI.SingleVariable, MOI.ZeroOne}(index(v).value)
end

"""
    set_binary(v::VariableRef)

Add a constraint on the variable `v` that it must take values in the set
``\\{0,1\\}``. See also [`unset_binary`](@ref).
"""
function set_binary(variable_ref::VariableRef)
    return _moi_set_binary(backend(owner_model(variable_ref)), variable_ref)
end

function _moi_set_binary(backend, variable_ref)
    if _moi_is_binary(backend, variable_ref)
        return
    elseif _moi_is_integer(backend, variable_ref)
        error("Cannot set the variable_ref $(variable_ref) to binary as it " *
              "is already integer.")
    end
    MOI.add_constraint(backend, MOI.SingleVariable(index(variable_ref)),
                       MOI.ZeroOne())
    return
end

"""
    unset_binary(variable_ref::VariableRef)

Remove the binary constraint on the variable `variable_ref`.
"""
function unset_binary(variable_ref::VariableRef)
    JuMP.delete(owner_model(variable_ref), BinaryRef(variable_ref))
    return
end

"""
    BinaryRef(v::VariableRef)

Return a constraint reference to the constraint constrainting `v` to be binary.
Errors if one does not exist.
"""
function BinaryRef(v::VariableRef)
    moi_bin = _MOICON{MOI.SingleVariable, MOI.ZeroOne}
    return ConstraintRef{Model, moi_bin, ScalarShape}(
        owner_model(v), _binary_index(v), ScalarShape())
end

"""
    start_value(v::VariableRef)

Return the start value (MOI attribute `VariablePrimalStart`) of the variable
`v`. See also [`set_start_value`](@ref).

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".
"""
function start_value(v::VariableRef)::Union{Nothing, Float64}
    return MOI.get(owner_model(v), MOI.VariablePrimalStart(), v)
end


"""
    set_start_value(variable::VariableRef, value::Number)

Set the start value (MOI attribute `VariablePrimalStart`) of the variable `v` to
`value`. See also [`start_value`](@ref).

Note: `VariablePrimalStart`s are sometimes called "MIP-starts" or "warmstarts".
"""
function set_start_value(variable::VariableRef, value::Number)
    MOI.set(owner_model(variable), MOI.VariablePrimalStart(), variable,
            Float64(value))
    return
end

"""
    value(v::VariableRef)

Get the value of this variable in the result returned by a solver. Use
[`has_values`](@ref) to check if a result exists before asking for values.
"""
function value(v::VariableRef)::Float64
    return MOI.get(owner_model(v), MOI.VariablePrimal(), v)
end

"""
    has_values(model::Model)

Return `true` if the solver has a primal solution available to query, otherwise
return `false`. See also [`value`](@ref).
"""
has_values(model::Model) = primal_status(model) != MOI.NO_SOLUTION

@Base.deprecate setvalue(v::VariableRef, val::Number) set_start_value(v, val)

"""
    add_variable(m::Model, v::AbstractVariable, name::String="")

Add a variable `v` to `Model m` and sets its name.
"""
function add_variable end

function add_variable(model::Model, v::ScalarVariable, name::String="")
    return _moi_add_variable(backend(model), model, v, name)
end

function _moi_add_variable(backend, model, v::ScalarVariable, name::String)
    index = MOI.add_variable(backend)
    var_ref = VariableRef(model, index)
    _moi_constrain_variable(backend, index, v.info)
    if !isempty(name)
        set_name(var_ref, name)
    end
    return var_ref
end

function _moi_constrain_variable(backend::MOI.ModelLike, index, info)
    # We don't call the _moi* versions (e.g., _moi_set_lower_bound) because they
    # have extra checks that are not necessary for newly created variables.
    if info.has_lb
        MOI.add_constraint(backend, MOI.SingleVariable(index),
                           MOI.GreaterThan{Float64}(info.lower_bound))
    end
    if info.has_ub
        MOI.add_constraint(backend, MOI.SingleVariable(index),
                           MOI.LessThan{Float64}(info.upper_bound))
    end
    if info.has_fix
        MOI.add_constraint(backend, MOI.SingleVariable(index),
                           MOI.EqualTo{Float64}(info.fixed_value))
    end
    if info.binary
        MOI.add_constraint(backend, MOI.SingleVariable(index),
                           MOI.ZeroOne())
    end
    if info.integer
        MOI.add_constraint(backend, MOI.SingleVariable(index), MOI.Integer())
    end
    if info.has_start
        MOI.set(backend, MOI.VariablePrimalStart(), index,
                Float64(info.start))
    end
end

"""
    ConstrainedVariables <: AbstractVariable

Vector of variables `scalar_variables` constrained to belong to `set`.
Adding this variable can be thought as doing:
```julia
function JuMP.add_variable(model::Model, variable::JuMP.ConstrainedVariables, names)
    var_refs = JuMP.add_variable.(model, variable.scalar_variables,
                                  JuMP.vectorize(names, variable.shape))
    JuMP.add_constraint(model, JuMP.VectorConstraint(var_refs, variable.set))
    return JuMP.reshape_vector(var_refs, variable.shape)
end
```
but adds the variables with `MOI.add_constrained_variables(model, variable.set)`
instead. See [the MOI documentation](http://www.juliaopt.org/MathOptInterface.jl/v0.9.3/apireference/#Variables-1)
for the difference between adding the variables with `MOI.add_constrained_variables`
and adding them with `MOI.add_variables` and adding the constraint separately.
"""
struct ConstrainedVariables{S <: MOI.AbstractVectorSet, Shape <: AbstractShape,
                            ScalarVarType <: AbstractVariable} <: AbstractVariable
    scalar_variables::Vector{ScalarVarType}
    set::S
    shape::Shape
end

function ConstrainedVariables(variables::Vector{<:AbstractVariable}, set::MOI.AbstractVectorSet)
    return ConstrainedVariables(variables, set, VectorShape())
end

function add_variable(model::Model, variable::ConstrainedVariables, names)
    var_indices = _moi_add_constrained_variables(
        backend(model), variable.scalar_variables, variable.set, vectorize(names, variable.shape))
    var_refs = [VariableRef(model, var_index) for var_index in var_indices]
    return reshape_vector(var_refs, variable.shape)
end

function _moi_add_constrained_variables(
    backend::MOI.ModelLike, scalar_variables::Vector{<:ScalarVariable},
    set::MOI.AbstractVectorSet, names::Vector{String})
    if set isa MOI.Reals
        var_indices = MOI.add_variables(backend, MOI.dimension(set))
    else
        var_indices, con_index = MOI.add_constrained_variables(backend, set)
    end
    for (index, variable) in zip(var_indices, scalar_variables)
        _moi_constrain_variable(backend, index, variable.info)
    end
    if names !== nothing
        for (var_index, name) in zip(var_indices, names)
            MOI.set(backend, MOI.VariableName(), var_index, name)
        end
    end
    return var_indices
end

"""
    all_variables(model::Model)::Vector{VariableRef}

Returns a list of all variables currently in the model. The variables are
ordered by creation time.

# Example
```jldoctest
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
    all_indices = MOI.get(
        model, MOI.ListOfVariableIndices())::Vector{MOI.VariableIndex}
    return VariableRef[VariableRef(model, idx) for idx in all_indices]
end

function dual(vref::VariableRef)
    error("To query the dual variables associated with a variable bound, first " *
          "obtain a constraint reference using one of `UpperBoundRef`, `LowerBoundRef`, " *
          "or `FixRef`, and then call `dual` on the returned constraint reference.\nFor " *
          "example, if `x <= 1`, instead of `dual(x)`, call `dual(UpperBoundRef(x))`.")
end

function value(::AbstractArray{<:AbstractJuMPScalar})
    error("`JuMP.value` is not defined for collections of JuMP types. Use" *
          " Julia's broadcast syntax instead: `JuMP.value.(x)`.")
end
