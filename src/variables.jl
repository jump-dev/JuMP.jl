#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

abstract type AbstractVariable end

# Any fields can usually be either a number or an expression
mutable struct VariableInfoExpr
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

function set_lower_bound_or_error(_error::Function, info::VariableInfoExpr, lower)
    info.has_lb && _error("Cannot specify variable lower_bound twice")
    info.has_lb = true
    info.lower_bound = lower
end
function set_upper_bound_or_error(_error::Function, info::VariableInfoExpr, upper)
    info.has_ub && _error("Cannot specify variable upper_bound twice")
    info.has_ub = true
    info.upper_bound = upper
end
function fix_or_error(_error::Function, info::VariableInfoExpr, value)
    info.has_fix && _error("Cannot specify variable fixed value twice")
    info.has_fix = true
    info.fixed_value = value
end
function set_binary_or_error(_error::Function, info::VariableInfoExpr)
    info.binary === false || _error("'Bin' and 'binary' keyword argument cannot both be specified.")
    info.binary = true
end
function set_integer_or_error(_error::Function, info::VariableInfoExpr)
    info.integer === false || _error("'Int' and 'integer' keyword argument cannot both be specified.")
    info.integer = true
end

function is_info_keyword(kw::Expr)
    kw.args[1] in [:lower_bound, :upper_bound, :start, :binary, :integer]
end
# :(start = 0)     -> (:start, 0)
# :(start = i + 1) -> (:start, :($(Expr(:escape, :(i + 1)))))
function keywordify(kw::Expr)
    (kw.args[1], esc_nonconstant(kw.args[2]))
end
function VariableInfoExpr(; lower_bound=NaN, upper_bound=NaN, start=NaN, binary=false, integer=false)
    # isnan(::Expr) is not defined so we need to do !== NaN
    VariableInfoExpr(lower_bound !== NaN, lower_bound, upper_bound !== NaN, upper_bound, false, NaN, start !== NaN, start, binary, integer)
end

mutable struct VariableInfo{S, T, U, V}
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

function constructor_expr(info::VariableInfoExpr)
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
    m::Model
    index::MOIVAR
end

# Should be implemented by all `AbstractVariableRef`
owner_model(v::VariableRef) = v.m

Base.iszero(::VariableRef) = false
Base.copy(v::VariableRef) = VariableRef(v.m, v.index)
if VERSION >= v"0.7-"
    Base.broadcastable(v::VariableRef) = Ref(v)
end

isequal_canonical(v::VariableRef, other::VariableRef) = isequal(v, other)

"""
    delete(model::Model, variable_ref::VariableRef)

Delete the variable associated with `variable_ref` from the model `model`.
"""
function delete(model::Model, variable_ref::VariableRef)
    @assert model === variable_ref.m
    MOI.delete!(model.moi_backend, variable_ref.index)
end

"""
    is_valid(model::Model, variable_ref::VariableRef)

Return `true` if `variable` refers to a valid variable in `model`.
"""
function is_valid(model::Model, variable_ref::VariableRef)
    return (model === variable_ref.m &&
            MOI.isvalid(model.moi_backend, variable_ref.index))
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/JuliaOpt/MathOptInterface.jl/issues/234#issuecomment-366868878
Base.hash(v::VariableRef, h::UInt) = hash(objectid(v.m), hash(v.index.value, h))
Base.isequal(v1::VariableRef, v2::VariableRef) = v1.m === v2.m && v1.index == v2.index


"""
    VariableToValueMap{T}

An object for storing a mapping from variables to a value of type `T`.
"""
struct VariableToValueMap{T}
    m::Model
    d::Dict{MOIVAR,T}
end

function VariableToValueMap{T}(m::Model) where T
    return VariableToValueMap{T}(m, Dict{MOIVAR,T}())
end

function Base.getindex(vm::VariableToValueMap, v::VariableRef)
    @assert v.m === vm.m # TODO: better error message
    return vm.d[index(v)]
end

function Base.setindex!(vm::VariableToValueMap{T}, value::T, v::VariableRef) where T
    @assert v.m === vm.m # TODO: better error message
    vm.d[index(v)] = value
end

Base.setindex!(vm::VariableToValueMap{T}, value, v::VariableRef) where T = setindex!(vm, convert(T, value), v)

function Base.delete!(vm::VariableToValueMap,v::VariableRef)
    delete!(vm.d, index(v))
    vm
end

Base.empty!(vm::VariableToValueMap) = empty!(vm.d)
Base.isempty(vm::VariableToValueMap) = isempty(vm.d)

Base.haskey(vm::VariableToValueMap, v::VariableRef) = (vm.m === v.m) && haskey(vm.d, index(v))



index(v::VariableRef) = v.index

function VariableRef(m::Model)
    index = MOI.addvariable!(m.moi_backend)
    return VariableRef(m, index)
end

# Name setter/getters
# These functions need to be implemented for all `AbstractVariableRef`s
"""
    name(v::VariableRef)::String

Get a variable's name.
"""
name(v::VariableRef) = MOI.get(v.m, MOI.VariableName(), v)

"""
    set_name(v::VariableRef,s::AbstractString)

Set a variable's name.
"""
set_name(v::VariableRef, s::String) = MOI.set!(v.m, MOI.VariableName(), v, s)

MOI.SingleVariable(v::VariableRef) = MOI.SingleVariable(index(v))

# Note: No validation is performed that the variables belong to the same model.
MOI.VectorOfVariables(vars::Vector{VariableRef}) = MOI.VectorOfVariables(index.(vars))

VariableRef(m::Model, f::MOI.SingleVariable) = VariableRef(m, f.variable)

function set_objective(m::Model, sense::Symbol, x::VariableRef)
    # TODO: This code is repeated here, for GenericAffExpr, and for GenericQuadExpr.
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.moi_backend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.moi_backend, MOI.ObjectiveFunction{MOI.SingleVariable}(),
             MOI.SingleVariable(x))
end

"""
    objective_function(m::Model, ::Type{VariableRef})

Return a `VariableRef` object representing the objective function.
Error if the objective is not a `SingleVariable`.
"""
function objective_function(m::Model, ::Type{VariableRef})
    f = MOI.get(m.moi_backend, MOI.ObjectiveFunction{MOI.SingleVariable}())::MOI.SingleVariable
    return VariableRef(m, f)
end

struct SingleVariableConstraint{V <: AbstractVariableRef,
                                S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::V
    set::S
end

moi_function_and_set(c::SingleVariableConstraint) = (MOI.SingleVariable(c.func), c.set)
shape(::SingleVariableConstraint) = ScalarShape()

struct VectorOfVariablesConstraint{V <: AbstractVariableRef, S <: MOI.AbstractVectorSet, Shape <: AbstractShape} <: AbstractConstraint
    func::Vector{V}
    set::S
    shape::Shape
end
function VectorOfVariablesConstraint(func::Vector{<:AbstractVariableRef},
                                     set::MOI.AbstractVectorSet)
    VectorOfVariablesConstraint(func, set, VectorShape())
end

moi_function_and_set(c::VectorOfVariablesConstraint) = (MOI.VectorOfVariables(c.func), c.set)
shape(c::VectorOfVariablesConstraint) = c.shape

function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.SingleVariable, SetType <: MOI.AbstractScalarSet}
    model = ref.m
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return SingleVariableConstraint(VariableRef(model, f), s)
end

function constraint_object(ref::ConstraintRef{Model, MOICON{FuncType, SetType}}) where
        {FuncType <: MOI.VectorOfVariables, SetType <: MOI.AbstractVectorSet}
    model = ref.m
    f = MOI.get(model, MOI.ConstraintFunction(), ref)::FuncType
    s = MOI.get(model, MOI.ConstraintSet(), ref)::SetType
    return VectorOfVariablesConstraint(map(v -> VariableRef(model, v), f.variables),
                                       s, ref.shape)
end


## Bound setter/getters

# lower bounds

has_lower_bound(v::VariableRef) = haskey(v.m.variable_to_lower_bound,index(v))

function lower_bound_index(v::VariableRef)
    @assert has_lower_bound(v) # TODO error message
    return v.m.variable_to_lower_bound[index(v)]
end
function set_lower_bound_index(v::VariableRef, cindex::MOILB)
    v.m.variable_to_lower_bound[index(v)] = cindex
end

"""
    set_lower_bound(v::VariableRef,lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower bound constraint.
"""
function set_lower_bound(v::VariableRef,lower::Number)
    newset = MOI.GreaterThan(convert(Float64,lower))
    # do we have a lower bound already?
    if has_lower_bound(v)
        cindex = lower_bound_index(v)
        MOI.set!(v.m.moi_backend, MOI.ConstraintSet(), cindex, newset)
    else
        @assert !is_fixed(v)
        cindex = MOI.addconstraint!(v.m.moi_backend, MOI.SingleVariable(index(v)), newset)
        set_lower_bound_index(v, cindex)
    end
    nothing
end

function LowerBoundRef(v::VariableRef)
    return ConstraintRef{Model, MOILB, ScalarShape}(v.m, lower_bound_index(v),
                                                    ScalarShape())
end

"""
    delete_lower_bound(v::VariableRef)

Delete the lower bound constraint of a variable.
"""
function delete_lower_bound(variable_ref::VariableRef)
    JuMP.delete(variable_ref.m, LowerBoundRef(variable_ref))
    delete!(variable_ref.m.variable_to_lower_bound, index(variable_ref))
    return nothing
end

"""
    lower_bound(v::VariableRef)

Return the lower bound of a variable. Error if one does not exist.
"""
function lower_bound(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), LowerBoundRef(v))::MOI.GreaterThan
    return cset.lower
end

# upper bounds

has_upper_bound(v::VariableRef) = haskey(v.m.variable_to_upper_bound,index(v))

function upper_bound_index(v::VariableRef)
    @assert has_upper_bound(v) # TODO error message
    return v.m.variable_to_upper_bound[index(v)]
end
function set_upper_bound_index(v::VariableRef, cindex::MOIUB)
    v.m.variable_to_upper_bound[index(v)] = cindex
end

"""
    set_upper_bound(v::VariableRef,upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound constraint.
"""
function set_upper_bound(v::VariableRef,upper::Number)
    newset = MOI.LessThan(convert(Float64,upper))
    # do we have an upper bound already?
    if has_upper_bound(v)
        cindex = upper_bound_index(v)
        MOI.set!(v.m.moi_backend, MOI.ConstraintSet(), cindex, newset)
    else
        @assert !is_fixed(v)
        cindex = MOI.addconstraint!(v.m.moi_backend, MOI.SingleVariable(index(v)), newset)
        set_upper_bound_index(v, cindex)
    end
    nothing
end

function UpperBoundRef(v::VariableRef)
    return ConstraintRef{Model, MOIUB, ScalarShape}(v.m, upper_bound_index(v),
                                                    ScalarShape())
end

"""
    delete_upper_bound(v::VariableRef)

Delete the upper bound constraint of a variable.
"""
function delete_upper_bound(variable_ref::VariableRef)
    JuMP.delete(variable_ref.m, UpperBoundRef(variable_ref))
    delete!(variable_ref.m.variable_to_upper_bound, index(variable_ref))
    return nothing
end

"""
    upper_bound(v::VariableRef)

Return the upper bound of a variable. Error if one does not exist.
"""
function upper_bound(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), UpperBoundRef(v))::MOI.LessThan
    return cset.upper
end

# fixed value

is_fixed(v::VariableRef) = haskey(v.m.variable_to_fix,index(v))

function fix_index(v::VariableRef)
    @assert is_fixed(v) # TODO error message
    return v.m.variable_to_fix[index(v)]
end
function set_fix_index(v::VariableRef, cindex::MOIFIX)
    v.m.variable_to_fix[index(v)] = cindex
end

"""
    fix(v::VariableRef,upper::Number)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise create a new one.
"""
function fix(v::VariableRef,upper::Number)
    newset = MOI.EqualTo(convert(Float64,upper))
    # are we already fixed?
    if is_fixed(v)
        cindex = fix_index(v)
        MOI.set!(v.m.moi_backend, MOI.ConstraintSet(), cindex, newset)
    else
        @assert !has_upper_bound(v) && !has_lower_bound(v) # Do we want to remove these instead of throwing an error?
        cindex = MOI.addconstraint!(v.m.moi_backend, MOI.SingleVariable(index(v)), newset)
        set_fix_index(v, cindex)
    end
    nothing
end

"""
    unfix(v::VariableRef)

Delete the fixing constraint of a variable.
"""
function unfix(variable_ref::VariableRef)
    JuMP.delete(variable_ref.m, FixRef(variable_ref))
    delete!(variable_ref.m.variable_to_fix, index(variable_ref))
    return nothing
end

"""
    fix_value(v::VariableRef)

Return the value to which a variable is fixed. Error if one does not exist.
"""
function fix_value(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), FixRef(v))::MOI.EqualTo
    return cset.value
end

function FixRef(v::VariableRef)
    return ConstraintRef{Model, MOIFIX, ScalarShape}(v.m, fix_index(v),
                                                     ScalarShape())
end

# integer and binary constraints

is_integer(v::VariableRef) = haskey(v.m.variable_to_integrality,index(v))

function integer_index(v::VariableRef)
    @assert is_integer(v) # TODO error message
    return v.m.variable_to_integrality[index(v)]
end
function set_integer_index(v::VariableRef, cindex::MOIINT)
    v.m.variable_to_integrality[index(v)] = cindex
end

"""
    set_integer(v::VariableRef)

Add an integrality constraint on the variable `v`.
"""
function set_integer(v::VariableRef)
    is_integer(v) && return
    @assert !is_binary(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moi_backend, MOI.SingleVariable(index(v)), MOI.Integer())
    set_integer_index(v, cindex)
    nothing
end

function unset_integer(variable_ref::VariableRef)
    JuMP.delete(variable_ref.m, IntegerRef(variable_ref))
    delete!(variable_ref.m.variable_to_integrality, index(variable_ref))
    return nothing
end

function IntegerRef(v::VariableRef)
    return ConstraintRef{Model, MOIINT, ScalarShape}(v.m, integer_index(v),
                                                     ScalarShape())
end

is_binary(v::VariableRef) = haskey(v.m.variable_to_zero_one,index(v))

function binary_index(v::VariableRef)
    @assert is_binary(v) # TODO error message
    return v.m.variable_to_zero_one[index(v)]
end
function set_binary_index(v::VariableRef, cindex::MOIBIN)
    v.m.variable_to_zero_one[index(v)] = cindex
end

"""
    set_binary(v::VariableRef)

Add a constraint on the variable `v` that it must take values in the set ``\\{0,1\\}``.
"""
function set_binary(v::VariableRef)
    is_binary(v) && return
    @assert !is_integer(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moi_backend, MOI.SingleVariable(index(v)), MOI.ZeroOne())
    set_binary_index(v, cindex)
    nothing
end

function unset_binary(variable_ref::VariableRef)
    JuMP.delete(variable_ref.m, BinaryRef(variable_ref))
    delete!(variable_ref.m.variable_to_zero_one, index(variable_ref))
    return nothing
end

function BinaryRef(v::VariableRef)
    return ConstraintRef{Model, MOIBIN, ScalarShape}(v.m, binary_index(v),
                                                     ScalarShape())
end


start_value(v::VariableRef) = MOI.get(v.m, MOI.VariablePrimalStart(), v)
set_start_value(v::VariableRef, val::Number) = MOI.set!(v.m, MOI.VariablePrimalStart(), v, val)

"""
    result_value(v::VariableRef)

Get the value of this variable in the result returned by a solver.
Use `has_result_values` to check if a result exists before asking for values.
Replaces `getvalue` for most use cases.
"""
result_value(v::VariableRef) = MOI.get(v.m, MOI.VariablePrimal(), v)
has_result_values(m::Model) = MOI.canget(m, MOI.VariablePrimal(), VariableRef)

@Base.deprecate setvalue(v::VariableRef, val::Number) set_start_value(v, val)

"""
    add_variable(m::Model, v::AbstractVariable, name::String="")

Add a variable `v` to `Model m` and sets its name.
"""
function add_variable end

function add_variable(m::Model, v::ScalarVariable, name::String="")
    info = v.info
    vref = VariableRef(m)
    if info.has_lb
        set_lower_bound(vref, info.lower_bound)
    end
    if info.has_ub
        set_upper_bound(vref, info.upper_bound)
    end
    if info.has_fix
        fix(vref, info.fixed_value)
    end
    if info.binary
        set_binary(vref)
    end
    if info.integer
        set_integer(vref)
    end
    if info.has_start
        set_start_value(vref, info.start)
    end
    if !isempty(name)
        set_name(vref, name)
    end
    return vref
end
