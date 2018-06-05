#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

abstract type AbstractVariable end

# Any fields can usually be either a number or an expression
mutable struct VariableInfoExpr
    haslb::Bool
    lowerbound::Any
    hasub::Bool
    upperbound::Any
    hasfix::Bool
    fixedvalue::Any
    hasstart::Bool
    start::Any
    binary::Any
    integer::Any
end

function setlowerbound(info::VariableInfoExpr, lower, _error::Function)
    info.haslb && _error("Cannot specify variable lowerbound twice")
    info.haslb = true
    info.lowerbound = lower
end
function setupperbound(info::VariableInfoExpr, upper, _error::Function)
    info.hasub && _error("Cannot specify variable lowerbound twice")
    info.hasub = true
    info.upperbound = upper
end
function fix(info::VariableInfoExpr, value, _error::Function)
    info.hasfix && _error("Cannot specify variable fixed value twice")
    info.hasfix = true
    info.fixedvalue = value
end
function setbinary(info::VariableInfoExpr, _error::Function)
    info.binary === false || _error("'Bin' and 'binary' keyword argument cannot both be specified.")
    info.binary = true
end
function setinteger(info::VariableInfoExpr, _error::Function)
    info.integer === false || _error("'Int' and 'integer' keyword argument cannot both be specified.")
    info.integer = true
end

function isinfokeyword(kw::Expr)
    kw.args[1] in [:lowerbound, :upperbound, :start, :binary, :integer]
end
# :(start = 0)     -> (:start, 0)
# :(start = i + 1) -> (:start, :($(Expr(:escape, :(i + 1)))))
function keywordify(kw::Expr)
    (kw.args[1], esc_nonconstant(kw.args[2]))
end
function VariableInfoExpr(; lowerbound=NaN, upperbound=NaN, start=NaN, binary=false, integer=false)
    # isnan(::Expr) is not defined so we need to do !== NaN
    VariableInfoExpr(lowerbound !== NaN, lowerbound, upperbound !== NaN, upperbound, false, NaN, start !== NaN, start, binary, integer)
end

mutable struct VariableInfo{S, T, U, V}
    haslb::Bool
    lowerbound::S
    hasub::Bool
    upperbound::T
    hasfix::Bool
    fixedvalue::U
    hasstart::Bool
    start::V
    binary::Bool
    integer::Bool
end

constructor_expr(info::VariableInfoExpr) = :(VariableInfo($(info.haslb), $(info.lowerbound), $(info.hasub), $(info.upperbound), $(info.hasfix), $(info.fixedvalue), $(info.hasstart), $(info.start), $(info.binary), $(info.integer)))

struct ScalarVariable{S, T, U, V} <: AbstractVariable
    info::VariableInfo{S, T, U, V}
end

"""
    AbstractVariableRef

Variable returned by [`addvariable`](@ref). Affine (resp. quadratic) operations with variables of type `V<:AbstractVariableRef` and coefficients of type `T` create a `GenericAffExpr{T,V}` (resp. `GenericQuadExpr{T,V}`).
"""
abstract type AbstractVariableRef <: AbstractJuMPScalar end

variablereftype(v::AbstractVariableRef) = typeof(v)

"""
    VariableRef <: AbstractVariableRef

Holds a reference to the model and the corresponding MOI.VariableIndex.
"""
struct VariableRef <: AbstractVariableRef
    m::Model
    index::MOIVAR
end

Base.iszero(::VariableRef) = false
Base.copy(v::VariableRef) = VariableRef(v.m, v.index)

isequal_canonical(v::VariableRef, other::VariableRef) = isequal(v, other)

function MOI.delete!(m::Model, v::VariableRef)
    @assert m === v.m
    MOI.delete!(m.moibackend, v.index)
end

MOI.isvalid(m::Model, v::VariableRef) = (v.m === m) && MOI.isvalid(m.moibackend, v.index)

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/JuliaOpt/MathOptInterface.jl/issues/234#issuecomment-366868878
Base.hash(v::VariableRef, h::UInt) = hash(object_id(v.m), hash(v.index.value, h))
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
    index = MOI.addvariable!(m.moibackend)
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
    setname(v::VariableRef,s::AbstractString)

Set a variable's name.
"""
setname(v::VariableRef, s::String) = MOI.set!(v.m, MOI.VariableName(), v, s)

MOI.SingleVariable(v::VariableRef) = MOI.SingleVariable(index(v))

# Note: No validation is performed that the variables belong to the same model.
MOI.VectorOfVariables(vars::Vector{VariableRef}) = MOI.VectorOfVariables(index.(vars))

VariableRef(m::Model, f::MOI.SingleVariable) = VariableRef(m, f.variable)

function setobjective(m::Model, sense::Symbol, x::VariableRef)
    # TODO: This code is repeated here, for GenericAffExpr, and for GenericQuadExpr.
    if sense == :Min
        moisense = MOI.MinSense
    else
        @assert sense == :Max
        moisense = MOI.MaxSense
    end
    MOI.set!(m.moibackend, MOI.ObjectiveSense(), moisense)
    MOI.set!(m.moibackend, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))
end

"""
    objectivefunction(m::Model, ::Type{VariableRef})

Return a `VariableRef` object representing the objective function.
Error if the objective is not a `SingleVariable`.
"""
function objectivefunction(m::Model, ::Type{VariableRef})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.SingleVariable}())::MOI.SingleVariable
    return VariableRef(m, f)
end

struct SingleVariableConstraint{V <: AbstractVariableRef, S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::V
    set::S
end

moi_function_and_set(c::SingleVariableConstraint) = (MOI.SingleVariable(c.func), c.set)

struct VectorOfVariablesConstraint{V <: AbstractVariableRef, S <: MOI.AbstractVectorSet} <: AbstractConstraint
    func::Vector{V}
    set::S
end

moi_function_and_set(c::VectorOfVariablesConstraint) = (MOI.VectorOfVariables(c.func), c.set)

function constraintobject(cr::ConstraintRef{Model}, ::Type{VariableRef}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    f = MOI.get(cr.m, MOI.ConstraintFunction(), cr)::MOI.SingleVariable
    s = MOI.get(cr.m, MOI.ConstraintSet(), cr)::SetType
    return SingleVariableConstraint(VariableRef(cr.m, f), s)
end

function constraintobject(cr::ConstraintRef{Model}, ::Type{Vector{VariableRef}}, ::Type{SetType}) where {SetType <: MOI.AbstractVectorSet}
    m = cr.m
    f = MOI.get(m, MOI.ConstraintFunction(), cr)::MOI.VectorOfVariables
    s = MOI.get(m, MOI.ConstraintSet(), cr)::SetType
    return VectorOfVariablesConstraint(map(v -> VariableRef(m, v), f.variables), s)
end


## Bound setter/getters

# lower bounds

haslowerbound(v::VariableRef) = haskey(v.m.variabletolowerbound,index(v))

function lowerboundindex(v::VariableRef)
    @assert haslowerbound(v) # TODO error message
    return v.m.variabletolowerbound[index(v)]
end
function setlowerboundindex(v::VariableRef, cindex::MOILB)
    v.m.variabletolowerbound[index(v)] = cindex
end

"""
    setlowerbound(v::VariableRef,lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower bound constraint.
"""
function setlowerbound(v::VariableRef,lower::Number)
    newset = MOI.GreaterThan(convert(Float64,lower))
    # do we have a lower bound already?
    if haslowerbound(v)
        cindex = lowerboundindex(v)
        MOI.modifyconstraint!(v.m.moibackend, cindex, newset)
    else
        @assert !isfixed(v)
        cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), newset)
        setlowerboundindex(v, cindex)
    end
    nothing
end

function LowerBoundRef(v::VariableRef)
    return ConstraintRef{Model, MOILB}(v.m, lowerboundindex(v))
end

"""
    deletelowerbound(v::VariableRef)

Delete the lower bound constraint of a variable.
"""
function deletelowerbound(v::VariableRef)
    MOI.delete!(v.m, LowerBoundRef(v))
    delete!(v.m.variabletolowerbound, index(v))
    nothing
end

"""
    lowerbound(v::VariableRef)

Return the lower bound of a variable. Error if one does not exist.
"""
function lowerbound(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), LowerBoundRef(v))::MOI.GreaterThan
    return cset.lower
end

# upper bounds

hasupperbound(v::VariableRef) = haskey(v.m.variabletoupperbound,index(v))

function upperboundindex(v::VariableRef)
    @assert hasupperbound(v) # TODO error message
    return v.m.variabletoupperbound[index(v)]
end
function setupperboundindex(v::VariableRef, cindex::MOIUB)
    v.m.variabletoupperbound[index(v)] = cindex
end

"""
    setupperbound(v::VariableRef,upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound constraint.
"""
function setupperbound(v::VariableRef,upper::Number)
    newset = MOI.LessThan(convert(Float64,upper))
    # do we have an upper bound already?
    if hasupperbound(v)
        cindex = upperboundindex(v)
        MOI.modifyconstraint!(v.m.moibackend, cindex, newset)
    else
        @assert !isfixed(v)
        cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), newset)
        setupperboundindex(v, cindex)
    end
    nothing
end

function UpperBoundRef(v::VariableRef)
    return ConstraintRef{Model, MOIUB}(v.m, upperboundindex(v))
end

"""
    deleteupperbound(v::VariableRef)

Delete the upper bound constraint of a variable.
"""
function deleteupperbound(v::VariableRef)
    MOI.delete!(v.m, UpperBoundRef(v))
    delete!(v.m.variabletoupperbound, index(v))
    nothing
end

"""
    upperbound(v::VariableRef)

Return the upper bound of a variable. Error if one does not exist.
"""
function upperbound(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), UpperBoundRef(v))::MOI.LessThan
    return cset.upper
end

# fixed value

isfixed(v::VariableRef) = haskey(v.m.variabletofix,index(v))

function fixindex(v::VariableRef)
    @assert isfixed(v) # TODO error message
    return v.m.variabletofix[index(v)]
end
function setfixindex(v::VariableRef, cindex::MOIFIX)
    v.m.variabletofix[index(v)] = cindex
end

"""
    fix(v::VariableRef,upper::Number)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise create a new one.
"""
function fix(v::VariableRef,upper::Number)
    newset = MOI.EqualTo(convert(Float64,upper))
    # are we already fixed?
    if isfixed(v)
        cindex = fixindex(v)
        MOI.modifyconstraint!(v.m.moibackend, cindex, newset)
    else
        @assert !hasupperbound(v) && !haslowerbound(v) # Do we want to remove these instead of throwing an error?
        cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), newset)
        setfixindex(v, cindex)
    end
    nothing
end

"""
    unfix(v::VariableRef)

Delete the fixing constraint of a variable.
"""
function unfix(v::VariableRef)
    MOI.delete!(v.m, FixRef(v))
    delete!(v.m.variabletofix, index(v))
    nothing
end

"""
    fixvalue(v::VariableRef)

Return the value to which a variable is fixed. Error if one does not exist.
"""
function fixvalue(v::VariableRef)
    cset = MOI.get(v.m, MOI.ConstraintSet(), FixRef(v))::MOI.EqualTo
    return cset.value
end

function FixRef(v::VariableRef)
    return ConstraintRef{Model,MOIFIX}(v.m, fixindex(v))
end

# integer and binary constraints

isinteger(v::VariableRef) = haskey(v.m.variabletointegrality,index(v))

function integerindex(v::VariableRef)
    @assert isinteger(v) # TODO error message
    return v.m.variabletointegrality[index(v)]
end
function setintegerindex(v::VariableRef, cindex::MOIINT)
    v.m.variabletointegrality[index(v)] = cindex
end

"""
    setinteger(v::VariableRef)

Add an integrality constraint on the variable `v`.
"""
function setinteger(v::VariableRef)
    isinteger(v) && return
    @assert !isbinary(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), MOI.Integer())
    setintegerindex(v, cindex)
    nothing
end

function unsetinteger(v::VariableRef)
    MOI.delete!(v.m, IntegerRef(v))
    delete!(v.m.variabletointegrality, index(v))
end

function IntegerRef(v::VariableRef)
    return ConstraintRef{Model,MOIINT}(v.m, integerindex(v))
end

isbinary(v::VariableRef) = haskey(v.m.variabletozeroone,index(v))

function binaryindex(v::VariableRef)
    @assert isbinary(v) # TODO error message
    return v.m.variabletozeroone[index(v)]
end
function setbinaryindex(v::VariableRef, cindex::MOIBIN)
    v.m.variabletozeroone[index(v)] = cindex
end

"""
    setbinary(v::VariableRef)

Add a constraint on the variable `v` that it must take values in the set ``\\{0,1\\}``.
"""
function setbinary(v::VariableRef)
    isbinary(v) && return
    @assert !isinteger(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), MOI.ZeroOne())
    setbinaryindex(v, cindex)
    nothing
end

function unsetbinary(v::VariableRef)
    MOI.delete!(v.m, BinaryRef(v))
    delete!(v.m.variabletozeroone, index(v))
end

function BinaryRef(v::VariableRef)
    return ConstraintRef{Model,MOIBIN}(v.m, binaryindex(v))
end


startvalue(v::VariableRef) = MOI.get(v.m, MOI.VariablePrimalStart(), v)
setstartvalue(v::VariableRef, val::Number) = MOI.set!(v.m, MOI.VariablePrimalStart(), v, val)

"""
    resultvalue(v::VariableRef)

Get the value of this variable in the result returned by a solver.
Use `hasresultvalues` to check if a result exists before asking for values.
Replaces `getvalue` for most use cases.
"""
resultvalue(v::VariableRef) = MOI.get(v.m, MOI.VariablePrimal(), v)
hasresultvalues(m::Model) = MOI.canget(m, MOI.VariablePrimal(), VariableRef)

@Base.deprecate setvalue(v::VariableRef, val::Number) setstartvalue(v, val)

"""
    addvariable(m::Model, v::AbstractVariable, name::String="")

Add a variable `v` to `Model m` and sets its name.
"""
function addvariable end

function addvariable(m::Model, v::ScalarVariable, name::String="")
    info = v.info
    vref = VariableRef(m)
    if info.haslb
        setlowerbound(vref, info.lowerbound)
    end
    if info.hasub
        setupperbound(vref, info.upperbound)
    end
    if info.hasfix
        fix(vref, info.fixedvalue)
    end
    if info.binary
        setbinary(vref)
    end
    if info.integer
        setinteger(vref)
    end
    if info.hasstart
        setstartvalue(vref, info.start)
    end
    if !isempty(name)
        setname(vref, name)
    end
    return vref
end
