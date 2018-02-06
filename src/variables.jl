#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


#############################################################################
# VariableIndex
# Holds a reference to the model and the corresponding MOI.VariableIndex.
struct Variable <: AbstractJuMPScalar
    m::Model
    index::MOIVAR
end

function MOI.delete!(m::Model, v::Variable)
    @assert m === v.m
    MOI.delete!(m.moibackend, v.index)
end

MOI.isvalid(m::Model, v::Variable) = (v.m === m) && MOI.isvalid(m.moibackend, v.index)


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

function Base.getindex(vm::VariableToValueMap, v::Variable)
    @assert v.m === vm.m # TODO: better error message
    return vm.d[index(v)]
end

function Base.setindex!(vm::VariableToValueMap{T}, value::T, v::Variable) where T
    @assert v.m === vm.m # TODO: better error message
    vm.d[index(v)] = value
end

Base.setindex!(vm::VariableToValueMap{T}, value, v::Variable) where T = setindex!(vm, convert(T, value), v)

function Base.delete!(vm::VariableToValueMap,v::Variable)
    delete!(vm.d, index(v))
    vm
end

Base.empty!(vm::VariableToValueMap) = empty!(vm.d)
Base.isempty(vm::VariableToValueMap) = isempty(vm.d)

Base.haskey(vm::VariableToValueMap, v::Variable) = (vm.m === v.m) && haskey(vm.d, index(v))



index(v::Variable) = v.index

function Variable(m::Model)
    index = MOI.addvariable!(m.moibackend)
    return Variable(m, index)
end

# Name setter/getters

"""
    name(v::Variable)::String

Get a variable's name.
"""
name(v::Variable) = MOI.get(v.m, MOI.VariableName(), v)

"""
    setname(v::Variable,s::AbstractString)

Set a variable's name.
"""
setname(v::Variable, s::String) = MOI.set!(v.m, MOI.VariableName(), v, s)

MOI.SingleVariable(v::Variable) = MOI.SingleVariable(index(v))

# Note: No validation is performed that the variables belong to the same model.
MOI.VectorOfVariables(vars::Vector{Variable}) = MOI.VectorOfVariables(index.(vars))

Variable(m::Model, f::MOI.SingleVariable) = Variable(m, f.variable)

function setobjective(m::Model, sense::Symbol, x::Variable)
    # TODO: This code is repeated here, for AffExpr, and for QuadExpr.
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
    objectivefunction(m::Model, ::Type{Variable})

Return a `Variable` object representing the objective function.
Error if the objective is not a `SingleVariable`.
"""
function objectivefunction(m::Model, ::Type{Variable})
    f = MOI.get(m.moibackend, MOI.ObjectiveFunction{MOI.SingleVariable}())::MOI.SingleVariable
    return Variable(m, f)
end

struct SingleVariableConstraint{S <: MOI.AbstractScalarSet} <: AbstractConstraint
    func::Variable
    set::S
end

moi_function_and_set(c::SingleVariableConstraint) = (MOI.SingleVariable(c.func), c.set)

struct VectorOfVariablesConstraint{S <: MOI.AbstractVectorSet} <: AbstractConstraint
    func::Vector{Variable}
    set::S
end

moi_function_and_set(c::VectorOfVariablesConstraint) = (MOI.VectorOfVariables(c.func), c.set)

function constraintobject(cr::ConstraintRef{Model}, ::Type{Variable}, ::Type{SetType}) where {SetType <: MOI.AbstractScalarSet}
    f = MOI.get(cr.m, MOI.ConstraintFunction(), cr)::MOI.SingleVariable
    s = MOI.get(cr.m, MOI.ConstraintSet(), cr)::SetType
    return SingleVariableConstraint(Variable(cr.m, f), s)
end

function constraintobject(cr::ConstraintRef{Model}, ::Type{Vector{Variable}}, ::Type{SetType}) where {SetType <: MOI.AbstractVectorSet}
    m = cr.m
    f = MOI.get(m, MOI.ConstraintFunction(), cr)::MOI.VectorOfVariables
    s = MOI.get(m, MOI.ConstraintSet(), cr)::SetType
    return VectorOfVariablesConstraint(map(v -> Variable(m, v), f.variables), s)
end


## Bound setter/getters

# lower bounds

haslowerbound(v::Variable) = haskey(v.m.variabletolowerbound,index(v))

function lowerboundindex(v::Variable)
    @assert haslowerbound(v) # TODO error message
    return v.m.variabletolowerbound[index(v)]
end
function setlowerboundindex(v::Variable, cindex::MOILB)
    v.m.variabletolowerbound[index(v)] = cindex
end

"""
    setlowerbound(v::Variable,lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower bound constraint.
"""
function setlowerbound(v::Variable,lower::Number)
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

function LowerBoundRef(v::Variable)
    return ConstraintRef{Model, MOILB}(v.m, lowerboundindex(v))
end

"""
    deletelowerbound(v::Variable)

Delete the lower bound constraint of a variable.
"""
function deletelowerbound(v::Variable)
    MOI.delete!(v.m, LowerBoundRef(v))
    delete!(v.m.variabletolowerbound, index(v))
    nothing
end

"""
    lowerbound(v::Variable)

Return the lower bound of a variable. Error if one does not exist.
"""
function lowerbound(v::Variable)
    cset = MOI.get(v.m, MOI.ConstraintSet(), LowerBoundRef(v))::MOI.GreaterThan
    return cset.lower
end

# upper bounds

hasupperbound(v::Variable) = haskey(v.m.variabletoupperbound,index(v))

function upperboundindex(v::Variable)
    @assert hasupperbound(v) # TODO error message
    return v.m.variabletoupperbound[index(v)]
end
function setupperboundindex(v::Variable, cindex::MOIUB)
    v.m.variabletoupperbound[index(v)] = cindex
end

"""
    setupperbound(v::Variable,upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound constraint.
"""
function setupperbound(v::Variable,upper::Number)
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

function UpperBoundRef(v::Variable)
    return ConstraintRef{Model, MOIUB}(v.m, upperboundindex(v))
end

"""
    deleteupperbound(v::Variable)

Delete the upper bound constraint of a variable.
"""
function deleteupperbound(v::Variable)
    MOI.delete!(v.m, UpperBoundRef(v))
    delete!(v.m.variabletoupperbound, index(v))
    nothing
end

"""
    upperbound(v::Variable)

Return the upper bound of a variable. Error if one does not exist.
"""
function upperbound(v::Variable)
    cset = MOI.get(v.m, MOI.ConstraintSet(), UpperBoundRef(v))::MOI.LessThan
    return cset.upper
end

# fixed value

isfixed(v::Variable) = haskey(v.m.variabletofix,index(v))

function fixindex(v::Variable)
    @assert isfixed(v) # TODO error message
    return v.m.variabletofix[index(v)]
end
function setfixindex(v::Variable, cindex::MOIFIX)
    v.m.variabletofix[index(v)] = cindex
end

"""
    fix(v::Variable,upper::Number)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise create a new one.
"""
function fix(v::Variable,upper::Number)
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
    unfix(v::Variable)

Delete the fixing constraint of a variable.
"""
function unfix(v::Variable)
    MOI.delete!(v.m, FixRef(v))
    delete!(v.m.variabletofix, index(v))
    nothing
end

"""
    fixvalue(v::Variable)

Return the value to which a variable is fixed. Error if one does not exist.
"""
function fixvalue(v::Variable)
    cset = MOI.get(v.m, MOI.ConstraintSet(), FixRef(v))::MOI.EqualTo
    return cset.value
end

function FixRef(v::Variable)
    return ConstraintRef{Model,MOIFIX}(v.m, fixindex(v))
end

# integer and binary constraints

isinteger(v::Variable) = haskey(v.m.variabletointegrality,index(v))

function integerindex(v::Variable)
    @assert isinteger(v) # TODO error message
    return v.m.variabletointegrality[index(v)]
end
function setintegerindex(v::Variable, cindex::MOIINT)
    v.m.variabletointegrality[index(v)] = cindex
end

"""
    setinteger(v::Variable)

Add an integrality constraint on the variable `v`.
"""
function setinteger(v::Variable)
    isinteger(v) && return
    @assert !isbinary(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), MOI.Integer())
    setintegerindex(v, cindex)
    nothing
end

function unsetinteger(v::Variable)
    MOI.delete!(v.m, IntegerRef(v))
    delete!(v.m.variabletointegrality, index(v))
end

function IntegerRef(v::Variable)
    return ConstraintRef{Model,MOIINT}(v.m, integerindex(v))
end

isbinary(v::Variable) = haskey(v.m.variabletozeroone,index(v))

function binaryindex(v::Variable)
    @assert isbinary(v) # TODO error message
    return v.m.variabletozeroone[index(v)]
end
function setbinaryindex(v::Variable, cindex::MOIBIN)
    v.m.variabletozeroone[index(v)] = cindex
end

"""
    setbinary(v::Variable)

Add a constraint on the variable `v` that it must take values in the set ``\{0,1\}``.
"""
function setbinary(v::Variable)
    isbinary(v) && return
    @assert !isinteger(v) # TODO error message
    cindex = MOI.addconstraint!(v.m.moibackend, MOI.SingleVariable(index(v)), MOI.ZeroOne())
    setbinaryindex(v, cindex)
    nothing
end

function unsetbinary(v::Variable)
    MOI.delete!(v.m, BinaryRef(v))
    delete!(v.m.variabletozeroone, index(v))
end

function BinaryRef(v::Variable)
    return ConstraintRef{Model,MOIBIN}(v.m, binaryindex(v))
end


startvalue(v::Variable) = MOI.get(v.m, MOI.VariablePrimalStart())
setstartvalue(v::Variable, val::Number) = MOI.set!(v.m, MOI.VariablePrimalStart(), v, val)

"""
    resultvalue(v::Variable)

Get the value of this variable in the result returned by a solver.
Use `hasresultvalues` to check if a result exists before asking for values.
Replaces `getvalue` for most use cases.
"""
resultvalue(v::Variable) = MOI.get(v.m, MOI.VariablePrimal(), v)
hasresultvalues(m::Model) = MOI.canget(m, MOI.VariablePrimal(), Variable)

@Base.deprecate setvalue(v::Variable, val::Number) setstart(v, val)
