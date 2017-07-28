#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


#############################################################################
# Variable class
# Holds a reference to the model and the corresponding
# VariableReference for the internal instance
struct Variable <: AbstractJuMPScalar
    m::Model
    instanceref::MOIVAR
end


"""
    VariableToValueMap{T}

An object for storing a mapping from variables to a value of type `T`.
"""
struct VariableToValueMap{T}
    m::Model
    d::Dict{MOIVAR,T}
end

function (::Type{VariableToValueMap{T}})(m::Model) where T
    return VariableToValueMap{T}(m, Dict{MOIVAR,T}())
end

function Base.getindex(vm::VariableToValueMap, v::Variable)
    @assert v.m === vm.m # TODO: better error message
    return vm.d[instancereference(v)]
end

function Base.setindex!(vm::VariableToValueMap{T}, value::T, v::Variable) where T
    @assert v.m === vm.m # TODO: better error message
    vm.d[instancereference(v)] = value
end

Base.setindex!(vm::VariableToValueMap{T}, value, v::Variable) where T = setindex!(vm, convert(T, value), v)

Base.empty!(vm::VariableToValueMap) = empty!(vm.d)
Base.isempty(vm::VariableToValueMap) = isempty(vm.d)

Base.haskey(vm::VariableToValueMap, v::Variable) = (vm.m === v.m) && haskey(vm.d, instancereference(v))



instancereference(v::Variable) = v.instanceref

# linearindex(x::Variable) = x.col

# Variable(m::Model, lower, upper, cat::Symbol, name::AbstractString="", value::Number=NaN) =
#     error("Attempt to create scalar Variable with lower bound of type $(typeof(lower)) and upper bound of type $(typeof(upper)). Bounds must be scalars in Variable constructor.")

function Variable(m::Model,name::AbstractString="")
    ref = MOI.addvariable!(m.instance)

    v = Variable(m, ref)

    if name != ""
        m.variablenames[v] = name
    end

    # TODO: try to update solver instance
    # if m.internalModelLoaded
    #     if method_exists(MathProgBase.addvar!, (typeof(m.internalModel),Vector{Int},Vector{Float64},Float64,Float64,Float64))
    #         MathProgBase.addvar!(m.internalModel,float(lower),float(upper),0.0)
    #     else
    #         Base.warn_once("Solver does not appear to support adding variables to an existing model. JuMP's internal model will be discarded.")
    #         m.internalModelLoaded = false
    #     end
    # end
    return Variable(m, ref)
end

# Name setter/getters
# """
#     setname(v::Variable,n::AbstractString)
#
# Set the variable's internal name.
# """
# function setname(v::Variable,n::AbstractString)
#     push!(v.m.customNames, v)
#     v.m.colNames[v.col] = n
#     v.m.colNamesIJulia[v.col] = n
# end

"""
    name(v::Variable)

Get a variable's internal name.
"""
name(v::Variable) = var_str(REPLMode, v)


## Bound setter/getters

# lower bounds

haslowerbound(v::Variable) = haskey(v.m.variabletolowerbound,instancereference(v))

function lowerboundreference(v::Variable)
    @assert haslowerbound(v) # TODO error message
    return v.m.variabletolowerbound[instancereference(v)]
end
function setlowerboundreference(v::Variable, cref::LBREF)
    v.m.variabletolowerbound[instancereference(v)] = cref
end

"""
    setlowerbound(v::Variable,lower::Number)

Set the lower bound of a variable. If one does not exist, create a new lower bound constraint.
"""
function setlowerbound(v::Variable,lower::Number)
    newset = MOI.GreaterThan(convert(Float64,lower))
    # do we have a lower bound already?
    if haslowerbound(v)
        cref = lowerboundreference(v)
        MOI.modifyconstraint!(v.m.instance, cref, newset)
        @assert !v.m.solverinstanceattached # TODO
    else
        @assert !isfixed(v)
        cref = MOI.addconstraint!(v.m.instance, MOI.SingleVariable(instancereference(v)), newset)
        setlowerboundreference(v, cref)
    end
    nothing
end

"""
    deletelowerbound(v::Variable)

Delete the lower bound constraint of a variable.
"""
function deletelowerbound(v::Variable)
    cref = lowerboundreference(v)
    delete!(v.m.instance, cref)
    delete!(v.m.variabletolowerbound, instancereference(v))
    @assert !v.m.solverinstanceattached # TODO
    nothing
end

"""
    lowerbound(v::Variable)

Return the lower bound of a variable. Error if one does not exist.
"""
function lowerbound(v::Variable)
    cref = lowerboundreference(v)
    cset = MOI.getattribute(v.m.instance, MOI.ConstraintSet(), cref)::MOI.GreaterThan
    return cset.lower
end

# upper bounds

hasupperbound(v::Variable) = haskey(v.m.variabletoupperbound,instancereference(v))

function upperboundreference(v::Variable)
    @assert hasupperbound(v) # TODO error message
    return v.m.variabletoupperbound[instancereference(v)]
end
function setupperboundreference(v::Variable, cref::UBREF)
    v.m.variabletoupperbound[instancereference(v)] = cref
end

"""
    setupperbound(v::Variable,upper::Number)

Set the upper bound of a variable. If one does not exist, create an upper bound constraint.
"""
function setupperbound(v::Variable,upper::Number)
    newset = MOI.LessThan(convert(Float64,upper))
    # do we have an upper bound already?
    if hasupperbound(v)
        cref = upperboundreference(v)
        MOI.modifyconstraint!(v.m.instance, cref, newset)
        @assert !v.m.solverinstanceattached # TODO
    else
        @assert !isfixed(v)
        cref = MOI.addconstraint!(v.m.instance, MOI.SingleVariable(instancereference(v)), newset)
        setupperboundreference(v, cref)
    end
    nothing
end

"""
    deleteupperbound(v::Variable)

Delete the upper bound constraint of a variable.
"""
function deleteupperbound(v::Variable)
    cref = upperboundreference(v)
    delete!(v.m.instance, cref)
    delete!(v.m.variabletoupperbound, instancereference(v))
    @assert !v.m.solverinstanceattached # TODO
    nothing
end

"""
    upperbound(v::Variable)

Return the upper bound of a variable. Error if one does not exist.
"""
function upperbound(v::Variable)
    cref = upperboundreference(v)
    cset = MOI.getattribute(v.m.instance, MOI.ConstraintSet(), cref)::MOI.LessThan
    return cset.upper
end

# fixed value

isfixed(v::Variable) = haskey(v.m.variabletofix,instancereference(v))

function fixreference(v::Variable)
    @assert isfixed(v) # TODO error message
    return v.m.variabletofix[instancereference(v)]
end
function setfixreference(v::Variable, cref::FIXREF)
    v.m.variabletofix[instancereference(v)] = cref
end

"""
    fix(v::Variable,upper::Number)

Fix a variable to a value. Update the fixing constraint if one exists, otherwise create a new one.
"""
function fix(v::Variable,upper::Number)
    newset = MOI.EqualTo(convert(Float64,upper))
    # are we already fixed?
    if isfixed(v)
        cref = fixreference(v)
        MOI.modifyconstraint!(v.m.instance, cref, newset)
        @assert !v.m.solverinstanceattached # TODO
    else
        @assert !hasupperbound(v) && !haslowerbound(v) # Do we want to remove these instead of throwing an error?
        cref = MOI.addconstraint!(v.m.instance, MOI.SingleVariable(instancereference(v)), newset)
        setfixreference(v, cref)
    end
    nothing
end

"""
    unfix(v::Variable)

Delete the fixing constraint of a variable.
"""
function unfix(v::Variable)
    cref = getfixreference(v)
    delete!(v.m.instance, cref)
    delete!(v.m.variabletofix, instancereference(v))
    @assert !v.m.solverinstanceattached # TODO
    nothing
end

"""
    fixvalue(v::Variable)

Return the value to which a variable is fixed. Error if one does not exist.
"""
function fixvalue(v::Variable)
    cref = fixreference(v)
    cset = MOI.getattribute(v.m.instance, MOI.ConstraintSet(), cref)::MOI.EqualTo
    return cset.value
end

# integer and binary constraints

isinteger(v::Variable) = haskey(v.m.variabletointegrality,instancereference(v))

function integerreference(v::Variable)
    @assert isinteger(v) # TODO error message
    return v.m.variabletointegrality[instancereference(v)]
end
function setintegerreference(v::Variable, cref::INTREF)
    v.m.variabletointegrality[instancereference(v)] = cref
end

"""
    setinteger(v::Variable)

Add an integrality constraint on the variable `v`.
"""
function setinteger(v::Variable)
    isinteger(v) && return
    @assert !isbinary(v) # TODO error message
    cref = MOI.addconstraint!(v.m.instance, MOI.SingleVariable(instancereference(v)), MOI.Integer())
    setintegerreference(v, cref)
    nothing
end

function unsetinteger(v::Variable)
    cref = integerreference(v)
    delete!(v.m.instance, cref)
    delete!(v.m.variabletointegrality, instancereference(v))
    @assert !v.m.solverinstanceattached # TODO
end

isbinary(v::Variable) = haskey(v.m.variabletozeroone,instancereference(v))

function binaryreference(v::Variable)
    @assert isbinary(v) # TODO error message
    return v.m.variabletozeroone[instancereference(v)]
end
function setbinaryreference(v::Variable, cref::BINREF)
    v.m.variabletozeroone[instancereference(v)] = cref
end

"""
    setbinary(v::Variable)

Add a constraint on the variable `v` that it must take values in the set ``\{0,1\\}``.
"""
function setbinary(v::Variable)
    isbinary(v) && return
    @assert !isinteger(v) # TODO error message
    cref = MOI.addconstraint!(v.m.instance, MOI.SingleVariable(instancereference(v)), MOI.ZeroOne())
    setbinaryreference(v, cref)
    nothing
end

function unsetbinary(v::Variable)
    cref = binaryreference(v)
    delete!(v.m.instance, cref)
    delete!(v.m.variabletozeroone, instancereference(v))
    @assert !v.m.solverinstanceattached # TODO
end


# solution objects

variablestart(m::Model) = m.variablestart::VariableToValueMap{Float64}
# TODO do we want these or should we have people use the variablestart object directly?
startvalue(v::Variable) = variablestart(v.m)[v]
setstartvalue(v::Variable, val::Number) = variablestart(v.m)[v] = val

variableresult(m::Model) = m.variableresult::VariableToValueMap{Float64}
hasvariableresult(m::Model) = !isempty(variableresult(m))

"""
    resultvalue(v::Variable)

Get the value of this variable in the result returned by a solver.
Use `hasvariableresult` to check if a result exists before asking for values.
Replaces `getvalue` for most use cases.
"""
resultvalue(v::Variable) = variableresult(v.m)[v]

@Base.deprecate setvalue(v::Variable, val::Number) setstart(v, val)



# function resultvalue(arr::Array{Variable})
#     ret = similar(arr, Float64)
#     # return immediately for empty array
#     if isempty(ret)
#         return ret
#     end
#     m = first(arr).m
#     # whether this was constructed via @variable, essentially
#     registered = haskey(m.varData, arr)
#     for I in eachindex(arr)
#         ret[I] = resultvalue(arr[I])
#     end
#     # Copy printing data from @variable for Array{Variable} to corresponding Array{Float64} of values
#     if registered
#         m.varData[ret] = m.varData[arr]
#     end
#     ret
# end

# Dual value (reduced cost) getter
#
# # internal method that doesn't print a warning if the value is NaN
# _getDual(v::Variable) = v.m.redCosts[v.col]
#
# getdualwarn(::Variable) = warn("Variable bound duals (reduced costs) not available. Check that the model was properly solved and no integer variables are present.")
#
# """
#     getdual(v::Variable)
#
# Get the reduced cost of this variable in the solution. Similar behavior to `getvalue` for indexable variables.
# """
# function getdual(v::Variable)
#     if length(v.m.redCosts) < MathProgBase.numvar(v.m)
#         getdualwarn(v)
#         NaN
#     else
#         _getDual(v)
#     end
# end

# const var_cats = [:Cont, :Int, :Bin, :SemiCont, :SemiInt]
#
# """
#     setcategory(v::Variable, cat::Symbol)
#
# Set the variable category for `v` after construction. Possible categories are `:Cont, :Int, :Bin, :SemiCont, :SemiInt`.
# """
# function setcategory(v::Variable, cat::Symbol)
#     cat in var_cats || error("Unrecognized variable category $cat. Should be one of:\n    $var_cats")
#     v.m.colCat[v.col] = cat
# end
#
# """
#     getcategory(v::Variable)
#
# Get the variable category for `v`
# """
# getcategory(v::Variable) = v.m.colCat[v.col]
