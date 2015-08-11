#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

export addLazyCallback, addCutCallback, addHeuristicCallback, addInfoCallback
@Base.deprecate setLazyCallback      addLazyCallback
@Base.deprecate setCutCallback       addCutCallback
@Base.deprecate setHeuristicCallback addHeuristicCallback

abstract JuMPCallback
type LazyCallback <: JuMPCallback
    f::Function
    fractional::Bool
end
type CutCallback <: JuMPCallback
    f::Function
end
type HeuristicCallback <: JuMPCallback
    f::Function
end
type InfoCallback <: JuMPCallback
    f::Function
end

type CallbackAbort <: Exception end
export CallbackAbort

Base.copy{T<:JuMPCallback}(c::T) = T(copy(c))
Base.copy(c::LazyCallback) = LazyCallback(copy(c.f), c.fractional)

function addLazyCallback(m::Model, f::Function; fractional::Bool=false)
    m.internalModelLoaded = false
    push!(m.callbacks, LazyCallback(f,fractional))
end
function addCutCallback(m::Model, f::Function)
    m.internalModelLoaded = false
    push!(m.callbacks, CutCallback(f))
end
function addHeuristicCallback(m::Model, f::Function)
    m.internalModelLoaded = false
    push!(m.callbacks, HeuristicCallback(f))
end
function addInfoCallback(m::Model, f::Function)
    m.internalModelLoaded = false
    push!(m.callbacks, InfoCallback(f))
end

function lazycallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{LazyCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPSol || state == :MIPNode
    anyfrac = mapreduce(|, cbs) do cb
        cb.fractional
    end
    if state == :MIPSol
        MathProgBase.cbgetmipsolution(d,m.colVal)
    elseif anyfrac
        MathProgBase.cbgetlpsolution(d,m.colVal)
    else
        return :Continue
    end
    try
        for cb in cbs
            if state == :MIPSol || cb.fractional
                cb.f(d)
            end
        end
    catch y
        if isa(y, CallbackAbort)
            return :Exit
        else
            rethrow(y)
        end
    end
    :Continue
end

attach_callbacks(m::Model, cbs::Vector{LazyCallback}) =
    MathProgBase.setlazycallback!(m.internalModel, d -> lazycallback(d,m,cbs))

function cutcallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{CutCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPNode
    MathProgBase.cbgetlpsolution(d,m.colVal)
    try
        for cb in cbs
            cb.f(d)
        end
    catch y
        if isa(y, CallbackAbort)
            return :Exit
        else
            rethrow(y)
        end
    end
    :Continue
end

attach_callbacks(m::Model, cbs::Vector{CutCallback}) =
    MathProgBase.setcutcallback!(m.internalModel, d -> cutcallback(d,m,cbs))

function heurcallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{HeuristicCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPNode
    MathProgBase.cbgetlpsolution(d,m.colVal)
    try
        for cb in cbs
            cb.f(d)
        end
    catch y
        if isa(y, CallbackAbort)
            return :Exit
        else
            rethrow(y)
        end
    end
    :Continue
end

attach_callbacks(m::Model, cbs::Vector{HeuristicCallback}) =
    MathProgBase.setheuristiccallback!(m.internalModel, d -> heurcallback(d,m,cbs))

function infocallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{InfoCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPInfo
    try
        for cb in cbs
            cb.f(d)
        end
    catch y
        if isa(y, CallbackAbort)
            return :Exit
        else
            rethrow(y)
        end
    end
    :Continue
end

attach_callbacks(m::Model, cbs::Vector{InfoCallback}) =
    MathProgBase.setinfocallback!(m.internalModel, d -> infocallback(d,m,cbs))

function registercallbacks(m::Model)
    isempty(m.callbacks) && return # might as well avoid allocating the indexedVector

    cbtypes = unique(map(typeof, m.callbacks))
    for typ in cbtypes
        cbs::Vector{typ} = filter(m.callbacks) do cb
            isa(cb, typ)
        end
        attach_callbacks(m, cbs)
    end

    # prepare storage for callbacks
    m.indexedVector = IndexedVector(Float64, m.numCols)
end


# TODO: Should this be somewhere else?
const sensemap = @compat Dict(:(<=) => '<', :(==) => '=', :(>=) => '>')


## Lazy constraints
export addLazyConstraint, @addLazyConstraint

macro addLazyConstraint(cbdata, x)
    cbdata = esc(cbdata)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        newaff, parsecode = parseExprToplevel(lhs, :aff)
        sense, vectorized = _canonicalize_sense(x.args[2])
        vectorized && error("Cannot add vectorized constraint in lazy callback")
        quote
            aff = zero(AffExpr)
            $parsecode
            constr = constructconstraint!($newaff, $(quot(sense)))
            addLazyConstraint($cbdata, constr)
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addLazyConstraint(cbdata::MathProgBase.MathProgCallbackData, constr::LinearConstraint)
    if length(constr.terms.vars) == 0
        MathProgBase.cbaddlazy!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    assert_isfinite(constr.terms)
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    MathProgBase.cbaddlazy!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
end

addLazyConstraint(cbdata::MathProgBase.MathProgCallbackData,constr::QuadConstraint) = error("Quadratic lazy constraints are not supported.")

## User cuts
export addUserCut, @addUserCut

macro addUserCut(cbdata, x)
    cbdata = esc(cbdata)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    if length(x.args) == 3 # simple comparison
        lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
        newaff, parsecode = parseExprToplevel(lhs, :aff)
        sense, vectorized = _canonicalize_sense(x.args[2])
        vectorized && error("Cannot add vectorized constraint in cut callback")
        quote
            aff = zero(AffExpr)
            $parsecode
            constr = constructconstraint!($newaff, $(quot(sense)))
            addUserCut($cbdata, constr)
        end
    else
        error("Syntax error (ranged constraints not permitted in callbacks)")
    end
end

function addUserCut(cbdata::MathProgBase.MathProgCallbackData, constr::LinearConstraint)
    if length(constr.terms.vars) == 0
        MathProgBase.cbaddcut!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    assert_isfinite(constr.terms)
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    MathProgBase.cbaddcut!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
end

## User heuristic
export addSolution, setSolutionValue!

addSolution(cbdata::MathProgBase.MathProgCallbackData) = MathProgBase.cbaddsolution!(cbdata)
function setSolutionValue!(cbdata::MathProgBase.MathProgCallbackData, v::Variable, x)
    MathProgBase.cbsetsolutionvalue!(cbdata, convert(Cint, v.col), x)
end
