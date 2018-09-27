#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

export addlazycallback, addcutcallback, addheuristiccallback, addinfocallback

abstract type JuMPCallback end
mutable struct LazyCallback <: JuMPCallback
    f::Function
    fractional::Bool
end
mutable struct CutCallback <: JuMPCallback
    f::Function
end
mutable struct HeuristicCallback <: JuMPCallback
    f::Function
end
mutable struct InfoCallback <: JuMPCallback
    f::Function
    when::Symbol
end

mutable struct StopTheSolver end

mutable struct CallbackAbort <: Exception end
export CallbackAbort

Base.copy(c::T) where {T<:JuMPCallback} = T(copy(c))
Base.copy(c::LazyCallback) = LazyCallback(copy(c.f), c.fractional)

function addlazycallback(m::Model, f::Function; fractional::Bool=false)
    m.internalModelLoaded = false
    push!(m.callbacks, LazyCallback(f,fractional))
end
function addcutcallback(m::Model, f::Function)
    m.internalModelLoaded = false
    push!(m.callbacks, CutCallback(f))
end
function addheuristiccallback(m::Model, f::Function)
    m.internalModelLoaded = false
    push!(m.callbacks, HeuristicCallback(f))
end

function _unspecifiedstate()
    Compat.@warn("""Info Callbacks without a when clause will currently default
         to fire only in the :Intermediate state to preserve its behavior
         with previous versions of JuMP.

         This behavior will be deprecated in subsequent versions of JuMP, so
         please rewrite all invocations of addinfocallbacks(m, f) to
         addinfocallbacks(m, f, when = :Intermediate) instead.
         """)
    :Intermediate
end

function addinfocallback(m::Model, f::Function; when::Symbol = _unspecifiedstate())
    m.internalModelLoaded = false
    push!(m.callbacks, InfoCallback(f, when))
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
                ret = cb.f(d)
                if ret === StopTheSolver
                    return :Exit
                end
            end
        end
    catch y
        if isa(y, CallbackAbort)
            @Compat.warn "Throwing CallbackAbort() from a callback is deprecated. Use \"return JuMP.StopTheSolver\" instead."
            return :Exit
        else
            rethrow(y)
        end
    end
    :Continue
end

function attach_callbacks(m::Model, cbs::Vector{LazyCallback})
    cb = d -> lazycallback(d,m,cbs)
    if applicable(MathProgBase.setlazycallback!, m.internalModel, cb)
        MathProgBase.setlazycallback!(m.internalModel, cb)
    else
        error("Solver does not support lazy callbacks")
    end
end

function cutcallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{CutCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPNode
    MathProgBase.cbgetlpsolution(d,m.colVal)
    try
        for cb in cbs
            ret = cb.f(d)
            if ret === StopTheSolver
                return :Exit
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

function attach_callbacks(m::Model, cbs::Vector{CutCallback})
    cb = d -> cutcallback(d,m,cbs)
    if applicable(MathProgBase.setcutcallback!, m.internalModel, cb)
        MathProgBase.setcutcallback!(m.internalModel, cb)
    else
        error("Solver does not support cut callbacks")
    end
end


function heurcallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{HeuristicCallback})
    state = MathProgBase.cbgetstate(d)
    @assert state == :MIPNode
    MathProgBase.cbgetlpsolution(d,m.colVal)
    try
        for cb in cbs
            ret = cb.f(d)
            if ret === StopTheSolver
                return :Exit
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

function attach_callbacks(m::Model, cbs::Vector{HeuristicCallback})
    cb = d -> heurcallback(d,m,cbs)
    if applicable(MathProgBase.setheuristiccallback!, m.internalModel, cb)
        MathProgBase.setheuristiccallback!(m.internalModel, cb)
    else
        error("Solver does not support heuristic callbacks")
    end
end

function infocallback(d::MathProgBase.MathProgCallbackData, m::Model, cbs::Vector{InfoCallback})
    state = MathProgBase.cbgetstate(d)
    if state == :MIPSol
        MathProgBase.cbgetmipsolution(d,m.colVal)
    elseif state == :MIPNode
        MathProgBase.cbgetlpsolution(d,m.colVal)
    end
    try
        for cb in cbs
            if cb.when == state
                ret = cb.f(d)
                if ret === StopTheSolver
                    return :Exit
                end
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

function attach_callbacks(m::Model, cbs::Vector{InfoCallback})
    cb = d -> infocallback(d,m,cbs)
    if applicable(MathProgBase.setinfocallback!, m.internalModel, cb)
        MathProgBase.setinfocallback!(m.internalModel, cb)
    else
        error("Solver does not support info callbacks")
    end
end

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
const sensemap = Dict(:(<=) => '<', :(==) => '=', :(>=) => '>')


## Lazy constraints
export @lazyconstraint

macro lazyconstraint(args...)
    cbdata = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)
    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:(=)), extra) # filtering expressions corresponding to kw args specs
    extra = filter(ex->!isexpr(ex,:(=)), extra) # others

    localcut_val = false # by default, the lazy constraint is global
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :localcut
            localcut_val = esc(ex.args[2])   # excepted if otherwise specified ...
        else error("in @lazyconstraint($(join(args,','))), invalid keyword argument: $(kwarg)")
        end
    end

    if isexpr(x, :call) && length(x.args) == 3 # simple comparison
        lhs = :($(x.args[2]) - $(x.args[3])) # move everything to the lhs
        newaff, parsecode = parseExprToplevel(lhs, :aff)
        sense, vectorized = _canonicalize_sense(x.args[1])
        code = quote
            aff = zero(AffExpr)
            $parsecode
            constr = constructconstraint!($newaff, $(quot(sense)))
        end
        if vectorized
            return quote
                $code
                for con in constr
                    addlazyconstraint($cbdata, con, localcut=$localcut_val)
                end
            end
        else
            return quote
                $code
                addlazyconstraint($cbdata, constr, localcut=$localcut_val)
            end
        end
    else
        error("Syntax error in @lazyconstraint, expected one-sided comparison.")
    end
end

function addlazyconstraint(cbdata::MathProgBase.MathProgCallbackData, constr::LinearConstraint; localcut::Bool=false)
    if length(constr.terms.vars) == 0
        MathProgBase.cbaddlazy!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    assert_isfinite(constr.terms)
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    if localcut
        MathProgBase.cbaddlazylocal!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
    else
        MathProgBase.cbaddlazy!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
    end
end

addlazyconstraint(cbdata::MathProgBase.MathProgCallbackData,constr::QuadConstraint; localcut::Bool=false) = error("Quadratic lazy constraints are not supported.")

## User cuts
export addusercut
export @usercut

macro usercut(args...)
    cbdata = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)
    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:(=)), extra) # filtering expressions corresponding to kw args specs
    extra = filter(ex->!isexpr(ex,:(=)), extra) # others

    localcut_val = false # by default, the user cut is global
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :localcut
            localcut_val = esc(ex.args[2])   # excepted if otherwise specified ...
        else error("in @usercut($(join(args,','))), invalid keyword argument: $(kwarg)")
        end
    end

    #cbdata = esc(cbdata)
    if isexpr(x, :call) && length(x.args) == 3 # simple comparison
        lhs = :($(x.args[2]) - $(x.args[3])) # move everything to the lhs
        newaff, parsecode = parseExprToplevel(lhs, :aff)
        sense, vectorized = _canonicalize_sense(x.args[1])
        vectorized && error("Cannot add vectorized constraint in cut callback")
        quote
            aff = zero(AffExpr)
            $parsecode
            constr = constructconstraint!($newaff, $(quot(sense)))
            addusercut($cbdata, constr, localcut=$localcut_val)
        end
    else
        error("Syntax error in @usercut, expected one-sided comparison.")
    end
end

function addusercut(cbdata::MathProgBase.MathProgCallbackData, constr::LinearConstraint; localcut::Bool=false)
    if length(constr.terms.vars) == 0
        MathProgBase.cbaddcut!(cbdata, Cint[], Float64[], sensemap[sense(constr)], rhs(constr))
        return
    end
    assert_isfinite(constr.terms)
    m::Model = constr.terms.vars[1].m
    indices, coeffs = merge_duplicates(Cint, constr.terms, m.indexedVector, m)
    if localcut
        MathProgBase.cbaddcutlocal!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
    else
        MathProgBase.cbaddcut!(cbdata, indices, coeffs, sensemap[sense(constr)], rhs(constr))
    end
end

## User heuristic
export addsolution, setsolutionvalue

addsolution(cbdata::MathProgBase.MathProgCallbackData) = MathProgBase.cbaddsolution!(cbdata)
function setsolutionvalue(cbdata::MathProgBase.MathProgCallbackData, v::Variable, x)
    MathProgBase.cbsetsolutionvalue!(cbdata, convert(Cint, v.col), x)
end
