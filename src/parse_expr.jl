#  Copyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function _parse_idx_set(arg::Expr)
    parse_done, idxvar, idxset = Containers._try_parse_idx_set(arg)
    if parse_done
        return idxvar, idxset
    end
    error("Invalid syntax: $arg")
end

"""
    destructive_add!(expression, terms...)

Returns the result of `expression + (*)(terms...)` and possibly destroys
`expression` (i.e., reusing its storage for the result, if possible). For
example `result = destructive_add!(expression, 2, x)` sets `result` to the value
of `expression + 2*x`. `expression` should not be referred to after this call.

Note: Discarding the result of `destructive_add!` suggests an improper usage.
"""
function destructive_add! end

destructive_add!(ex::Constant, c::Constant, x::Constant) = ex + c*x

destructive_add!(ex::Constant, c::Constant, x::AbstractVariableRef) = GenericAffExpr(_float(ex), x => _float(c))

function destructive_add!(ex::Constant, c::Constant, x::T) where T<:GenericAffExpr
    # It's only safe to mutate the first argument.
    if iszero(c)
        convert(T, ex)
    else
        x = c*x
        x.constant += ex
        x
    end
end

function destructive_add!(ex::Constant, c::Constant, x::T) where T<:GenericQuadExpr
    # It's only safe to mutate the first argument.
    if iszero(c)
        convert(T, ex)
    else
        result = c*x
        result.aff.constant += ex
        result
    end
end

function destructive_add!(ex::Constant, c::AbstractVariableRef, x::AbstractVariableRef)
    result = c*x
    result.aff.constant += ex
    result
end

function destructive_add!(ex::Constant, c::T, x::T) where T<:GenericAffExpr
    q = c*x
    q.aff.constant += ex
    q
end

function destructive_add!(ex::Constant, c::GenericAffExpr{C,V}, x::V) where {C,V}
    q = c*x
    q.aff.constant += ex
    q
end

function destructive_add!(ex::Constant, c::T, x::Constant) where T<:GenericQuadExpr
    if iszero(x)
        convert(T, ex)
    else
        q = c*x
        q.aff.constant += ex
        q
    end
end

function destructive_add!(ex::AbstractVariableRef, c::Constant, x::AbstractVariableRef)
    return GenericAffExpr(0.0, ex => 1.0, x => _float(c))
end
function destructive_add!(ex::AbstractVariableRef, c::Constant, x::Constant)
    return GenericAffExpr(_float(c * x), ex => 1.0)
end

function destructive_add!(aff::GenericAffExpr, c::Constant, x::Constant)
    aff.constant += c*x
    aff
end

function destructive_add!(aff::GenericAffExpr{C,V}, c::Constant, x::V) where {C,V}
    add_to_expression!(aff, convert(C, c), x)
    aff
end

function destructive_add!(aff::GenericAffExpr{C,V},c::Constant,x::GenericAffExpr{C,V}) where {C,V}
    if !iszero(c)
        sizehint!(aff, length(linear_terms(aff)) + length(linear_terms(x)))
        for (coef, var) in linear_terms(x)
            add_to_expression!(aff, c*coef, var)
        end
        aff.constant += c*x.constant
    end
    aff
end

# help w/ ambiguity
destructive_add!(aff::GenericAffExpr{C,V}, c::Constant, x::Constant) where {C,V<:Constant} = aff + c*x

function destructive_add!(aff::GenericAffExpr{C,V}, c::V, x::Constant) where {C,V}
    if !iszero(x)
        add_to_expression!(aff, convert(C, x), c)
    end
    aff
end

destructive_add!(aff::GenericAffExpr{C,V},c::V,x::V) where {C,V} =
    GenericQuadExpr{C,V}(aff, UnorderedPair(c,x) => 1.0)

# TODO: add generic versions of following two methods
function destructive_add!(aff::GenericAffExpr, c::GenericAffExpr, x::AbstractVariableRef)
    quad = c*x
    quad.aff = destructive_add!(quad.aff, 1.0, aff)
    quad
end

function destructive_add!(aff::GenericAffExpr, c::AbstractVariableRef, x::GenericAffExpr)
    quad = c*x
    # TODO: Consider implementing the add_to_expression! method for cases like
    # this and below.
    quad.aff = destructive_add!(quad.aff, 1.0, aff)
    quad
end

function destructive_add!(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Constant) where {C,V}
    destructive_add!(aff, x, c)
end

function destructive_add!(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Constant) where {C,V}
    if iszero(x)
        GenericQuadExpr{C,V}(aff)
    else
        result = c*x
        add_to_expression!(result, aff)
        return result
    end
end

function destructive_add!(aff::GenericAffExpr{C,V}, c::Constant, x::GenericQuadExpr{C,V}) where {C,V}
    if iszero(c)
        GenericQuadExpr{C,V}(aff)
    else
        result = c*x
        add_to_expression!(result, aff)
        return result
    end
end

function destructive_add!(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}
    q = convert(GenericQuadExpr{C,V}, ex)
    destructive_add!(q, one(C), c*x)
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::Constant,x::V) where {C,V}
    if !iszero(c)
        add_to_expression!(quad.aff, convert(C, c), x)
    end
    quad
end

function destructive_add!(quad::GenericQuadExpr,c::Constant,x::Constant)
    quad.aff.constant += c*x
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},x::V,y::V) where {C,V}
    add_to_expression!(quad, one(C), x, y)
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::Constant,x::GenericAffExpr{C,V}) where {C,V}
    quad.aff = destructive_add!(quad.aff, c, x)
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Constant) where {C,V}
    quad.aff = destructive_add!(quad.aff,c,x)
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::V) where {C,V}
    for (coef, var) in linear_terms(c)
        add_to_expression!(quad, coef, var, x)
    end
    if !iszero(c.constant)
        quad.aff = destructive_add!(quad.aff,c.constant,x)
    end
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V}) where {C,V}
    for (coef, var) in linear_terms(x)
        add_to_expression!(quad, coef, c, var)
    end
    if !iszero(x.constant)
        quad.aff = destructive_add!(quad.aff,c,x.constant)
    end
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Constant) where {C,V}
    if !iszero(x)
        for (coef, var1, var2) in quad_terms(c)
            add_to_expression!(quad, x*coef, var1, var2)
        end
        quad.aff = destructive_add!(quad.aff,c.aff,x)
    end
    quad
end

function destructive_add!(quad::GenericQuadExpr{C,V},c::Constant,x::GenericQuadExpr{C,V}) where {C,V}
    if !iszero(c)
        for (coef, var1, var2) in quad_terms(x)
            add_to_expression!(quad, c*coef, var1, var2)
        end
        quad.aff = destructive_add!(quad.aff,c,x.aff)
    end
    quad
end

function destructive_add!(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V}) where {C,V}
    destructive_add!(ex, 1.0, c*x)
end

# Catch nonlinear expressions and parameters being used in add_constraint, etc.

const _NLExpr = Union{NonlinearExpression,NonlinearParameter}
_nl_expr_err() = error("""Cannot use nonlinear expression or parameter in @constraint or @objective.
                          Use @NLconstraint or @NLobjective instead.""")
# Following three definitions avoid ambiguity warnings
destructive_add!(expr::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::V) where {C,V<:_NLExpr} = _nl_expr_err()
destructive_add!(expr::GenericQuadExpr{C,V}, c::Constant, x::V) where {C,V<:_NLExpr} = _nl_expr_err()
destructive_add!(expr::GenericQuadExpr{C,V}, c::V, x::GenericAffExpr{C,V}) where {C,V<:_NLExpr} = _nl_expr_err()
for T1 in (GenericAffExpr,GenericQuadExpr), T2 in (Constant,VariableRef,GenericAffExpr,GenericQuadExpr)
    @eval destructive_add!(::$T1, ::$T2, ::_NLExpr) = _nl_expr_err()
    @eval destructive_add!(::$T1, ::_NLExpr, ::$T2) = _nl_expr_err()
end

function destructive_add!(ex::AbstractArray{T}, c::AbstractArray, x::AbstractArray) where {T<:GenericAffExpr}
    add_to_expression!.(ex, c*x)
end
function destructive_add!(ex::AbstractArray{T}, c::AbstractArray, x::Constant) where {T<:GenericAffExpr}
    add_to_expression!.(ex, c*x)
end
function destructive_add!(ex::AbstractArray{T}, c::Constant, x::AbstractArray) where {T<:GenericAffExpr}
    add_to_expression!.(ex, c*x)
end
function destructive_add!(ex::AbstractArray{T}, c::Number, x::Number) where {T<:GenericAffExpr}
    add_to_expression!.(ex, c*x)
end

function destructive_add!(ex::AbstractArray{<:AbstractJuMPScalar}, c::Number, x::UniformScaling)
    return ex + c*x
end
function destructive_add!(ex::AbstractArray{<:AbstractJuMPScalar}, c::UniformScaling, x::Number)
    return ex + c*x
end
function destructive_add!(ex::AbstractArray{<:AbstractJuMPScalar}, c::UniformScaling, x::UniformScaling)
    return ex + c*x
end
function destructive_add!(c::Number, x::UniformScaling, ex::AbstractArray{<:AbstractJuMPScalar})
    add_to_expression!.(x * ex, c)
end
function destructive_add!(c::UniformScaling, x::Constant, ex::AbstractArray{<:AbstractJuMPScalar})
    return c + x * ex
end

function destructive_add!(ex::AbstractArray{<:GenericAffExpr},
                          c::AbstractArray{<:GenericQuadExpr},
                          x::Constant)
    result = c*x
    add_to_expression!.(result, ex)
    return result
end
function destructive_add!(ex::AbstractArray{<:GenericAffExpr}, c::Constant,
                          x::AbstractArray{<:GenericQuadExpr})
    result = c*x
    add_to_expression!.(result, ex)
    return result
end

# For some reason, the broadcast syntax `ex .+ c * x` fails if `x` is an
# `Adjoint`. But if we explicitly call `broadcast` it seems to work.
# See JuMP PR #1698 for more discussion.
destructive_add!(ex, c, x) = broadcast(+, ex, c * x)

_destructive_add_with_reorder!(ex, arg) = destructive_add!(ex, 1.0, arg)
# Special case because "Val{false}()" is used as the default empty expression.
_destructive_add_with_reorder!(ex::Val{false}, arg) = copy(arg)
# Calling `copy` on the matrix will not copy the entries
_destructive_add_with_reorder!(ex::Val{false}, arg::AbstractArray) = copy.(arg)
function _destructive_add_with_reorder!(ex::Val{false}, arg::Symmetric)
    Symmetric(copy.(arg))
end
_destructive_add_with_reorder!(ex::Val{false}, args...) = (*)(args...)


@generated function _destructive_add_with_reorder!(ex, x, y)
    if x <: Union{AbstractVariableRef,GenericAffExpr} && y <: Constant
        :(destructive_add!(ex, y, x))
    else
        :(destructive_add!(ex, x, y))
    end
end

@generated function _destructive_add_with_reorder!(ex, args...)
    n = length(args)
    @assert n ≥ 3
    varidx = findall(t -> (t <: AbstractVariableRef || t <: GenericAffExpr), collect(args))
    allscalar = all(t -> (t <: Constant), args[setdiff(1:n, varidx)])
    idx = (allscalar && length(varidx) == 1) ? varidx[1] : n
    coef = Expr(:call, :*, [:(args[$i]) for i in setdiff(1:n,idx)]...)
    :(destructive_add!(ex, $coef, args[$idx]))
end

# takes a generator statement and returns a properly nested for loop
# with nested filters as specified
function _parse_gen(ex, atleaf)
    if isexpr(ex, :flatten)
        return _parse_gen(ex.args[1], atleaf)
    end
    if !isexpr(ex, :generator)
        return atleaf(ex)
    end
    function itrsets(sets)
        if isa(sets, Expr)
            return sets
        elseif length(sets) == 1
            return sets[1]
        else
            return Expr(:block, sets...)
        end
    end

    idxvars = []
    if isexpr(ex.args[2], :filter) # if condition
        loop = Expr(:for, esc(itrsets(ex.args[2].args[2:end])),
                    Expr(:if, esc(ex.args[2].args[1]),
                          _parse_gen(ex.args[1], atleaf)))
        for idxset in ex.args[2].args[2:end]
            idxvar, s = _parse_idx_set(idxset)
            push!(idxvars, idxvar)
        end
    else
        loop = Expr(:for, esc(itrsets(ex.args[2:end])),
                         _parse_gen(ex.args[1], atleaf))
        for idxset in ex.args[2:end]
            idxvar, s = _parse_idx_set(idxset)
            push!(idxvars, idxvar)
        end
    end
    return loop
end

function _parse_generator(x::Expr, aff::Symbol, lcoeffs, rcoeffs, newaff=gensym())
    @assert isexpr(x,:call)
    @assert length(x.args) > 1
    @assert isexpr(x.args[2],:generator) || isexpr(x.args[2],:flatten)
    header = x.args[1]
    if _is_sum(header)
        _parse_generator_sum(x.args[2], aff, lcoeffs, rcoeffs, newaff)
    else
        error("Expected sum outside generator expression; got $header")
    end
end

function _parse_generator_sum(x::Expr, aff::Symbol, lcoeffs, rcoeffs, newaff)
    # We used to preallocate the expression at the lowest level of the loop.
    # When rewriting this some benchmarks revealed that it actually doesn't
    # seem to help anymore, so might as well keep the code simple.
    code = _parse_gen(x, t -> _parse_expr(t, aff, lcoeffs, rcoeffs, aff)[2])
    return :($code; $newaff=$aff)
end

_is_complex_expr(ex) = isa(ex,Expr) && !isexpr(ex,:ref)

_parse_expr_toplevel(x, aff::Symbol) = _parse_expr(x, aff, [], [])

function _is_comparison(ex::Expr)
    if isexpr(ex, :comparison)
        return true
    elseif isexpr(ex, :call)
        if ex.args[1] in (:<=, :≤, :>=, :≥, :(==))
            return true
        else
            return false
        end
    else
        return false
    end
end

# x[i=1] <= 2 is a somewhat common user error. Catch it here.
function _has_assignment_in_ref(ex::Expr)
    if isexpr(ex, :ref)
        return any(x -> isexpr(x, :(=)), ex.args)
    else
        return any(_has_assignment_in_ref, ex.args)
    end
end
_has_assignment_in_ref(other) = false

# output is assigned to newaff
function _parse_expr(x, aff::Symbol, lcoeffs::Vector, rcoeffs::Vector, newaff::Symbol=gensym())
    if !isa(x,Expr)
        # at the lowest level
        callexpr = Expr(:call, :_destructive_add_with_reorder!, aff,
                        lcoeffs..., esc(x), rcoeffs...)
        return newaff, :($newaff = $callexpr)
    else
        if x.head == :call && x.args[1] == :+
            b = Expr(:block)
            aff_ = aff
            for arg in x.args[2:(end-1)]
                aff_, code = _parse_expr(arg, aff_, lcoeffs, rcoeffs)
                push!(b.args, code)
            end
            newaff, code = _parse_expr(x.args[end], aff_, lcoeffs, rcoeffs, newaff)
            push!(b.args, code)
            return newaff, b
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                return _parse_expr(x.args[2], aff, vcat(-1.0, lcoeffs), rcoeffs, newaff)
            else # a - b - c ...
                b = Expr(:block)
                aff_, code = _parse_expr(x.args[2], aff, lcoeffs, rcoeffs)
                push!(b.args, code)
                for arg in x.args[3:(end-1)]
                    aff_,code = _parse_expr(arg, aff_, vcat(-1.0, lcoeffs), rcoeffs)
                    push!(b.args, code)
                end
                newaff,code = _parse_expr(x.args[end], aff_, vcat(-1.0, lcoeffs), rcoeffs, newaff)
                push!(b.args, code)
                return newaff, b
            end
        elseif x.head == :call && x.args[1] == :*
            # we might need to recurse on multiple arguments, e.g.,
            # (x+y)*(x+y)
            n_expr = mapreduce(_is_complex_expr, +, x.args)
            if n_expr == 1 # special case, only recurse on one argument and don't create temporary objects
                which_idx = 0
                for i in 2:length(x.args)
                    if _is_complex_expr(x.args[i])
                        which_idx = i
                    end
                end
                return _parse_expr(
                    x.args[which_idx], aff,
                    vcat(lcoeffs, [esc(x.args[i]) for i in 2:(which_idx - 1)]),
                    vcat(rcoeffs, [esc(x.args[i]) for i in (which_idx + 1):length(x.args)]),
                    newaff)
            else
                blk = Expr(:block)
                for i in 2:length(x.args)
                    if _is_complex_expr(x.args[i])
                        s = gensym()
                        newaff_, parsed = _parse_expr_toplevel(x.args[i], s)
                        push!(blk.args, :($s = 0.0; $parsed))
                        x.args[i] = newaff_
                    else
                        x.args[i] = esc(x.args[i])
                    end
                end
                callexpr = Expr(:call, :_destructive_add_with_reorder!, aff,
                                lcoeffs..., x.args[2:end]..., rcoeffs...)
                push!(blk.args, :($newaff = $callexpr))
                return newaff, blk
            end
        elseif x.head == :call && x.args[1] == :^ && _is_complex_expr(x.args[2])
            if x.args[3] == 2
                blk = Expr(:block)
                s = gensym()
                newaff_, parsed = _parse_expr_toplevel(x.args[2], s)
                push!(blk.args, :($s = Val(false); $parsed))
                push!(blk.args, :($newaff = _destructive_add_with_reorder!(
                    $aff, $(Expr(:call, :*, lcoeffs..., newaff_, newaff_,
                                 rcoeffs...)))))
                return newaff, blk
            elseif x.args[3] == 1
                return _parse_expr(:(JuMP.GenericQuadExpr($(x.args[2]))), aff, lcoeffs, rcoeffs)
            elseif x.args[3] == 0
                return _parse_expr(:(JuMP.GenericQuadExpr(one($(x.args[2])))), aff, lcoeffs, rcoeffs)
            else
                blk = Expr(:block)
                s = gensym()
                newaff_, parsed = _parse_expr_toplevel(x.args[2], s)
                push!(blk.args, :($s = Val(false); $parsed))
                push!(blk.args, :($newaff = _destructive_add_with_reorder!(
                    $aff, $(Expr(:call, :*, lcoeffs...,
                                 Expr(:call, :^, newaff_, esc(x.args[3])),
                                      rcoeffs...)))))
                return newaff, blk
            end
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            return _parse_expr(numerator, aff, lcoeffs,vcat(esc(:(1/$denom)),rcoeffs),newaff)
        elseif isexpr(x,:call) && length(x.args) >= 2 && (isexpr(x.args[2],:generator) || isexpr(x.args[2],:flatten))
            return newaff, _parse_generator(x,aff,lcoeffs,rcoeffs,newaff)
        elseif x.head == :curly
            _error_curly(x)
        else # at lowest level?
            if _is_comparison(x)
                error("Unexpected comparison in expression $x.")
            end
            if _has_assignment_in_ref(x)
                @warn "Unexpected assignment in expression $x. This will" *
                             " become a syntax error in a future release."
            end
            callexpr = Expr(:call, :_destructive_add_with_reorder!, aff,
                            lcoeffs..., esc(x), rcoeffs...)
            return newaff, :($newaff = $callexpr)
        end
    end
end
