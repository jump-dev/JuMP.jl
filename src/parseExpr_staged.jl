#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# `a in b` is a comparison after JuliaLang/julia#13078
const in_is_compare = VERSION >= v"0.5.0-dev+901"

function tryParseIdxSet(arg::Expr)
    if arg.head === :(=) || (!in_is_compare && arg.head === :in)
        @assert length(arg.args) == 2
        return true, arg.args[1], arg.args[2]
    elseif (in_is_compare && arg.head === :comparison &&
            length(arg.args) == 3 && arg.args[2] === :in)
        return true, arg.args[1], arg.args[3]
    end
    return false, nothing, nothing
end

function parseIdxSet(arg::Expr)
    parse_done, idxvar, idxset = tryParseIdxSet(arg)
    if parse_done
        return idxvar, idxset
    end
    error("Invalid syntax: $arg")
end

# addToExpression(ex::Number, c::Number) = ex + c

addToExpression(ex::Number, c::Number, x::Number) = ex + c*x

addToExpression(ex::Number, c::Number, x::Variable) = AffExpr([x],[c],ex)

function addToExpression{T<:GenericAffExpr}(ex::Number, c::Number, x::T)
    # It's only safe to mutate the first argument.
    if c == 0
        T(ex)
    else
        x = copy(x)
        scale!(x.coeffs, c)
        x.constant *= c
        x.constant += ex
        x
    end
end

function addToExpression{T<:GenericQuadExpr}(ex::Number, c::Number, x::T)
    # It's only safe to mutate the first argument.
    if c == 0
        T(ex)
    else
        x = copy(x)
        scale!(x.qcoeffs, c)
        scale!(x.aff.coeffs, c)
        x.aff.constant *= c
        x.aff.constant += ex
        x
    end
end

addToExpression(ex::Number, c::Variable, x::Variable) = QuadExpr([c],[x],[1.0],zero(AffExpr))

function addToExpression(ex::Number, c::GenericAffExpr, x::GenericAffExpr)
    q = c*x
    q.aff.constant += ex
    q
end

function addToExpression{C,V}(ex::Number, c::GenericAffExpr{C,V}, x::V)
    q = c*x
    q.aff.constant += ex
    q
end

function addToExpression{T<:GenericQuadExpr}(ex::Number, c::T, x::Number)
    if x == 0
        T(ex)
    else
        q = c*x
        q.aff.constant += ex
        q
    end
end

function addToExpression(aff::GenericAffExpr, c::Number, x::Number)
    aff.constant += c*x
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::Number, x::V)
    if c != 0
        push!(aff.vars,   x)
        push!(aff.coeffs, c)
    end
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    if c != 0
        append!(aff.vars, x.vars)
        sizehint!(aff.coeffs, length(aff.coeffs)+length(x.coeffs))
        for i in 1:length(x.coeffs)
            push!(aff.coeffs, c*x.coeffs[i])
        end
        aff.constant += c*x.constant
    end
    aff
end

# help w/ ambiguity
addToExpression{C,V<:Number}(aff::GenericAffExpr{C,V}, c::Number, x::Number) = aff + c*x

function addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::V, x::Number)
    if x != 0
        push!(aff.vars,   c)
        push!(aff.coeffs, x)
    end
    aff
end

addToExpression{C,V}(aff::GenericAffExpr{C,V},c::V,x::V) =
    GenericQuadExpr{C,V}([c],[x],[one(C)],aff)

# TODO: add generic versions of following two methods
addToExpression(aff::AffExpr,c::AffExpr,x::Variable) =
    QuadExpr(c.vars,
             fill(x,length(c.vars)),
             c.coeffs,
             addToExpression(aff,c.constant,x))

addToExpression(aff::AffExpr,c::Variable,x::AffExpr) =
    QuadExpr(fill(c,length(x.vars)),
             x.vars,
             x.coeffs,
             addToExpression(aff,c,x.constant))

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::GenericAffExpr{C,V},x::Number)
    if x != 0
        append!(aff.vars, c.vars)
        append!(aff.coeffs, c.coeffs * x)
        aff.constant += c.constant * x
    end
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number)
    if x == 0
        GenericQuadExpr{C,V}(aff)
    else
        GenericQuadExpr{C,V}(copy(c.qvars1),
                            copy(c.qvars2),
                            c*qcoeffs*x,
                            addToExpression(aff,c.aff,x))
    end
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V})
    if c == 0
        GenericQuadExpr{C,V}(aff)
    else
        GenericQuadExpr{C,V}(copy(x.qvars1),
                            copy(x.qvars2),
                            c*x.qcoeffs,
                            addToExpression(aff,c,x.aff))
    end
end

function addToExpression{C,V}(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V})
    q = convert(GenericQuadExpr{C,V}, ex)
    addToExpression(q, one(C), c*x)
    q
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::V)
    if c != 0
        push!(quad.aff, convert(C,c), x)
    end
    quad
end

function addToExpression(quad::GenericQuadExpr,c::Number,x::Number)
    quad.aff.constant += c*x
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},x::V,y::V)
    push!(quad.qvars1, x)
    push!(quad.qvars2, y)
    push!(quad.qcoeffs, one(C))
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    if c != 0
        append!(quad.aff.vars, x.vars)
        sizehint!(quad.aff.coeffs, length(quad.aff.coeffs)+length(x.coeffs))
        for i in 1:length(x.coeffs)
            push!(quad.aff.coeffs, c*x.coeffs[i])
        end
        quad.aff.constant += c*x.constant
    end
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number)
    if x != 0
        addToExpression(quad.aff,c,x)
    end
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::V)
    append!(quad.qvars1, c.vars)
    append!(quad.qvars2, fill(x,length(c.vars)))
    append!(quad.qcoeffs, c.coeffs)
    addToExpression(quad.aff,c.constant,x)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::V,x::GenericAffExpr{C,V})
    append!(quad.qvars1, fill(c,length(x.vars)))
    append!(quad.qvars2, x.vars)
    append!(quad.qcoeffs, x.coeffs)
    addToExpression(quad.aff,c,x.constant)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number)
    if x != 0
        append!(quad.qvars1,c.qvars1)
        append!(quad.qvars2,c.qvars2)
        sizehint!(quad.qcoeffs, length(quad.qcoeffs)+length(c.qcoeffs))
        for i in 1:length(c.qcoeffs)
            push!(quad.qcoeffs, c.qcoeffs[i]*x)
        end
        addToExpression(quad,c.aff,x)
    end
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericQuadExpr{C,V})
    if c != 0
        append!(quad.qvars1,x.qvars1)
        append!(quad.qvars2,x.qvars2)
        sizehint!(quad.qcoeffs, length(quad.qcoeffs)+length(x.qcoeffs))
        for i in 1:length(x.qcoeffs)
            push!(quad.qcoeffs, c*x.qcoeffs[i])
        end
        addToExpression(quad,c,x.aff)
    end
    quad
end

function addToExpression{C,V}(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V})
    q = c*x
    addToExpression(ex, 1.0, q)
    ex
end

# Catch nonlinear expressions being used in addConstraint, etc.
typealias _NLExpr ReverseDiffSparse.ParametricExpression
_nlexprerr() = error("""Cannot use nonlinear expression in @addConstraint or @setObjective.
                        Use @addNLConstraint or @setNLObjective instead.""")
# Following three definitions avoid ambiguity warnings
addToExpression{C,V<:_NLExpr}(expr::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::V) = _nlexprerr()
addToExpression{C,V<:_NLExpr}(expr::GenericQuadExpr{C,V}, c::Number, x::V) = _nlexprerr()
addToExpression{C,V<:_NLExpr}(expr::GenericQuadExpr{C,V}, c::V, x::GenericAffExpr{C,V}) = _nlexprerr()
for T1 in (GenericAffExpr,GenericQuadExpr), T2 in (Number,Variable,GenericAffExpr,GenericQuadExpr)
    @eval addToExpression(::$T1, ::$T2, ::_NLExpr) = _nlexprerr()
    @eval addToExpression(::$T1, ::_NLExpr, ::$T2) = _nlexprerr()
end

addToExpression(ex, c, x) = ex + c*x

@generated addToExpression_reorder(ex, arg) = :(addToExpression(ex, 1.0, arg))

@generated function addToExpression_reorder(ex, x, y)
    if x <: Union{Variable,AffExpr} && y <: Number
        :(addToExpression(ex, y, x))
    else
        :(addToExpression(ex, x, y))
    end
end

@generated function addToExpression_reorder(ex, args...)
    n = length(args)
    @assert n ≥ 3
    varidx = find(t -> (t == Variable || t == AffExpr), collect(args))
    allscalar = all(t -> (t <: Number), args[setdiff(1:n, varidx)])
    idx = (allscalar && length(varidx) == 1) ? varidx[1] : n
    coef = Expr(:call, :*, [:(args[$i]) for i in setdiff(1:n,idx)]...)
    :(addToExpression(ex, $coef, args[$idx]))
end

function parseCurly(x::Expr, aff::Symbol, lcoeffs, rcoeffs, newaff=gensym())
    header = x.args[1]
    if length(x.args) < 3
        error("Need at least two arguments for $header")
    end
    if header ∈ [:sum, :∑, :Σ]
        parseSum(x, aff, lcoeffs, rcoeffs, newaff)
    elseif header ∈ [:norm1, :norm2, :norminf, :norm∞]
        parseNorm(header, x, aff, lcoeffs, rcoeffs, newaff)
    else
        error("Expected sum or norm2 outside curly braces; got $header")
    end
end

function parseSum(x::Expr, aff::Symbol, lcoeffs, rcoeffs, newaff)
    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        inneraff, innercode = parseExpr(x.args[3], aff, lcoeffs, rcoeffs, aff)
        code = quote
            if $(esc(cond.args[1]))
                $innercode
            end
        end
        for level in length(x.args):-1:4
            _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
            idxvar = esc(_idxvar)
            code = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $code
                end
            end)
        end
    else # no condition
        inneraff, code = parseExpr(x.args[2], aff, lcoeffs, rcoeffs, aff)
        for level in length(x.args):-1:3
            _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
            idxvar = esc(_idxvar)
            code = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $code
                end
            end)
        end
        len = :len
        # precompute the number of elements to add
        # this is unncessary if we're just summing constants
        _, lastidxset = parseIdxSet(x.args[length(x.args)]::Expr)
        preblock = :($len += length($(esc(lastidxset))))
        for level in (length(x.args)-1):-1:3
            _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
            idxvar = esc(_idxvar)
            preblock = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $preblock
                end
            end)
        end
        preblock = quote
            $len = 0
            $preblock
            _sizehint_expr!($aff,$len)
        end
        code = :($preblock;$code)
    end
    :($code; $newaff=$aff)
end

function parseNorm(normp::Symbol, x::Expr, aff::Symbol, lcoeffs, rcoeffs, newaff)
    @assert string(x.args[1])[1:4] == "norm"
    # we have a filter condition
    finalaff = gensym()
    gennorm  = gensym()
    len      = gensym()
    normexpr = gensym()
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        inneraff, innercode = parseExprToplevel(x.args[3], :normaff)
        code = quote
            if $(esc(cond.args[1]))
                normaff = 0.0
                $innercode
                push!($normexpr, $inneraff)
            end
        end
        for level in length(x.args):-1:4
            _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
            idxvar = esc(_idxvar)
            code = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $code
                end
            end)
        end
        preblock = :($normexpr = GenericAffExpr[])
    else # no condition
        inneraff, code = parseExprToplevel(x.args[2], :normaff)
        code = :(normaff = 0.0; $code; push!($normexpr, $inneraff))
        _, lastidxset = parseIdxSet(x.args[length(x.args)]::Expr)
        preblock = :($len += length($(esc(lastidxset))))
        for level in length(x.args):-1:3
            _idxvar, idxset = parseIdxSet(x.args[level]::Expr)
            idxvar = esc(_idxvar)
            code = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $code
                end
            end)
            preblock = :(let
                $(localvar(idxvar))
                for $idxvar in $(esc(idxset))
                    $preblock
                end
            end)
        end
        preblock = quote
            $len = 0
            $preblock
            $normexpr = GenericAffExpr[]
            _sizehint_expr!($normexpr, $len)
        end
    end
    if normp == :norm2
        param = 2
    elseif normp == :norm1
        param = 1
    elseif normp ∈ [:norminf, :norm∞]
        param = Inf
    else
        error("Unrecognized norm: $normp")
    end
    quote
        $preblock
        $code
        $gennorm = _build_norm($param,$normexpr)
        $newaff = $(Expr(:call, :addToExpression_reorder, aff,lcoeffs...,gennorm,rcoeffs...))
    end
end

is_complex_expr(ex) = isa(ex,Expr) && !isexpr(ex,:ref)

parseExprToplevel(x, aff::Symbol) = parseExpr(x, aff, [], [])

# output is assigned to newaff
function parseExpr(x, aff::Symbol, lcoeffs::Vector, rcoeffs::Vector, newaff::Symbol=gensym())
    if !isa(x,Expr)
        # at the lowest level
        callexpr = Expr(:call,:addToExpression_reorder,aff,lcoeffs...,esc(x),rcoeffs...)
        return newaff, :($newaff = $callexpr)
    else
        if x.head == :call && x.args[1] == :+
            b = Expr(:block)
            aff_ = aff
            for arg in x.args[2:(end-1)]
                aff_, code = parseExpr(arg, aff_, lcoeffs, rcoeffs)
                push!(b.args, code)
            end
            newaff, code = parseExpr(x.args[end], aff_, lcoeffs, rcoeffs, newaff)
            push!(b.args, code)
            return newaff, b
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                return parseExpr(x.args[2], aff, vcat(-1.0, lcoeffs), rcoeffs, newaff)
            else # a - b - c ...
                b = Expr(:block)
                aff_, code = parseExpr(x.args[2], aff, lcoeffs, rcoeffs)
                push!(b.args, code)
                for arg in x.args[3:(end-1)]
                    aff_,code = parseExpr(arg, aff_, vcat(-1.0, lcoeffs), rcoeffs)
                    push!(b.args, code)
                end
                newaff,code = parseExpr(x.args[end], aff_, vcat(-1.0, lcoeffs), rcoeffs, newaff)
                push!(b.args, code)
                return newaff, b
            end
        elseif x.head == :call && x.args[1] == :*
            # we might need to recurse on multiple arguments, e.g.,
            # (x+y)*(x+y)
            n_expr = mapreduce(is_complex_expr, +, x.args)
            if n_expr == 1 # special case, only recurse on one argument and don't create temporary objects
                which_idx = 0
                for i in 2:length(x.args)
                    if is_complex_expr(x.args[i])
                        which_idx = i
                    end
                end
                return parseExpr(x.args[which_idx], aff, vcat(lcoeffs, [esc(x.args[i]) for i in 2:(which_idx-1)]),
                                                         vcat(rcoeffs, [esc(x.args[i]) for i in (which_idx+1):length(x.args)]),
                                                         newaff)
            else
                blk = Expr(:block)
                for i in 2:length(x.args)
                    if is_complex_expr(x.args[i])
                        s = gensym()
                        newaff_, parsed = parseExprToplevel(x.args[i], s)
                        push!(blk.args, :($s = 0.0; $parsed))
                        x.args[i] = newaff_
                    else
                        x.args[i] = esc(x.args[i])
                    end
                end
                callexpr = Expr(:call,:addToExpression_reorder,aff,lcoeffs...,x.args[2:end]...,rcoeffs...)
                push!(blk.args, :($newaff = $callexpr))
                return newaff, blk
            end
        elseif x.head == :call && x.args[1] == :^ && is_complex_expr(x.args[2])
            x.args[3] == 2 || error("Only exponents of 2 are currently supported")
            blk = Expr(:block)
            s = gensym()
            newaff_, parsed = parseExprToplevel(x.args[2], s)
            push!(blk.args, :($s = 0.0; $parsed))
            push!(blk.args, :($newaff = $aff + $(Expr(:call,:*,lcoeffs...,newaff_,newaff_,rcoeffs...))))
            return newaff, blk
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            return parseExpr(numerator, aff, lcoeffs,vcat(esc(:(1/$denom)),rcoeffs),newaff)
        elseif x.head == :curly
            return newaff, parseCurly(x,aff,lcoeffs,rcoeffs,newaff)
        else # at lowest level?
            !isexpr(x,:comparison) || error("Unexpected comparison in expression $x")
            callexpr = Expr(:call,:addToExpression_reorder,aff,lcoeffs...,esc(x),rcoeffs...)
            return newaff, :($newaff = $callexpr)
        end
    end
end

# Semi-automatically generated precompilation hints
include("precompile.jl")
