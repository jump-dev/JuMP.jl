#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#addToExpression(ex::Number, c::Number) = ex + c

addToExpression(ex::Number, c::Number, x::Number) = ex + c*x

addToExpression(ex::Number, c::Number, x::Variable) = AffExpr([x],[c],ex)

function addToExpression(ex::Number, c::Number, x::GenericAffExpr)
    # It's only safe to mutate the first argument.
    x = copy(x)
    scale!(x.coeffs, c)
    x.constant *= c
    x.constant += ex
    x
end

function addToExpression(ex::Number, c::Number, x::GenericQuadExpr)
    # It's only safe to mutate the first argument.
    x = copy(x)
    scale!(x.qcoeffs, c)
    scale!(x.aff.coeffs, c)
    x.aff.constant *= c
    x.aff.constant += ex
    x
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

function addToExpression(ex::Number, c::GenericQuadExpr, x::Number)
    q = c*x
    q.aff.constant += ex
    q
end

function addToExpression(aff::GenericAffExpr, c::Number, x::Number)
    aff.constant += c*x
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::Number,x::V)
    push!(aff,convert(C,c),x)
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::V,x::V)
    GenericQuadExpr{C,V}([c],[x],[one(C)],aff)
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    append!(aff.vars, x.vars)
    sizehint!(aff.coeffs, length(aff.coeffs)+length(x.coeffs))
    for i in 1:length(x.coeffs)
        push!(aff.coeffs, c*x.coeffs[i])
    end
    aff.constant += c*x.constant
    aff
end

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

addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::GenericQuadExpr{C,V}, x::Number) =
    GenericQuadExpr{C,V}(copy(c.qvars1),
                         copy(c.qvars2),
                         c*qcoeffs*x,
                         addToExpression(aff,c.aff,x))

addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) =
    GenericQuadExpr{C,V}(copy(x.qvars1),
                         copy(x.qvars2),
                         c*x.qcoeffs,
                         addToExpression(aff,c,x.aff))

function addToExpression{C,V}(ex::GenericAffExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V})
    q = convert(GenericQuadExpr{C,V}, ex)
    addToExpression(q, one(C), c*x)
    q
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::V)
    push!(quad.aff, convert(C,c), x)
    quad
end

function addToExpression(quad::GenericQuadExpr,c::Number,x::Number)
    quad.aff.constant += c*x
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::V,x::V)
    push!(quad.qvars1, x)
    push!(quad.qvars2, c)
    push!(quad.qcoeffs, one(C))
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    append!(quad.aff.vars, x.vars)
    sizehint!(quad.aff.coeffs, length(quad.aff.coeffs)+length(x.coeffs))
    for i in 1:length(x.coeffs)
        push!(quad.aff.coeffs, c*x.coeffs[i])
    end
    quad.aff.constant += c*x.constant
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
    append!(quad.qvars1,c.qvars1)
    append!(quad.qvars2,c.qvars2)
    sizehint!(quad.qcoeffs, length(quad.qcoeffs)+length(c.qcoeffs))
    for i in 1:length(c.qcoeffs)
        push!(quad.qcoeffs, c.qcoeffs[i]*x)
    end
    addToExpression(quad,c.aff,x)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericQuadExpr{C,V})
    append!(quad.qvars1,x.qvars1)
    append!(quad.qvars2,x.qvars2)
    sizehint!(quad.qcoeffs, length(quad.qcoeffs)+length(x.qcoeffs))
    for i in 1:length(x.qcoeffs)
        push!(quad.qcoeffs, c*x.qcoeffs[i])
    end
    addToExpression(quad,c,x.aff)
    quad
end

function addToExpression{C,V}(ex::GenericQuadExpr{C,V}, c::GenericAffExpr{C,V}, x::GenericAffExpr{C,V})
    q = c*x
    addToExpression(ex, 1.0, q)
    ex
end

# Fallback to help in the simple case x::Array{AffExpr} + y::Array{AffExpr}
function addToExpression{N}(ex::Array{AffExpr,N}, c::Number, x::Array{AffExpr,N})
    size(ex) == size(x) || error("Incompatible sizes: $(size(ex)) + $(size(x))")
    for I in eachindex(ex)
        append!(ex[I], c * x[I])
    end
    ex
end

function addToExpression{N}(ex::Array{AffExpr,N}, c::Array{AffExpr,N}, x::Number)
    size(ex) == size(x) || error("Incompatible sizes: $(size(ex)) + $(size(x))")
    for I in eachindex(ex)
        append!(ex[I], c[I] * x)
    end
    ex
end

_lift(x::OneIndexedArray) = x.innerArray
_lift(x) = x

addToExpression(aff, c, x) = _lift(aff) + _lift(c) * _lift(x)

@generated addToExpression_reorder(ex, arg) = :(addToExpression(ex, 1.0, arg))

@generated function addToExpression_reorder(ex, args...)
    if !isleaftype(ex) || mapreduce(t -> !isleaftype(t), |, args)
        error("Can't process abstract types")
    end
    # how many of the multiplicands are variables?
    n_var = mapreduce(t -> (t == Variable || t == AffExpr), +, args)
    has_quad = mapreduce(t -> (t == QuadExpr), |, args)
    if n_var == 0 && !has_quad
        #println("No variables")
        return :(addToExpression(ex, 1.0, (*)(args...)))
    elseif n_var == 1 && !has_quad # linear
        #println("Linear")
        coef_expr = Expr(:call, :*)
        coef_idx = Int[]
        var_idx = 0
        for i in 1:length(args)
            if args[i] == Variable || args[i] == AffExpr
                var_idx = i
            else
                push!(coef_expr.args, :(args[$i]))
            end
        end
        return :(addToExpression(ex, $coef_expr, args[$var_idx]))
    else
        #println("Nonlinear")
        coef_expr = Expr(:call, :*, [:(args[$i]) for i in 1:(length(args)-1)]...)
        return :(addToExpression(ex, $coef_expr, args[$(length(args))]))
    end
end

function parseCurly(x::Expr, aff::Symbol, lcoeffs, rcoeffs)
    header = x.args[1]
    if length(x.args) < 3
        error("Need at least two arguments for $header")
    end
    if (header == :sum || header == :∑ || header == :Σ)
        parseSum(x, aff, lcoeffs, rcoeffs)
    elseif header == :norm2
        parseNorm(x, aff, lcoeffs, rcoeffs)
    else
        error("Expected sum or norm outside curly braces; got $header")
    end
end

function parseSum(x::Expr, aff::Symbol, lcoeffs, rcoeffs)
    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        newaff, innercode = parseExpr(x.args[3], aff, lcoeffs, rcoeffs, aff)
        @assert aff == newaff
        code = quote
            if $(esc(cond.args[1]))
                $innercode
            end
        end
        for level in length(x.args):-1:4
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
    else # no condition
        newaff, code = parseExpr(x.args[2], aff, lcoeffs, rcoeffs, aff)
        for level in length(x.args):-1:3
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
        len = :len
        # precompute the number of elements to add
        # this is unncessary if we're just summing constants
        preblock = :($len += length($(esc(x.args[length(x.args)].args[2]))))
        for level in (length(x.args)-1):-1:3
            preblock = Expr(:for, esc(x.args[level]),preblock)
        end
        preblock = quote
            $len = 0
            $preblock
            _sizehint_expr!($aff,$len)
        end
        code = :($preblock;$code)
    end
    code
end

function parseNorm(x::Expr, aff::Symbol, lcoeffs, rcoeffs)
    @assert string(x.args[1])[1:4] == "norm"
    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        newaff, innercode = parseExpr(x.args[3], :normaff, lcoeffs, rcoeffs)
        code = quote
            if $(esc(cond.args[1]))
                normaff = AffExpr()
                $innercode
                push!(normexpr, $newaff)
            end
        end
        for level in length(x.args):-1:4
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
        code = :(normexpr = AffExpr[]; $code; $aff = Norm{2}(normexpr))
    else # no condition
        newaff, code = parseExpr(x.args[2], :normaff, lcoeffs, rcoeffs)
        for level in length(x.args):-1:3
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                normaff = AffExpr();
                $code
                push!(normexpr, $newaff);
            end
            )
        end
        code = :(normexpr = AffExpr[]; $code; $aff = Norm{2}(normexpr))
        len = :len
        # precompute the number of elements to add
        # this is unncessary if we're just summing constants
        preblock = :($len += length($(esc(x.args[length(x.args)].args[2]))))
        for level in (length(x.args)-1):-1:3
            preblock = Expr(:for, esc(x.args[level]),preblock)
        end
        preblock = quote
            $len = 0
            $preblock
            _sizehint_expr!($aff, $len)
        end
        code = :($preblock;$code)
    end
    code
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
                        newaff_, parsed = parseExpr(x.args[i], s, [], [])
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
            newaff_, parsed = parseExpr(x.args[2], s, [], [])
            push!(blk.args, :($s = 0.0; $parsed))
            push!(blk.args, :($newaff = $newaff_*$newaff_))
            return newaff, blk
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            return parseExpr(numerator, aff, vcat(esc(:(1/$denom)),lcoeffs),rcoeffs,newaff)
        elseif x.head == :curly
            return newaff, :($(parseCurly(x,aff,lcoeffs,rcoeffs)); $newaff = $aff)
        else # at lowest level?
            !isexpr(x,:comparison) || error("Unexpected comparison in expression $x")
            callexpr = Expr(:call,:addToExpression_reorder,aff,lcoeffs...,esc(x),rcoeffs...)
            return newaff, :($newaff = $callexpr)
        end
    end
end
