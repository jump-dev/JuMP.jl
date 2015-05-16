#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# extracts left operands of *
# a*b*c -> a*b
function timescoef(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    x2 = copy(x)
    # just delete last argument
    splice!(x2.args,length(x2.args))
    # hacky change: if the call is e.g. (*)(a), just return :a
    if length(x2.args) == 2
        x2.args[2]
    else
        x2
    end
end
# extracts last operand of *
# a*b*c -> c
function timesvar(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    return x.args[end]
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::Number,x::V)
    push!(aff, convert(C,c), x)
    aff
end

function addToExpression(aff::GenericAffExpr,c::Number,x::Number)
    aff.constant += c*x
    aff
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::V,x::V)
    GenericQuadExpr{C,V}([c],[x],[one(C)],aff)
end

function addToExpression{C,V}(aff::GenericAffExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    append!(aff.vars, x.vars)
    append!(aff.coeffs, c*x.coeffs)
    aff.constant += c*x.constant
    aff
end

# TODO: add generic versions of following two methods
addToExpression{C}(aff::GenericAffExpr{C,Variable},c::GenericAffExpr{C,Variable},x::Variable) =
    GenericQuadExpr{C,Variable}(c.vars,
                                fill(x,length(c.vars)),
                                c.coeffs,
                                addToExpression(aff,c.constant,x))

addToExpression{C}(aff::GenericAffExpr{C,Variable},c::Variable,x::GenericAffExpr{C,Variable}) =
    GenericQuadExpr{C,Variable}(fill(c,length(x.vars)),
                                x.vars,
                                x.coeffs,
                                addToExpression(aff,x.constant,c))

addToExpression{C,V}(aff::GenericAffExpr{C,V}, c::Number, x::GenericQuadExpr{C,V}) =
    GenericQuadExpr{C,V}(copy(x.qvars1),
                         copy(x.qvars2),
                         c*x.qcoeffs,
                         addToExpression(aff,c,x.aff))

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::V)
    push!(quad.aff,convert(C,c),x)
    quad
end

function addToExpression(quad::GenericQuadExpr,c::Number,x::Number)
    quad.aff.constant += c*x
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::V,x::V)
    push!(quad.qvars1, x)
    push!(quad.qvars2, c)
    push!(quad.qcoeffs, 1.0)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericAffExpr{C,V})
    append!(quad.aff.vars, x.vars)
    append!(quad.aff.coeffs, c*x.coeffs)
    quad.aff.constant += c*x.constant
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::GenericAffExpr{C,V},x::Number)
    append!(quad.aff.vars, c.vars)
    append!(quad.aff.coeffs, c.coeffs*x)
    quad.aff.constant += c.constant*x
    quad
end

# TODO: add generic versions of following two methods
function addToExpression{C}(quad::GenericQuadExpr{C,Variable},c::GenericAffExpr{C,Variable},x::Variable)
    append!(quad.qvars1, c.vars)
    append!(quad.qvars2, fill(x,length(c.vars)))
    append!(quad.qcoeffs, c.coeffs)
    addToExpression(quad.aff,c.constant,x)
    quad
end

function addToExpression{C}(quad::GenericQuadExpr{C,Variable},c::Variable,x::GenericAffExpr{C,Variable})
    append!(quad.qvars1, fill(c,length(x.vars)))
    append!(quad.qvars2, x.vars)
    append!(quad.qcoeffs, x.coeffs)
    addToExpression(quad.aff,c,x.constant)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::GenericQuadExpr{C,V},x::Number)
    append!(quad.qvars1,c.qvars1)
    append!(quad.qvars2,c.qvars2)
    append!(quad.qcoeffs,c.qcoeffs*x)
    addToExpression(quad,c.aff,x)
    quad
end

function addToExpression{C,V}(quad::GenericQuadExpr{C,V},c::Number,x::GenericQuadExpr{C,V})
    append!(quad.qvars1,x.qvars1)
    append!(quad.qvars2,x.qvars2)
    append!(quad.qcoeffs,c*x.qcoeffs)
    addToExpression(quad,c,x.aff)
    quad
end

_lift(x::OneIndexedArray) = x.innerArray
_lift(x) = x

addToExpression(aff, c, x) = _lift(aff) + _lift(c) * _lift(x)

function parseCurly(x::Expr, aff::Symbol, constantCoef)
    if !(x.args[1] == :sum || x.args[1] == :∑ || x.args[1] == :Σ) # allow either N-ARY SUMMATION or GREEK CAPITAL LETTER SIGMA
        error("Expected sum outside curly braces")
    end
    if length(x.args) < 3
        error("Need at least two arguments for sum")
    end

    # we have a filter condition
    if isexpr(x.args[2],:parameters)
        cond = x.args[2]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        code = quote
            if $(esc(cond.args[1]))
                $(parseExpr(x.args[3], aff, constantCoef)[2])
            end
        end
        for level in length(x.args):-1:4
            code = :(
            for $(esc(x.args[level].args[1])) in $(esc(x.args[level].args[2]))
                $code
            end)
        end
    else # no condition
        code = parseExpr(x.args[2], aff, constantCoef)[2]
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
            if isa($aff,GenericAffExpr)
                sizehint!($aff.vars,length($aff.vars)+$len)
                sizehint!($aff.coeffs,length($aff.coeffs)+$len)
            else
                sizehint!($aff.qvars1,length($aff.qvars1)+$len)
                sizehint!($aff.qvars2,length($aff.qvars2)+$len)
                sizehint!($aff.qcoeffs,length($aff.qcoeffs)+$len)
                sizehint!($aff.aff.vars,length($aff.aff.vars)+$len)
                sizehint!($aff.aff.coeffs,length($aff.aff.coeffs)+$len)
            end
        end
        code = :($preblock;$code)
    end


    return code
end

parseExprToplevel(x, aff::Symbol) = parseExpr(x, aff, [1.0])
parseExpr(x, aff::Symbol, constantCoef::Vector) = parseExpr(x, aff, Expr(:call,:*,constantCoef...))

function parseExpr(x, aff::Symbol, constantCoef::Union(Number, Expr))
    if !isa(x,Expr)
        # at the lowest level
        return aff, :($aff = addToExpression($aff, $(esc(constantCoef)), $(esc(x))))
    else
        if x.head == :call && x.args[1] == :+
            return aff, Expr(:block,[parseExpr(arg,aff,constantCoef)[2] for arg in x.args[2:end]]...)
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                return parseExpr(x.args[2], aff, :(-1.0*$constantCoef))
            else # a - b - c ...
                return aff, Expr(:block,vcat(parseExpr(x.args[2], aff, constantCoef)[2],
                     {parseExpr(arg, aff, :(-1.0*$constantCoef))[2] for arg in x.args[3:end]})...)
            end
        elseif x.head == :call && x.args[1] == :*
            coef = timescoef(x)
            var = timesvar(x)
            return parseExpr(var, aff, :($constantCoef*$coef))
        elseif x.head == :call && x.args[1] == :/
            @assert length(x.args) == 3
            numerator = x.args[2]
            denom = x.args[3]
            return parseExpr(numerator, aff, :((1/$denom)*$constantCoef))
        elseif x.head == :curly
            return aff, parseCurly(x,aff,constantCoef)
        else # at lowest level?
            if isexpr(x,:comparison)
                error("Unexpected comparison in expression $x")
            end
            return aff, :($aff = addToExpression($aff, $(esc(constantCoef)), $(esc(x))))
        end
    end
end
