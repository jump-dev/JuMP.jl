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
    return x2
end
# extracts last operand of *
# a*b*c -> c
function timesvar(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    return x.args[end]
end

function addToExpression(aff::AffExpr,c::Number,x::Variable)
    push!(aff,convert(Float64,c),x)
    return aff
end

function addToExpression(aff::AffExpr,c::Number,x::Number)
    (aff.constant += c*x)
    return aff
end

function addToExpression(aff::AffExpr,c::Variable,x::Variable)
    q = QuadExpr([c],[x],[1.0],aff)
    return q
end

function addToExpression(aff::AffExpr,c::Number,x::AffExpr)
    append!(aff.vars, x.vars)
    append!(aff.coeffs, c*x.coeffs)
    aff.constant += c*x.constant
    return aff
end

addToExpression(aff::AffExpr,x::AffExpr,c::Variable) = addToExpression(aff,c,x)
addToExpression(aff::AffExpr,c::Variable,x::AffExpr) = QuadExpr(fill(c,length(x.vars)),
                                                                x.vars,
                                                                x.coeffs,
                                                                addToExpression(aff,x.constant,c))

addToExpression(aff::AffExpr, c::Number, x::QuadExpr) = QuadExpr(copy(x.qvars1),
                                                                 copy(x.qvars2),
                                                                 c*x.qcoeffs,
                                                                 addToExpression(aff,c,x.aff))

function addToExpression(quad::QuadExpr,c::Number,x::Variable)
    push!(quad.aff,convert(Float64,c),x)
    return quad
end

function addToExpression(quad::QuadExpr,c::Number,x::Number)
    (quad.aff.constant += c*x)
    return quad
end

function addToExpression(quad::QuadExpr,c::Variable,x::Variable)
    push!(quad.qvars1, x)
    push!(quad.qvars2, c)
    push!(quad.qcoeffs, 1.0)
    return quad
end

function addToExpression(quad::QuadExpr,c::Number,x::AffExpr)
    append!(quad.aff.vars, x.vars)
    append!(quad.aff.coeffs, c*x.coeffs)
    quad.aff.constant += c*x.constant
    return quad
end

addToExpression(quad::QuadExpr,x::AffExpr,c::Variable) = addToExpression(quad,c,x)
function addToExpression(quad::QuadExpr,c::Variable,x::AffExpr)
    append!(quad.qvars1, fill(c,length(x.vars)))
    append!(quad.qvars2, x.vars)
    append!(quad.qcoeffs, x.coeffs)
    addToExpression(quad.aff,x.constant,c)
    return quad
end

addToExpression(quad::QuadExpr,x::QuadExpr,c::Number) = addToExpression(quad,c,x)
function addToExpression(quad::QuadExpr,c::Number,x::QuadExpr)
    append!(quad.qvars1,x.qvars1)
    append!(quad.qvars2,x.qvars2)
    append!(quad.qcoeffs,c*x.qcoeffs)
    addToExpression(quad,c,x.aff)
    return quad
end

addToExpression(aff, c, x) = error("Cannot construct an affine expression with a term of type ($(typeof(c)))*($(typeof(x)))")

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
            return parseExpr(var, aff, :($coef*$constantCoef))
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
