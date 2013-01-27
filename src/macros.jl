

# extracts left operands of *
# a*b*c -> a*b
function timescoef(x::Expr)
    if x.head != :call || x.args[1] != :*
        error("In expression $x expected multiplication")
    end
    x2 = copy(x)
    # just delete last argument
    delete!(x2.args,length(x2.args))
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

# returns expression for array of the coefficient terms in the comprehension
function comprehensioncoef(x::Expr)
    if x.head != :comprehension
        error("In expression $x expected comprehension")
    end
    comp = expr(symbol("typed-comprehension"),{Float64,timescoef(x.args[1]),x.args[2:end]...})
    containers = Expr[ ex.args[2] for ex in x.args[2:end] ]
    out = quote l = 1 end
    for c in containers
        push!(out.args,:(l *= length($c)))
    end
    push!(out.args,:(reshape($comp,(l,))))
    return out
end

# returns expression for array of the variable terms in the comprehension
function comprehensionvar(x::Expr)
    if x.head != :comprehension
        error("In expression $x expected comprehension")
    end
    comp = expr(symbol("typed-comprehension"),{Variable,timesvar(x.args[1]),x.args[2:end]...})
    containers = Expr[ ex.args[2] for ex in x.args[2:end] ]
    out = quote l = 1 end
    for c in containers
        push!(out.args,:(l *= length($c)))
    end
    push!(out.args,:(reshape($comp,(l,))))
    return out
end

# parses top-level expression and returns expression for array of the coefficient terms
function topcoef(x::Expr)
    if x.head == :call && x.args[1] == :+
        return expr(:vcat,map(topcoef,x.args[2:end]))
    elseif x.head == :comprehension
        return comprehensioncoef(x)
    elseif x.head == :call && x.args[1] == :*
        println("timescoef of $x is $(timescoef(x))")
        return timescoef(x)
    else
        error("Unable to parse expression $x")
    end
end

function topvar(x::Expr)
    if x.head == :call && x.args[1] == :+
        return expr(:vcat,map(topvar,x.args[2:end]))
    elseif x.head == :comprehension
        return comprehensionvar(x)
    elseif x.head == :call && x.args[1] == :*
        return timesvar(x)
    else
        error("Unable to parse expression $x")
    end
end

function addToExpression(aff::AffExpr,c::Number,x::Variable)
    push!(aff.vars,x)
    push!(aff.coeffs,c)
end

function addToExpression(aff::AffExpr,c::Number,x::Number)
    aff.constant += c*x
end

function parseCurly(x::Expr, aff::Symbol, constantCoef)
    if (x.args[1] != :sum)
        error("Expected sum outside curly braces")
    end
    if length(x.args) < 3
        error("Need at least two arguments for sum")
    end

    # we have a filter condition
    if (isa(x.args[end],Expr) && x.args[end].head == :parameters)
        cond = x.args[end]
        if length(cond.args) != 1
            error("No commas after semicolon allowed in sum expression, use && for multiple conditions")
        end
        # generate inner loop code first and then wrap in for loops
        code = quote
            if $(cond.args[1])
                $(parseExpr(x.args[2], aff, constantCoef))
            end
        end
        for level in (length(x.args)-1):-1:3
            code = expr(:for, {x.args[level],expr(:block,{code})})
            # for $(x.args[level]) $code end
        end
    else # no condition
        code = parseExpr(x.args[2], aff, constantCoef)
        for level in length(x.args):-1:3
            code = expr(:for, {x.args[level],expr(:block,{code})})
            # for $(x.args[level]) $code end
        end
    end

    
    return code
end

function parseExpr(x, aff::Symbol, constantCoef)
    if !isa(x,Expr)
        # at the lowest level
        quote
            addToExpression($aff, $constantCoef, $x)
        end
    else
        if x.head == :call && x.args[1] == :+
            expr(:block,{parseExpr(arg,aff,constantCoef) for arg in x.args[2:end]})
        elseif x.head == :call && x.args[1] == :-
            if length(x.args) == 2 # unary subtraction
                parseExpr(x.args[2], aff, :(-1.0*$constantCoef))
            else # a - b - c ...
                expr(:block,vcat(parseExpr(x.args[2], aff, constantCoef),
                   {parseExpr(arg, aff, :(-1.0*$constantCoef)) for arg in x.args[3:end]}))
            end
        elseif x.head == :call && x.args[1] == :*
            coef = timescoef(x)
            var = timesvar(x)
            parseExpr(var, aff, :($coef*$constantCoef))
        elseif x.head == :curly
            parseCurly(x,aff,constantCoef)
        else # at lowest level?
            quote
                addToExpression($aff, $constantCoef, $x)
            end
        end
    end
end

macro addConstraint(m, x)
    if (x.head != :comparison)
        error("Expected comparison operator in constraint $x")
    end
    aff = gensym()
    lhs = :($(x.args[1]) - $(x.args[3])) # move everything to the lhs
    esc(quote
        $aff = AffExpr()
        $(parseExpr(lhs, aff, 1.0))
        addConstraint($m, Constraint($aff, $(string(x.args[2]))))
    end)
end
        

macro sumExpr(x)
    esc(:(AffExpr($(topvar(x)),convert(Vector{Float64},$(topcoef(x))),0.)))
end


